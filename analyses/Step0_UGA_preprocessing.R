# =============================================================================
# Stage 0: Load, align, and pre-process Uganda rasters for the 
#           spillover-pandemic pipeline
#
# Inputs (3 raster files):
#   1. uga_pd_2020_1km.tif                        - Population density (~1km)
#   2. GLW4-2020.D-DA.CTL.tif                     - Cattle density (~10km)
#   3. 201501_Global_Travel_Time_to_Cities_UGA.tiff - Travel time to nearest 
#                                                     city >=50k pop (~1km)
#
# What this script does:
#   Step 1: Load and inspect all three rasters
#   Step 2: Align them to a common grid (the coarsest resolution, ~10km)
#   Step 3: Build spillover probability surface: pi(j) ~ L(j) * P(j)^gamma
#   Step 4: Identify urban hub cells from the travel time raster
#   Step 5: Sample spillover locations and extract 2D swaths
#   Step 6: Pre-compute gravity kernel weights for each swath
#   Step 7: Save everything and visualise
#
# Required packages: terra, sf
# =============================================================================

# --- Load packages -----------------------------------------------------------
# 'terra' is the modern replacement for 'raster' package, used for all 
# raster operations (reading, reprojecting, resampling, etc.)
library(terra)
library(ggplot2)

# --- Configuration -----------------------------------------------------------

# Path to the directory containing the three input raster files.
# UPDATE THIS to wherever you've saved them.
data_dir <- "data"

# File names
pop_file    <- file.path(data_dir, "uga_pd_2020_1km.tif")
cattle_file <- file.path(data_dir, "GLW4-2020.D-DA.CTL.tif")
travel_file <- file.path(data_dir, "201501_Global_Travel_Time_to_Cities_2015.tif")

# Output directory for processed data
out_dir <- "outputs/step0_UGA_preprocessing"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Gamma parameter for the spillover probability surface.
# pi(j) proportional to L(j) * P(j)^gamma
# gamma = 1 means simple product of cattle and population density.
# gamma < 1 would downweight very dense urban areas where per-capita 
# livestock contact is low. We keep it at 1 for now but it's explicit
# here so it can be varied later.
gamma <- 0.001

# Travel time threshold (minutes) for defining "hub" cells.
# Cells with travel time below this are considered to be at/near a city.
# 30 minutes means roughly within the urban footprint of a city >= 50k people.
hub_travel_time_threshold <- 30

# Number of spillover locations to pre-sample for the swath library.
# Each sampled location gets a 2D swath extracted. At simulation time,
# we draw from this library rather than querying rasters on the fly.
n_spillover_samples <- 500

# Half-width of the 2D swath (km). The swath extends this far on either
# side of the spillover-to-hub axis. 50km captures lateral diffusion
# and alternative corridors without requiring a full national grid.
swath_half_width_km <- 20

# Hub expansion radius (km) for swath extraction. After building the basic 
# transmission corridor from spillover to target hub cell, we expand the swath 
# to include all hub patch cells within this radius of the target. This gives 
# outbreaks access to broader urban connectivity once they reach any part of a 
# city - reflecting the reality that reaching the edge of Kampala provides 
# access to all of Kampala's transport networks, not just the narrow corridor 
# the outbreak traveled through. 50km captures most individual urban areas 
# while avoiding including distant cities.
hub_expansion_radius_km <- 50

# Random seed for reproducibility
set.seed(42)


# =============================================================================
# STEP 1: Load and inspect all three rasters
# =============================================================================
cat("\n========== STEP 1: Loading and inspecting rasters ==========\n\n")

# --- 1a. Population density --------------------------------------------------
# This is from WorldPop, giving estimated persons per ~1km grid cell (or per km²).
# Resolution is approximately 0.00833° (~1km at the equator).
cat("Loading population density...\n")
r_pop <- rast(pop_file)

# Print key properties
cat("  Dimensions:  ", nrow(r_pop), "x", ncol(r_pop), "\n")
cat("  Resolution:  ", res(r_pop), "degrees\n")
cat("  Resolution:  ~", round(res(r_pop)[1] * 111, 2), "km\n")
cat("  Extent:      ", as.vector(ext(r_pop)), "\n")
cat("  CRS:         ", crs(r_pop, describe = TRUE)$name, "\n")

# Quick summary of values (ignoring NAs)
pop_vals <- values(r_pop, na.rm = TRUE)
cat("  Value range: ", range(pop_vals), "\n")
cat("  Mean:        ", mean(pop_vals), "\n")
cat("  Valid cells: ", sum(!is.na(values(r_pop))), "/", ncell(r_pop), "\n\n")

# --- 1b. Cattle density ------------------------------------------------------
# From FAO GLW4 (Gridded Livestock of the World, version 4, 2020).
# Dasymetric (DA) version: cattle density has been statistically downscaled
# within census polygons using environmental covariates.
# Resolution is approximately 0.0833° (~10km at the equator).
# Units: number of cattle per km².
cat("Loading cattle density...\n")
r_cattle <- rast(cattle_file)

cat("  Dimensions:  ", nrow(r_cattle), "x", ncol(r_cattle), "\n")
cat("  Resolution:  ", res(r_cattle), "degrees\n")
cat("  Resolution:  ~", round(res(r_cattle)[1] * 111, 2), "km\n")
cat("  Extent:      ", as.vector(ext(r_cattle)), "\n")

cattle_vals <- values(r_cattle, na.rm = TRUE)
cat("  Value range: ", range(cattle_vals), "\n")
cat("  Mean:        ", mean(cattle_vals), "\n")
cat("  Valid cells: ", sum(!is.na(values(r_cattle))), "/", ncell(r_cattle), "\n\n")

# --- 1c. Travel time to nearest city -----------------------------------------
# From Weiss et al. (2018), updated to 2015. Gives travel time in MINUTES
# from each ~1km pixel to the nearest city with population >= 50,000.
# This serves as our measure of connectivity / remoteness.
cat("Loading travel time to cities...\n")
r_travel <- rast(travel_file)

# cat("  Dimensions:  ", nrow(r_travel), "x", ncol(r_travel), "\n")
# cat("  Resolution:  ", res(r_travel), "degrees\n")
# cat("  Resolution:  ~", round(res(r_travel)[1] * 111, 2), "km\n")
# cat("  Extent:      ", as.vector(ext(r_travel)), "\n")
# 
# travel_vals <- values(r_travel, na.rm = TRUE)
# cat("  Value range: ", range(travel_vals), " minutes\n")
# cat("  Mean:        ", mean(travel_vals), " minutes\n")
# cat("  Valid cells: ", sum(!is.na(values(r_travel))), "/", ncell(r_travel), "\n\n")


# =============================================================================
# STEP 2: Align all rasters to a common grid
# =============================================================================
cat("\n========== STEP 2: Aligning rasters to common grid ==========\n\n")

# Strategy: we resample everything to match the CATTLE raster's grid,
# because it is the coarsest (~10km). Resampling the 1km rasters down 
# to 10km is aggregation (averaging ~100 source pixels per target pixel),
# which is safe — we're just computing mean density within larger cells.
# Going the other way (upsampling cattle to 1km) would fabricate 
# sub-pixel detail that doesn't exist in the data.

# First, we need to crop all rasters to their common spatial extent
# (the intersection of all three bounding boxes).

# Find the intersection of all three extents
common_ext <- intersect(ext(r_pop), intersect(ext(r_cattle), ext(r_travel)))
cat("Common extent:\n")
cat("  xmin:", common_ext[1], " xmax:", common_ext[2], "\n")
cat("  ymin:", common_ext[3], " ymax:", common_ext[4], "\n\n")

# Crop the cattle raster to the common extent.
# This is our "template" — all other rasters will be resampled to match it.
r_cattle_crop <- crop(r_cattle, common_ext)
cat("Target grid (cattle resolution):\n")
cat("  Dimensions: ", nrow(r_cattle_crop), "x", ncol(r_cattle_crop), "\n")
cat("  Resolution: ~", round(res(r_cattle_crop)[1] * 111, 1), "km\n")
cat("  Total cells:", ncell(r_cattle_crop), "\n\n")

# Resample population density to the cattle grid.
# Method = "average" computes the mean of all source (1km) pixels that 
# fall within each target (~10km) cell. This gives us mean population 
# density (persons/km²) at the ~10km scale.
cat("Resampling population density to target grid (mean)...\n")
r_pop_aligned <- resample(r_pop, r_cattle_crop, method = "average")
cat("  Done. Dimensions: ", nrow(r_pop_aligned), "x", ncol(r_pop_aligned), "\n")

# Resample travel time to the cattle grid.
# Method = "average" gives the mean travel time within each ~10km cell.
# Alternative: "min" would give the best-case travel time (i.e., the 
# most connected point within each cell). Mean is more conservative.
cat("Resampling travel time to target grid (mean)...\n")
r_travel_aligned <- resample(r_travel, r_cattle_crop, method = "average")
cat("  Done. Dimensions: ", nrow(r_travel_aligned), "x", ncol(r_travel_aligned), "\n\n")

# The cattle raster is already on the target grid, but let's assign it
# a clean name for consistency.
r_cattle_aligned <- r_cattle_crop

# Build a land mask: cells where ALL THREE rasters have valid (non-NA) data
# and population is positive.
r_land_mask <- (!is.na(r_pop_aligned)) & 
  (!is.na(r_cattle_aligned)) & 
  (!is.na(r_travel_aligned)) &
  (r_pop_aligned > 0)

n_land <- sum(values(r_land_mask), na.rm = TRUE)
cat("Land mask: ", n_land, " cells with valid data in all rasters\n\n")

# Compute approximate totals as a sanity check.
# Cell area in km² (approximate, using the resolution in degrees * 111 km/deg).
cell_res_km <- res(r_cattle_crop)[1] * 111
cell_area_km2 <- cell_res_km^2
cat("Cell size: ~", round(cell_res_km, 1), "km, area: ~", 
    round(cell_area_km2, 0), "km²\n")

# Total population ~ sum of (density * area) across all land cells
total_pop <- sum(values(r_pop_aligned * r_land_mask) * cell_area_km2, na.rm = TRUE)
cat("Approximate total population: ", round(total_pop / 1e6, 1), "M\n")

# Total cattle ~ sum of (density * area) across all land cells
total_cattle <- sum(values(r_cattle_aligned * r_land_mask) * cell_area_km2, na.rm = TRUE)
cat("Approximate total cattle: ", round(total_cattle / 1e6, 1), "M\n\n")

###
# plot(r_cattle_aligned)
# plot(r_pop_aligned)
# plot(r_travel_aligned)

r_pop_aligned[!r_land_mask] <- NA
r_cattle_aligned[!r_land_mask] <- NA
r_travel_aligned[!r_land_mask] <- NA

# =============================================================================
# STEP 3: Build spillover probability surface
# =============================================================================
cat("\n========== STEP 3: Building spillover probability surface ==========\n\n")

# The probability that a spillover occurs in cell j is proportional to:
#
#   pi(j)  ~  L(j) * P(j)^gamma
#
# where:
#   L(j) = cattle density in cell j (animals per km²)
#   P(j) = human population density in cell j (persons per km²)
#   gamma = exponent on population density (default 1.0)
#
# Rationale: spillover requires contact between an infected animal and a 
# human. The rate of such contact events scales with both how many animals
# are present and how many humans are available to be exposed. The product
# naturally downweights cells that are high in one but zero in the other.

# Compute the unnormalised spillover weight for each cell.
# We clamp negative values to zero (shouldn't happen, but defensive).
r_spill_weight <- max(r_cattle_aligned, 0) * (max(r_pop_aligned, 0))^gamma

# Zero out cells outside the land mask
r_spill_weight <- r_spill_weight * r_land_mask

# Normalise to a probability distribution (sums to 1 over all cells)
total_weight <- sum(values(r_spill_weight), na.rm = TRUE)

if (total_weight <= 0) {
  stop("ERROR: Total spillover weight is zero. Check input rasters.")
}

r_spill_prob <- r_spill_weight / total_weight

cat("Spillover probability surface:\n")
cat("  Gamma:          ", gamma, "\n")
cat("  Non-zero cells: ", sum(values(r_spill_prob) > 0, na.rm = TRUE), "\n")
cat("  Sum:            ", sum(values(r_spill_prob), na.rm = TRUE), "\n")

# Show where the highest spillover probabilities are
spill_vals <- values(r_spill_prob)
top_quantiles <- quantile(spill_vals[spill_vals > 0], 
                          probs = c(0.5, 0.9, 0.95, 0.99), na.rm = TRUE)
cat("  Quantiles of non-zero values:\n")
print(top_quantiles)
cat("\n")


# =============================================================================
# STEP 4: Identify urban hub cells from travel time raster
# =============================================================================
cat("\n========== STEP 4: Identifying hub cells ==========\n\n")

# Hub cells are those with travel time below the threshold. These represent
# locations that are effectively "at" a city (within its urban footprint).
# In the pipeline, a spillover outbreak must reach a hub cell and establish
# N_crit infections there before being handed to PanDORA.

r_hub_mask <- r_land_mask & (r_travel_aligned < hub_travel_time_threshold)
n_hub_cells <- sum(values(r_hub_mask), na.rm = TRUE)
cat("Hub cells (travel time < ", hub_travel_time_threshold, " min): ", 
    n_hub_cells, "\n")

# To identify distinct hub clusters (separate cities), we use the patches()
# function which finds connected components of TRUE cells.
# directions = 8 means cells are connected if they touch diagonally too.
r_hub_patches <- patches(r_hub_mask, directions = 8, zeroAsNA = TRUE)

# Get the actual hub IDs (excluding NA/0)
actual_hub_ids <- unique(values(r_hub_patches))
actual_hub_ids <- actual_hub_ids[!is.na(actual_hub_ids) & actual_hub_ids > 0]
n_hubs <- length(actual_hub_ids)
cat("Distinct hub clusters: ", n_hubs, "\n\n")

# For each hub cluster, compute its centre and approximate population.
# We store this in a data frame for later use.
hub_info <- data.frame(
  hub_id   = integer(),
  n_cells  = integer(),
  pop_k    = numeric(),
  lat      = numeric(),
  lon      = numeric()
)

for (h in actual_hub_ids) {
  # Find which cells belong to this hub cluster
  hub_cells <- which(values(r_hub_patches) == h)
  
  # Get their row/col positions in the raster grid
  rc <- rowColFromCell(r_hub_patches, hub_cells)
  
  # Compute the geographic centre (mean of lat/lon of hub cells)
  coords <- xyFromCell(r_hub_patches, hub_cells)  # returns lon, lat columns
  center_lon <- mean(coords[, 1])
  center_lat <- mean(coords[, 2])
  
  # Approximate population of the hub
  hub_pop <- sum(values(r_pop_aligned)[hub_cells] * cell_area_km2, na.rm = TRUE)
  
  hub_info <- rbind(hub_info, data.frame(
    hub_id  = h,
    n_cells = length(hub_cells),
    pop_k   = round(hub_pop / 1e3),
    lat     = round(center_lat, 3),
    lon     = round(center_lon, 3)
  ))
  
  cat(sprintf("  Hub %d: %d cells, pop ~%dk, centre (%.2f, %.2f)\n",
              h, length(hub_cells), round(hub_pop / 1e3), center_lat, center_lon))
}

cat("\n")


# =============================================================================
# STEP 5: Sample spillover locations and extract 2D swaths
# =============================================================================
cat("\n========== STEP 5: Sampling spillovers & extracting swaths ==========\n\n")

# --- 5a. Sample spillover cell indices from the probability surface ----------

# Get the spillover probability for every cell as a flat vector
spill_probs_flat <- values(r_spill_prob)

# Replace NAs with 0 (non-land cells can't have spillovers)
spill_probs_flat[is.na(spill_probs_flat)] <- 0

# Sample cell indices, weighted by spillover probability.
# These are linear indices into the raster (1-based, as is R convention).
sampled_cells <- sample(
  x       = 1:ncell(r_spill_prob),
  size    = n_spillover_samples,
  replace = TRUE,
  prob    = spill_probs_flat
)

# Convert cell indices to geographic coordinates (lon, lat)
sampled_coords <- xyFromCell(r_spill_prob, sampled_cells)
sampled_lons <- sampled_coords[, 1]
sampled_lats <- sampled_coords[, 2]

# Look up the travel time at each sampled location
sampled_travel <- values(r_travel_aligned)[sampled_cells]

cat("Sampled", n_spillover_samples, "spillover locations\n")
cat("  Travel time to nearest city: median =", 
    round(median(sampled_travel, na.rm = TRUE)), "min,",
    "range = [", round(min(sampled_travel, na.rm = TRUE)), ",",
    round(max(sampled_travel, na.rm = TRUE)), "] min\n\n")

# --- 5b. For each spillover, determine routing and at-hub status -------------

# We need to know which hub each spillover should route toward, so we can 
# extract the swath connecting the spillover to that hub. We now route toward
# the nearest hub PATCH CELL (not hub center) for more realistic geography.

# Pre-compute a helper function: great-circle-ish distance in km
# (equirectangular approximation, fine for distances within Uganda)
dist_km <- function(lat1, lon1, lat2, lon2) {
  # lat/lon in degrees, returns distance in km
  dlat <- (lat2 - lat1) * 111.0
  dlon <- (lon2 - lon1) * 111.0 * cos(lat1 * pi / 180)
  sqrt(dlat^2 + dlon^2)
}

# Pre-compute all hub patch cells and their coordinates for efficient lookup
cat("Pre-computing hub patch cells for routing...\n")
all_hub_cells <- which(values(r_hub_patches) > 0)
all_hub_coords <- xyFromCell(r_hub_patches, all_hub_cells)
all_hub_patch_ids <- values(r_hub_patches)[all_hub_cells]

cat("  Total hub patch cells:", length(all_hub_cells), "\n")
cat("  Hub patches represented:", length(unique(all_hub_patch_ids)), "\n\n")

# For each sampled spillover, find the nearest hub PATCH CELL (for routing direction)
sampled_hub_id   <- integer(n_spillover_samples)
sampled_hub_lat  <- numeric(n_spillover_samples)
sampled_hub_lon  <- numeric(n_spillover_samples)
sampled_hub_dist <- numeric(n_spillover_samples)

cat("Finding nearest hub patch cell for each spillover...\n")
for (i in 1:n_spillover_samples) {
  # Compute distance from this spillover to every hub patch cell
  dists <- dist_km(sampled_lats[i], sampled_lons[i],
                   all_hub_coords[, 2], all_hub_coords[, 1])  # lat, lon
  
  # Pick the closest hub patch cell
  nearest_idx <- which.min(dists)
  nearest_cell <- all_hub_cells[nearest_idx]
  
  # Store the target information
  sampled_hub_id[i]   <- all_hub_patch_ids[nearest_idx]  # which patch this cell belongs to
  sampled_hub_lat[i]  <- all_hub_coords[nearest_idx, 2]   # lat of target cell
  sampled_hub_lon[i]  <- all_hub_coords[nearest_idx, 1]   # lon of target cell
  sampled_hub_dist[i] <- dists[nearest_idx]               # distance to target
  
  if (i %% 100 == 0) {
    cat(sprintf("  Processed %d / %d spillovers\n", i, n_spillover_samples))
  }
}

cat("\nDistance to nearest hub patch cell (for routing):\n")
cat("  Median:", round(median(sampled_hub_dist)), "km\n")
cat("  Range: [", round(min(sampled_hub_dist)), ",",
    round(max(sampled_hub_dist)), "] km\n")

# Determine which spillovers are "at-hub" based on hub PATCH areas
# If a spillover is already in any hub patch, there's no spatial bottleneck - 
# it can reach global transport immediately.
at_hub <- logical(n_spillover_samples)

for (i in 1:n_spillover_samples) {
  # Get the raster cell index for this spillover location
  spill_cell <- cellFromXY(r_hub_mask, cbind(sampled_lons[i], sampled_lats[i]))
  
  # Check if this cell is in any hub patch area (travel time < 30 min)
  if (!is.na(spill_cell) && spill_cell > 0 && spill_cell <= ncell(r_hub_mask)) {
    at_hub[i] <- values(r_hub_mask)[spill_cell]
    # Handle NA values in hub mask
    if (is.na(at_hub[i])) at_hub[i] <- FALSE
  } else {
    at_hub[i] <- FALSE
  }
}

cat("  Spillovers at-hub (in hub patch areas):", sum(at_hub), "\n")
cat("  Spillovers requiring spatial transmission:", sum(!at_hub), "\n")

# Show which hub patches are being targeted
target_patch_counts <- table(sampled_hub_id[!at_hub])
cat("  Target hub patches for spatial spillovers:\n")
for (patch_id in names(target_patch_counts)) {
  cat(sprintf("    Hub %s: %d spillovers\n", patch_id, target_patch_counts[patch_id]))
}
cat("\n")

# Visualize spillover locations, hub areas, and routing targets
cat("Plotting spillover locations and routing targets...\n")

plot(r_cattle_aligned, main = "Spillover Routing to Nearest Hub Patch Cells",
     col = hcl.colors(50, "Greens", rev = TRUE))

# Add hub patches as semi-transparent overlay
r_hub_outline <- r_hub_patches
r_hub_outline[r_hub_outline == 0 | is.na(r_hub_outline)] <- NA
plot(r_hub_outline, add = TRUE, 
     col = rainbow(length(actual_hub_ids), alpha = 0.4),
     legend = FALSE)

# Plot spillovers with different colors and routing lines
cell_width_deg <- res(r_cattle_aligned)[1]
jitter_amount <- cell_width_deg * 0.25
jittered_lons <- sampled_lons + runif(length(sampled_lons), -jitter_amount, jitter_amount)
jittered_lats <- sampled_lats + runif(length(sampled_lats), -jitter_amount, jitter_amount)

# At-hub spillovers in red (no spatial bottleneck)
if (sum(at_hub) > 0) {
  points(jittered_lons[at_hub], jittered_lats[at_hub], 
         pch = 16, col = "red", cex = 0.8)
}

# Spatial spillovers in black with routing lines to targets
if (sum(!at_hub) > 0) {
  # Plot spillover locations
  points(jittered_lons[!at_hub], jittered_lats[!at_hub], 
         pch = 16, col = "black", cex = 0.6)
  
  # Draw routing lines (sample a few to avoid clutter)
  n_lines_to_show <- min(20, sum(!at_hub))
  spatial_indices <- which(!at_hub)
  show_indices <- sample(spatial_indices, n_lines_to_show)
  
  for (idx in show_indices) {
    lines(c(sampled_lons[idx], sampled_hub_lon[idx]),
          c(sampled_lats[idx], sampled_hub_lat[idx]),
          col = "blue", lwd = 1, lty = 2)
  }
  
  # Mark target hub patch cells
  points(sampled_hub_lon[show_indices], sampled_hub_lat[show_indices],
         pch = 17, col = "blue", cex = 0.8)
}

legend("bottomright", 
       legend = c("Hub patches", "At-hub spillovers", "Spatial spillovers", 
                  "Routing lines", "Target cells"),
       col = c("purple", "red", "black", "blue", "blue"), 
       pch = c(15, 16, 16, NA, 17),
       lty = c(NA, NA, NA, 2, NA),
       cex = 0.8, bg = "white")

# --- 5c. Extract 2D swaths --------------------------------------------------

# For each spillover that is NOT at-hub, we extract a rectangular swath of 
# cells connecting the spillover to its nearest hub patch cell. The swath is oriented
# along the spillover-to-hub-target axis and extends swath_half_width_km on each side.
# We then EXPAND the swath to include nearby hub areas for realistic urban connectivity.

# Get all grid cell coordinates and data as vectors (for fast subsetting)
all_coords  <- xyFromCell(r_cattle_crop, 1:ncell(r_cattle_crop))
all_lons    <- all_coords[, 1]
all_lats    <- all_coords[, 2]
all_pop     <- values(r_pop_aligned)
all_cattle  <- values(r_cattle_aligned)
all_travel  <- values(r_travel_aligned)
all_land    <- values(r_land_mask)
# Replace NAs in the land mask with FALSE
all_land[is.na(all_land)] <- FALSE

cat("Extracting 2D swaths with hub expansion...\n")

# Hub expansion parameter for step 5c (swath extraction)
hub_expansion_radius_km <- 50  # include additional hub cells within this radius of target

# We'll store each swath as a list. The full collection goes in a list of lists.
swath_list <- vector("list", n_spillover_samples)

for (i in 1:n_spillover_samples) {
  
  # --- Handle at-hub spillovers ---
  # If the spillover is already at a hub, there's no spatial bottleneck.
  # We store a minimal record and the branching process will run non-spatially.
  if (at_hub[i]) {
    swath_list[[i]] <- list(
      at_hub        = TRUE,
      spill_lat     = sampled_lats[i],
      spill_lon     = sampled_lons[i],
      hub_id        = sampled_hub_id[i],
      hub_lat       = sampled_hub_lat[i],
      hub_lon       = sampled_hub_lon[i],
      travel_time   = sampled_travel[i],
      hub_dist_km   = sampled_hub_dist[i]
    )
    next
  }
  
  # --- Compute the swath geometry ---
  
  # Vector from spillover to TARGET HUB CELL, in km
  # (Note: now routing to nearest hub patch cell, not hub center)
  dx <- (sampled_hub_lon[i] - sampled_lons[i]) * 111 * 
    cos(sampled_lats[i] * pi / 180)
  dy <- (sampled_hub_lat[i] - sampled_lats[i]) * 111
  axis_length <- sqrt(dx^2 + dy^2)
  
  # Unit vector along the spillover-to-hub-target axis
  ux <- dx / axis_length
  uy <- dy / axis_length
  
  # Unit vector perpendicular to the axis (rotated 90°)
  px <- -uy
  py <- ux
  
  # For every grid cell, compute its displacement from the spillover in km
  cell_dx <- (all_lons - sampled_lons[i]) * 111 * 
    cos(sampled_lats[i] * pi / 180)
  cell_dy <- (all_lats - sampled_lats[i]) * 111
  
  # Project each cell's displacement onto the axis and perpendicular directions
  # "along" = how far along the spillover-to-hub-target axis (0 = at spillover,
  #           axis_length = at target hub cell)
  # "perp"  = how far off to the side of the axis
  along <- cell_dx * ux + cell_dy * uy
  perp  <- cell_dx * px + cell_dy * py
  
  # Select cells that fall within the basic swath:
  # - Along: from swath_half_width_km behind the spillover to swath_half_width_km past the target
  # - Perpendicular: within swath_half_width_km of the axis
  # - Must be on land (valid data in all rasters)
  in_swath <- which(
    all_land &
      (along >= -swath_half_width_km ) &
      (along <= axis_length + swath_half_width_km ) &
      (abs(perp) <= swath_half_width_km)
  )
  
  # --- Expand swath to include nearby hub patch areas ---
  # Once we have the basic transmission corridor, add nearby hub cells
  # to give the outbreak access to broader urban connectivity
  
  # Calculate distances from the TARGET hub cell to all other hub patch cells
  target_to_hub_dists <- dist_km(rep(sampled_hub_lat[i], nrow(all_hub_coords)),
                                 rep(sampled_hub_lon[i], nrow(all_hub_coords)),
                                 all_hub_coords[, 2], all_hub_coords[, 1])
  
  # Find hub cells within expansion radius of the target
  nearby_hub_cells <- all_hub_cells[target_to_hub_dists < hub_expansion_radius_km]
  
  # Find which nearby hub cells are NOT already included in the basic swath
  # Convert swath cell indices to cell numbers for comparison
  swath_cell_numbers <- cellFromXY(r_cattle_crop, cbind(all_lons[in_swath], all_lats[in_swath]))
  new_hub_cells <- setdiff(nearby_hub_cells, swath_cell_numbers)
  
  # Combine basic swath with expansion cells
  all_swath_cells <- c(in_swath, which(1:ncell(r_cattle_crop) %in% new_hub_cells))
  
  # If the swath has too few cells, skip (shouldn't happen in practice)
  if (length(all_swath_cells) < 5) {
    swath_list[[i]] <- list(
      at_hub    = TRUE,  # treat as at-hub if we can't build a swath
      spill_lat = sampled_lats[i],
      spill_lon = sampled_lons[i],
      hub_id    = sampled_hub_id[i],
      hub_lat   = sampled_hub_lat[i],
      hub_lon   = sampled_hub_lon[i],
      travel_time = sampled_travel[i],
      hub_dist_km = sampled_hub_dist[i],
      note      = "swath too small, treating as at-hub"
    )
    next
  }
  
  # Report expansion
  n_expansion_cells <- length(new_hub_cells)
  if (n_expansion_cells > 0 && i <= 10) {  # Report for first few swaths
    cat(sprintf("  Swath %d: added %d hub expansion cells (%.0f km radius)\n", 
                i, n_expansion_cells, hub_expansion_radius_km))
  }
  
  # --- Extract cell-level data for the expanded swath ---
  
  sw_lons   <- all_lons[all_swath_cells]
  sw_lats   <- all_lats[all_swath_cells]
  sw_pop    <- all_pop[all_swath_cells]         # population density (persons/km²)
  sw_cattle <- all_cattle[all_swath_cells]      # cattle density (cattle/km²)
  sw_travel <- all_travel[all_swath_cells]      # travel time to nearest city (min)
  sw_n      <- length(all_swath_cells)
  
  # Total population per cell (density * cell area)
  sw_pop_total <- sw_pop * cell_area_km2
  
  # Distance from each swath cell to the TARGET hub cell (km)
  sw_dist_to_hub <- dist_km(sw_lats, sw_lons,
                            sampled_hub_lat[i], sampled_hub_lon[i])
  
  # *** KEY FIX: Use hub PATCH areas, not distance to target ***
  # Boolean: which swath cells are "hub cells" based on travel time criterion
  # (same criterion used for at-hub detection: travel time < 30 min)
  sw_hub_mask <- logical(sw_n)
  for (j in 1:sw_n) {
    # Get the raster cell index for this swath cell
    cell_idx <- cellFromXY(r_hub_mask, cbind(sw_lons[j], sw_lats[j]))
    
    # Check if this cell is in any hub patch area
    if (!is.na(cell_idx) && cell_idx > 0 && cell_idx <= ncell(r_hub_mask)) {
      sw_hub_mask[j] <- values(r_hub_mask)[cell_idx]
      # Handle NA values in hub mask
      if (is.na(sw_hub_mask[j])) sw_hub_mask[j] <- FALSE
    } else {
      sw_hub_mask[j] <- FALSE
    }
  }
  
  # Find the index (within the swath) of the cell closest to the spillover.
  # This is where the branching process will start (I_{j0}(0) = 1).
  sw_dist_to_spill <- dist_km(sw_lats, sw_lons,
                              sampled_lats[i], sampled_lons[i])
  spill_cell_idx <- which.min(sw_dist_to_spill)
  
  # --- Compute pairwise distance matrix between all cells in the swath ---
  # This is needed for the gravity kernel in the branching process.
  # For a swath with N cells, this is an N x N matrix.
  # We use the equirectangular approximation (fast and fine for ~100km scales).
  
  # Coordinates in km (relative to an arbitrary origin; we use the spillover)
  sw_x_km <- (sw_lons - sampled_lons[i]) * 111 * 
    cos(sampled_lats[i] * pi / 180)
  sw_y_km <- (sw_lats - sampled_lats[i]) * 111
  
  # Pairwise distance matrix: dist_matrix[a, b] = distance in km between 
  # swath cell a and swath cell b.
  # We only compute this if the swath isn't too large (< 3000 cells).
  if (sw_n <= 3000) {
    dist_matrix <- as.matrix(dist(cbind(sw_x_km, sw_y_km)))
  } else {
    # For very large swaths, we'll compute distances on-the-fly during 
    # the branching process. Store coordinates instead.
    dist_matrix <- NULL
  }
  
  # --- Store the swath ---
  swath_list[[i]] <- list(
    at_hub          = FALSE,
    n_cells         = sw_n,
    n_expansion_cells = n_expansion_cells,  # Track how many cells added by expansion
    cell_lons       = sw_lons,
    cell_lats       = sw_lats,
    cell_pop_density = sw_pop,          # persons per km²
    cell_pop_total  = sw_pop_total,     # total persons in cell
    cell_cattle     = sw_cattle,        # cattle per km²
    cell_travel     = sw_travel,        # minutes to nearest city
    cell_dist_to_hub = sw_dist_to_hub,  # km to TARGET hub cell (for reference)
    hub_cell_mask   = sw_hub_mask,      # TRUE for cells in ANY hub patch area
    spill_cell_idx  = spill_cell_idx,   # index of spillover origin cell
    spill_lat       = sampled_lats[i],
    spill_lon       = sampled_lons[i],
    hub_id          = sampled_hub_id[i],
    hub_lat         = sampled_hub_lat[i],  # TARGET hub cell coordinates
    hub_lon         = sampled_hub_lon[i],  # TARGET hub cell coordinates
    length_km       = axis_length,      # straight-line distance to target
    travel_time     = sampled_travel[i],
    dist_matrix_km  = dist_matrix,      # N x N pairwise distances (or NULL)
    coords_km       = cbind(sw_x_km, sw_y_km)  # Nx2 coordinates in km
  )
  
  # Progress reporting
  if (i %% 100 == 0) {
    cat(sprintf("  Processed %d / %d spillovers\n", i, n_spillover_samples))
  }
}

# Summarise valid swaths with expansion info
valid_swaths <- Filter(function(s) !s$at_hub, swath_list)
at_hub_swaths <- Filter(function(s) s$at_hub, swath_list)

cat(sprintf("\nTotal swaths:   %d\n", length(swath_list)))
cat(sprintf("Valid (spatial): %d\n", length(valid_swaths)))
cat(sprintf("At-hub:          %d\n", length(at_hub_swaths)))

if (length(valid_swaths) > 0) {
  n_cells_vec <- sapply(valid_swaths, function(s) s$n_cells)
  length_vec  <- sapply(valid_swaths, function(s) s$length_km)
  hub_cells_vec <- sapply(valid_swaths, function(s) sum(s$hub_cell_mask))
  expansion_vec <- sapply(valid_swaths, function(s) s$n_expansion_cells)
  
  cat(sprintf("Cells per swath:      median = %d, range = [%d, %d]\n",
              median(n_cells_vec), min(n_cells_vec), max(n_cells_vec)))
  cat(sprintf("Distance to hub:      median = %.0f km, range = [%.0f, %.0f] km\n",
              median(length_vec), min(length_vec), max(length_vec)))
  cat(sprintf("Hub cells/swath:      median = %d, range = [%d, %d]\n",
              median(hub_cells_vec), min(hub_cells_vec), max(hub_cells_vec)))
  cat(sprintf("Expansion cells/swath: median = %d, range = [%d, %d]\n",
              median(expansion_vec), min(expansion_vec), max(expansion_vec)))
}

# Find the swath with the largest distance to hub target
max_distance <- 0
max_distance_idx <- 0

for (i in 1:length(swath_list)) {
  sw <- swath_list[[i]]

  # Skip at-hub spillovers (they have no length_km)
  if (!sw$at_hub && sw$length_km > max_distance) {
    max_distance <- sw$length_km
    max_distance_idx <- i
  }
}

cat(sprintf("\nSwath with largest distance:\n"))
cat(sprintf("  Index: %d\n", max_distance_idx))
cat(sprintf("  Distance: %.1f km\n\n", max_distance))

# Get the longest swath data
sw_longest <- swath_list[[max_distance_idx]]

cat("Longest swath details:\n")
cat("  Spillover:", round(sw_longest$spill_lat, 3), ",", round(sw_longest$spill_lon, 3), "\n")
cat("  Target:", round(sw_longest$hub_lat, 3), ",", round(sw_longest$hub_lon, 3), "\n")
cat("  Distance:", round(sw_longest$length_km, 1), "km\n")
cat("  Swath cells:", sw_longest$n_cells, "\n")
cat("  Hub cells in swath:", sum(sw_longest$hub_cell_mask), "\n")
cat("  Travel time:", round(sw_longest$travel_time, 0), "minutes\n\n")

# Plot the longest swath
plot(r_cattle_aligned, main = paste("Longest Swath (#", max_distance_idx, "): ",
                                    round(max_distance, 0), "km to Hub Patch #", sw_longest$hub_id),
     col = hcl.colors(50, "Greens", rev = TRUE))

# Highlight all swath cells in red
points(sw_longest$cell_lons, sw_longest$cell_lats, pch = 15, col = "red", cex = 0.3)

# Highlight hub cells within the swath in dark red (using the corrected hub_cell_mask)
hub_cells_in_swath <- sw_longest$hub_cell_mask
if (sum(hub_cells_in_swath) > 0) {
  points(sw_longest$cell_lons[hub_cells_in_swath], sw_longest$cell_lats[hub_cells_in_swath],
         pch = 15, col = "darkred", cex = 0.4)
}

# Mark spillover location with a star
points(sw_longest$spill_lon, sw_longest$spill_lat, pch = 8, col = "blue", cex = 2, lwd = 2)

# Mark TARGET hub cell with a triangle
points(sw_longest$hub_lon, sw_longest$hub_lat, pch = 17, col = "purple", cex = 2)

# Draw a line from spillover to target
lines(c(sw_longest$spill_lon, sw_longest$hub_lon),
      c(sw_longest$spill_lat, sw_longest$hub_lat),
      col = "black", lwd = 2, lty = 2)

# Add distance text along the line
mid_lon <- (sw_longest$spill_lon + sw_longest$hub_lon) / 2
mid_lat <- (sw_longest$spill_lat + sw_longest$hub_lat) / 2
text(mid_lon, mid_lat, paste(round(max_distance, 0), "km"),
     col = "white", font = 2, cex = 0.9,
     bg = "black", pos = 3)

# Add legend
legend("topright",
       legend = c("Spillover", "Target hub cell", "Swath cells", "Hub cells", "Axis"),
       col = c("blue", "purple", "red", "darkred", "black"),
       pch = c(8, 17, 15, 15, NA),
       lty = c(NA, NA, NA, NA, 2),
       cex = 0.8, bg = "white")

plot(r_hub_outline, add = TRUE,
     col = rainbow(length(actual_hub_ids), alpha = 0.4),
     legend = FALSE)

## Zooming in on the swath area
# Calculate zoom bounds around the longest swath
buffer <- 0.5  # degrees buffer around the swath (adjust as needed)

# Get the bounding box of the spillover and target
min_lon <- min(sw_longest$spill_lon, sw_longest$hub_lon) - buffer
max_lon <- max(sw_longest$spill_lon, sw_longest$hub_lon) + buffer  
min_lat <- min(sw_longest$spill_lat, sw_longest$hub_lat) - buffer
max_lat <- max(sw_longest$spill_lat, sw_longest$hub_lat) + buffer

# Optionally, expand to include all swath cells for even better context
if (length(sw_longest$cell_lons) > 0) {
  min_lon <- min(min_lon, min(sw_longest$cell_lons)) - buffer/2
  max_lon <- max(max_lon, max(sw_longest$cell_lons)) + buffer/2
  min_lat <- min(min_lat, min(sw_longest$cell_lats)) - buffer/2
  max_lat <- max(max_lat, max(sw_longest$cell_lats)) + buffer/2
}

# Plot with zoom
plot(r_cattle_aligned, main = paste("Longest Swath (#", max_distance_idx, ") - Zoomed: ", 
                                    round(max_distance, 0), "km to Hub Patch #", sw_longest$hub_id),
     col = hcl.colors(50, "Greens", rev = TRUE),
     xlim = c(min_lon, max_lon),
     ylim = c(min_lat, max_lat))

# Highlight all swath cells in red
points(sw_longest$cell_lons, sw_longest$cell_lats, pch = 15, col = "red", cex = 0.3)

# Highlight hub cells within the swath in dark red (using the corrected hub_cell_mask)
hub_cells_in_swath <- sw_longest$hub_cell_mask
if (sum(hub_cells_in_swath) > 0) {
  points(sw_longest$cell_lons[hub_cells_in_swath], sw_longest$cell_lats[hub_cells_in_swath], 
         pch = 15, col = "darkred", cex = 0.4)
}

# Mark spillover location with a star
points(sw_longest$spill_lon, sw_longest$spill_lat, pch = 8, col = "blue", cex = 2, lwd = 2)

# Mark TARGET hub cell with a triangle  
points(sw_longest$hub_lon, sw_longest$hub_lat, pch = 17, col = "purple", cex = 2)

# Draw a line from spillover to target
lines(c(sw_longest$spill_lon, sw_longest$hub_lon), 
      c(sw_longest$spill_lat, sw_longest$hub_lat), 
      col = "black", lwd = 2, lty = 2)

plot(r_hub_outline, add = TRUE, 
     col = rainbow(length(actual_hub_ids), alpha = 0.4),
     legend = FALSE)


# Add legend
legend("right", 
       legend = c("Spillover", "Target hub cell", "Swath cells", "Hub cells", "Axis"),
       col = c("blue", "purple", "red", "darkred", "black"),
       pch = c(8, 17, 15, 15, NA),
       lty = c(NA, NA, NA, NA, 2),
       cex = 0.8, bg = "white")

# Summarise
valid_swaths <- Filter(function(s) !s$at_hub, swath_list)
at_hub_swaths <- Filter(function(s) s$at_hub, swath_list)

cat(sprintf("\nTotal swaths:   %d\n", length(swath_list)))
cat(sprintf("Valid (spatial): %d\n", length(valid_swaths)))
cat(sprintf("At-hub:          %d\n", length(at_hub_swaths)))

if (length(valid_swaths) > 0) {
  n_cells_vec <- sapply(valid_swaths, function(s) s$n_cells)
  length_vec  <- sapply(valid_swaths, function(s) s$length_km)
  hub_cells_vec <- sapply(valid_swaths, function(s) sum(s$hub_cell_mask))
  
  cat(sprintf("Cells per swath:  median = %d, range = [%d, %d]\n",
              median(n_cells_vec), min(n_cells_vec), max(n_cells_vec)))
  cat(sprintf("Distance to hub:  median = %.0f km, range = [%.0f, %.0f] km\n",
              median(length_vec), min(length_vec), max(length_vec)))
  cat(sprintf("Hub cells/swath:  median = %d, range = [%d, %d]\n",
              median(hub_cells_vec), min(hub_cells_vec), max(hub_cells_vec)))
}


# =============================================================================
# STEP 6: Pre-compute gravity kernel weights for each swath
# =============================================================================
cat("\n\n========== STEP 6: Pre-computing gravity kernel weights ==========\n\n")

# The gravity kernel determines where each offspring case lands, given
# that the parent is in cell j. The probability that an offspring from 
# cell j lands in cell l is:
#
#   P(l | j) = w_{jl} / sum_m(w_{jm})
#
# where:
#   w_{jl} = N_l^alpha * exp(-beta * d_{jl})        for l != j
#   w_{jj} = N_j^alpha * exp(0) + w_self             for l == j
#          = N_j^alpha + w_self
#
# Parameters (with defaults — these will be varied/calibrated later):
alpha_gravity  <- 1.0     # population attractiveness exponent
beta_gravity   <- 0.1     # distance decay rate (per km). 1/beta = 10km decay length
w_self_gravity <- 100     # extra weight for within-cell transmission

cat(sprintf("Gravity kernel parameters:\n"))
cat(sprintf("  alpha  = %.2f  (population exponent)\n", alpha_gravity))
cat(sprintf("  beta   = %.3f  (decay rate, decay length = %.0f km)\n", 
            beta_gravity, 1/beta_gravity))
cat(sprintf("  w_self = %.0f  (self-transmission weight)\n\n", w_self_gravity))

# For each valid swath, compute the gravity weight matrix W[j, l] and 
# the normalised kernel P[j, l] = W[j, l] / sum_l(W[j, l]).
# We store P as a matrix since it's reused every generation of the 
# branching process (it doesn't change — susceptible depletion enters 
# through the effective R0, not the kernel).

for (i in seq_along(swath_list)) {
  sw <- swath_list[[i]]
  
  # Skip at-hub swaths
  
  if (sw$at_hub) next
  
  n <- sw$n_cells
  
  # If we have the distance matrix, compute kernel directly.
  # Otherwise skip (will compute on-the-fly during branching process).
  if (is.null(sw$dist_matrix_km)) {
    swath_list[[i]]$kernel_matrix <- NULL
    next
  }
  
  D <- sw$dist_matrix_km  # N x N distance matrix
  N_pop <- sw$cell_pop_total  # population (total persons) per cell
  
  # Compute the unnormalised weight matrix.
  # W[j, l] = N_pop[l]^alpha * exp(-beta * D[j, l])
  # We compute this as an outer product-like operation.
  
  # N_pop^alpha as a vector (one value per cell)
  pop_attract <- N_pop^alpha_gravity
  
  # Gravity weights: each row j gets pop_attract * exp(-beta * D[j,])
  # This is an N x N matrix where W[j, l] = pop_attract[l] * exp(-beta * D[j, l])
  W <- sweep(exp(-beta_gravity * D), 2, pop_attract, FUN = "*")
  
  # Add self-transmission weight to the diagonal
  diag(W) <- diag(W) + w_self_gravity
  
  # Normalise each row to get the probability kernel
  # P[j, l] = W[j, l] / sum(W[j, ])
  row_sums <- rowSums(W)
  P <- W / row_sums
  
  # Store the kernel matrix back in the swath
  swath_list[[i]]$kernel_matrix <- P
  
  if (i %% 100 == 0) {
    cat(sprintf("  Computed kernel for swath %d / %d (%d cells)\n", 
                i, n_spillover_samples, n))
  }
}

cat("Done computing gravity kernels.\n\n")

## Plotting to visualise and check the outputs given this kernel formulation
# Distance range to examine (0 to 50 km)
distances <- seq(0, 50, by = 0.5)

# Population scenarios to compare
pop_scenarios <- c(500, 2000, 5000, 10000)  # people per cell

# Create a data frame for plotting
kernel_data <- data.frame()

for (source_pop in pop_scenarios) {
  for (target_pop in pop_scenarios) {
    
    # Within-cell transmission (d = 0)
    w_self <- source_pop^alpha_gravity + w_self_gravity
    
    # Between-cell transmission for each distance
    w_between <- target_pop^alpha_gravity * exp(-beta_gravity * distances[-1])  # exclude d=0
    
    # Combine for normalization (approximate - assumes only these two cells exist)
    total_weight <- w_self + sum(w_between)
    
    # Transmission probabilities
    p_self <- w_self / total_weight
    p_between <- w_between / total_weight
    
    # Store results
    kernel_data <- rbind(kernel_data, data.frame(
      source_pop = source_pop,
      target_pop = target_pop,
      distance = 0,
      prob = p_self,
      type = "within-cell"
    ))
    
    kernel_data <- rbind(kernel_data, data.frame(
      source_pop = source_pop,
      target_pop = target_pop, 
      distance = distances[-1],
      prob = p_between,
      type = "between-cell"
    ))
  }
}

# Plot 1: Distance decay for different target populations (fixed source = 2000)
p1 <- ggplot(subset(kernel_data, source_pop == 2000 & type == "between-cell"), 
             aes(x = distance, y = prob, color = factor(target_pop))) +
  geom_line(size = 1) +
  # scale_y_log10() +
  labs(title = "Transmission Probability vs Distance\n(Source cell: 2000 people)",
       x = "Distance (km)", 
       y = "Transmission probability (log scale)",
       color = "Target\npopulation") +
  theme_minimal()

print(p1)

# Plot 2: Within-cell dominance across population sizes
within_cell_data <- subset(kernel_data, type == "within-cell" & target_pop == source_pop)
p2 <- ggplot(within_cell_data, aes(x = source_pop, y = prob)) +
  geom_line(size = 1, color = "red") +
  geom_point(size = 3, color = "red") +
  labs(title = "Within-Cell Transmission Probability\n(Effect of w_self = 500)",
       x = "Cell population", 
       y = "Probability of within-cell transmission") +
  theme_minimal()

print(p2)

# Plot 3: 2D visualization around a source cell
create_2d_kernel <- function(source_pop = 2000, grid_size = 21, cell_size_km = 5) {
  
  # Create a grid centered on the source
  center <- (grid_size + 1) / 2
  grid <- expand.grid(x = 1:grid_size, y = 1:grid_size)
  
  # Calculate distances from center
  grid$dist_km <- sqrt((grid$x - center)^2 + (grid$y - center)^2) * cell_size_km
  
  # Assume uniform target population
  target_pop <- 2000
  
  # Calculate weights
  grid$weight <- ifelse(grid$dist_km == 0, 
                        source_pop^alpha_gravity + w_self_gravity,  # within-cell
                        target_pop^alpha_gravity * exp(-beta_gravity * grid$dist_km))  # between-cell
  
  # Normalize to probabilities  
  grid$prob <- grid$weight / sum(grid$weight)
  
  return(grid)
}

grid_2d <- create_2d_kernel()

p3 <- ggplot(grid_2d, aes(x = x, y = y, fill = prob)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log10", name = "Transmission\nprobability") +
  geom_point(x = (21+1)/2, y = (21+1)/2, color = "red", size = 4, shape = 8) +  # source
  labs(title = "2D Gravity Kernel\n(Red star = source cell, 5km resolution)",
       x = "Grid X", y = "Grid Y") +
  theme_minimal() +
  coord_equal()

print(p3)

# Summary statistics
cat("\nGravity kernel summary:\n")
cat(sprintf("Distance decay length (1/beta): %.1f km\n", 1/beta_gravity))
cat(sprintf("w_self effect in small cell (1000 people): %.0f%% within-cell\n", 
            100 * (1000 + w_self_gravity) / (1000 + w_self_gravity + 1000)))
cat(sprintf("w_self effect in large cell (10000 people): %.0f%% within-cell\n", 
            100 * (10000 + w_self_gravity) / (10000 + w_self_gravity + 10000)))

# =============================================================================
# STEP 7: Save everything
# =============================================================================
cat("\n========== STEP 7: Saving outputs ==========\n\n")

# --- Save aligned rasters as GeoTIFFs ---
writeRaster(r_pop_aligned,    file.path(out_dir, "pop_aligned.tif"), overwrite = TRUE)
writeRaster(r_cattle_aligned, file.path(out_dir, "cattle_aligned.tif"), overwrite = TRUE)
writeRaster(r_travel_aligned, file.path(out_dir, "travel_aligned.tif"), overwrite = TRUE)
writeRaster(r_spill_prob,     file.path(out_dir, "spillover_prob.tif"), overwrite = TRUE)
writeRaster(r_hub_mask,       file.path(out_dir, "hub_mask.tif"), overwrite = TRUE)
writeRaster(r_hub_patches,    file.path(out_dir, "hub_patches.tif"), overwrite = TRUE)
cat("Saved aligned rasters to", out_dir, "\n")

# --- Save swath library as RDS (R's native serialisation format) ---
# This is the main output: a list of pre-computed swaths, each containing
# cell data, distance matrices, and gravity kernel matrices.
saveRDS(
  list(
    swaths   = swath_list,
    hub_info = hub_info,
    params   = list(
      gamma                   = gamma,
      hub_travel_time_threshold = hub_travel_time_threshold,
      n_spillover_samples     = n_spillover_samples,
      swath_half_width_km     = swath_half_width_km,
      cell_area_km2           = cell_area_km2,
      cell_res_km             = cell_res_km,
      alpha_gravity           = alpha_gravity,
      beta_gravity            = beta_gravity,
      w_self_gravity          = w_self_gravity
    )
  ),
  file = file.path(out_dir, "swath_library.rds")
)
cat("Saved swath library to", file.path(out_dir, "swath_library.rds"), "\n\n")


# =============================================================================
# STEP 8: Visualise
# =============================================================================
cat("\n========== STEP 8: Generating plots ==========\n\n")

png(file.path(out_dir, "raster_overview.png"), width = 1800, height = 1200, res = 150)
par(mfrow = c(2, 3), mar = c(3, 3, 3, 4))

# 1. Population density (log scale)
plot(log10(max(r_pop_aligned * r_land_mask, 0.1)), 
     main = "Population Density (log10)",
     col = hcl.colors(50, "YlOrRd", rev = TRUE))

# 2. Cattle density
plot(r_cattle_aligned * r_land_mask,
     main = "Cattle Density",
     col = hcl.colors(50, "Greens", rev = TRUE))

# 3. Travel time to cities
plot(r_travel_aligned * r_land_mask / 60,
     main = "Travel Time to Nearest City (hours)",
     col = hcl.colors(50, "RdYlGn"))

# 4. Spillover probability (log scale)
plot(log10(max(r_spill_prob, 1e-10)),
     main = sprintf("Spillover Probability (gamma=%.1f)", gamma),
     col = hcl.colors(50, "Reds", rev = TRUE))

# 5. Hub cells
plot(r_hub_patches, main = sprintf("Hub Clusters (n=%d)", n_hubs),
     col = hcl.colors(n_hubs, "Dark 2"))

# 6. Example swaths
plot(r_land_mask, col = c("white", "grey90"),
     main = sprintf("Example Swaths (n=%d)", min(20, length(valid_swaths))),
     legend = FALSE)
swath_cols <- rainbow(min(20, length(valid_swaths)), alpha = 0.5)
for (idx in seq_len(min(20, length(valid_swaths)))) {
  sw <- valid_swaths[[idx]]
  points(sw$cell_lons, sw$cell_lats, pch = ".", col = swath_cols[idx])
  points(sw$spill_lon, sw$spill_lat, pch = 8, col = swath_cols[idx], cex = 1.2)
  points(sw$hub_lon, sw$hub_lat, pch = 17, col = swath_cols[idx], cex = 1)
}

dev.off()
cat("Saved raster overview to", file.path(out_dir, "raster_overview.png"), "\n")

cat("\n✓ Stage 0 processing complete!\n")
cat("  All outputs in:", out_dir, "\n\n")