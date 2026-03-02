# =============================================================================
# Step5_UGA_actual_branchingProcess_run.R
#
# Runs Stages 3a, 3b, 3c of the spillover-pandemic pipeline on all 500
# Uganda spillover events using the pre-computed swath library (Step 0) and
# the viral property draws + route assignments (Step 1).
#
# For each event, dispatches to the appropriate branching process:
#   single_case    → no branching, chain_size = 1
#   nonspatial_sub → 1-cell subcritical (R₀ < 1), always goes extinct
#   spatial        → full multi-cell swath with gravity kernel
#   nonspatial_hub → 1-cell hub, R₀ ≥ 1, can establish
#
# Inputs:
#   - outputs/step0_UGA_preprocessing/swath_library.rds  (from Step 0)
#   - outputs/step1_viral_draws/event_table.rds          (from Step 1)
#   - functions/branching_process.R                       (simulation engine)
#
# Outputs (saved to outputs/step5_branching_results/):
#   - results_table.rds        : full results data frame (one row per event)
#   - results_table.csv        : same, as CSV for inspection
#   - hub_established_details.rds : detailed gen_log for hub-established events
#   - step5_params.rds         : parameter values used
#   - diagnostic plots (outcome breakdown, chain sizes, spatial outcomes, etc.)
#
# Usage:
#   source("analyses/Step5_UGA_actual_branchingProcess_run.R")
#   # or: Rscript analyses/Step5_UGA_actual_branchingProcess_run.R
# =============================================================================

library(ggplot2)

# --- Source function definitions ---------------------------------------------
source("functions/branching_process.R")

# --- Configuration -----------------------------------------------------------

# Input files
swath_file <- "outputs/step0_UGA_preprocessing/swath_library.rds"
event_file <- "outputs/step1_viral_draws/event_table.rds"

# Output directory
out_dir <- "outputs/step5_branching_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Branching process parameters
N_crit  <- 500L     # cumulative hub infections to trigger PanDORA handoff
g_max   <- 100L     # maximum generations before forced termination

# Log detailed generation history for spatial and hub events with R₀ ≥ 1?
# (We always log gen_log for hub-established events; this controls whether
#  we also log it for all spatial/hub events regardless of outcome.)
log_all_spatial <- FALSE   # set TRUE for debugging; FALSE for production speed

# Reproducibility
random_seed <- 123

# IFR for mortality tallying (Stage 5 of the pipeline)
# Used here for completeness — can be overridden downstream
IFR <- 0.01

# =============================================================================
# LOAD INPUTS
# =============================================================================
cat("\n================================================================\n")
cat("  Step 5: Branching Process Simulation (Uganda)\n")
cat("================================================================\n\n")

# --- Load swath library (Step 0) ---
if (!file.exists(swath_file)) {
  stop("Swath library not found at: ", swath_file,
       "\n  Run Step 0 first (Step0_UGA_preprocessing.R)")
}
cat("Loading swath library:", swath_file, "\n")
swath_library <- readRDS(swath_file)
swaths   <- swath_library$swaths
hub_info <- swath_library$hub_info
params_s0 <- swath_library$params

n_swaths <- length(swaths)
cat(sprintf("  %d swaths loaded\n", n_swaths))
cat(sprintf("  Hub info: %d hubs\n", nrow(hub_info)))
cat(sprintf("  Gravity kernel params: alpha=%.1f, beta=%.3f, w_self=%.0f\n",
            params_s0$alpha_gravity, params_s0$beta_gravity, params_s0$w_self_gravity))
cat(sprintf("  N_crit = %d, g_max = %d\n\n", N_crit, g_max))

# --- Load event table (Step 1) ---
if (!file.exists(event_file)) {
  stop("Event table not found at: ", event_file,
       "\n  Run Step 1 first (Step1_UGA_branchingProcess.R)")
}
cat("Loading event table:", event_file, "\n")
event_table <- readRDS(event_file)
n_events <- nrow(event_table)

cat(sprintf("  %d events loaded\n", n_events))
cat(sprintf("  Route breakdown:\n"))
route_tab <- table(event_table$route)
for (r in names(route_tab)) {
  cat(sprintf("    %-18s  %4d  (%5.1f%%)\n", r, route_tab[r],
              100 * route_tab[r] / n_events))
}
cat("\n")

# Sanity check: event count matches swath count
if (n_events != n_swaths) {
  stop(sprintf("Mismatch: %d events but %d swaths. Check inputs.", n_events, n_swaths))
}


# =============================================================================
# PRE-FLIGHT CHECKS
# =============================================================================
cat("--- Pre-flight checks ---\n")

# Check that spatial events have valid kernel matrices
n_spatial <- sum(event_table$route == "spatial")
n_kernel_ok <- 0
n_kernel_missing <- 0

for (i in which(event_table$route == "spatial")) {
  sw <- swaths[[i]]
  if (!is.null(sw$kernel_matrix) && !is.null(sw$cell_pop_total)) {
    n_kernel_ok <- n_kernel_ok + 1
  } else {
    n_kernel_missing <- n_kernel_missing + 1
    cat(sprintf("  WARNING: Event %d (spatial) has no kernel matrix\n", i))
  }
}
cat(sprintf("  Spatial events: %d total, %d with valid kernels, %d missing\n",
            n_spatial, n_kernel_ok, n_kernel_missing))

# Check hub info coverage for nonspatial_hub events
n_hub_events <- sum(event_table$route == "nonspatial_hub")
hub_ids_needed <- unique(event_table$hub_id[event_table$route == "nonspatial_hub"])
hub_ids_available <- hub_info$hub_id
missing_hubs <- setdiff(hub_ids_needed, hub_ids_available)
if (length(missing_hubs) > 0) {
  cat(sprintf("  WARNING: %d hub IDs needed but missing from hub_info: %s\n",
              length(missing_hubs), paste(missing_hubs, collapse = ", ")))
} else {
  cat(sprintf("  Hub events: %d total, all hub IDs found in hub_info\n", n_hub_events))
}

cat("  Pre-flight checks passed.\n\n")


# =============================================================================
# RUN BRANCHING PROCESSES
# =============================================================================
cat("================================================================\n")
cat("  Running branching processes for all events...\n")
cat("================================================================\n\n")

set.seed(random_seed)

# Pre-allocate results storage
results_list <- vector("list", n_events)

# Store detailed gen_logs for hub-established events (for downstream analysis)
hub_established_logs <- list()

# Timing
t_start <- Sys.time()

# Counters for progress reporting
n_done <- 0
n_single <- 0
n_sub    <- 0
n_spat   <- 0
n_hub    <- 0

for (i in 1:n_events) {
  
  # --- Build event list (simulate_event expects a named list) ---
  ev <- list(
    event_id = event_table$event_id[i],
    R0       = event_table$R0[i],
    k        = event_table$k[i],
    route    = event_table$route[i]
  )
  
  sw <- swaths[[i]]
  
  # --- Determine logging for this event ---
  # Always log gen_log for spatial and hub events (relatively cheap;
  # invaluable for diagnostics). Skip cell_log in production (expensive
  # for large swaths).
  do_log_gens  <- (ev$route %in% c("spatial", "nonspatial_hub")) || log_all_spatial
  do_log_cells <- FALSE  # set TRUE only for targeted debugging
  
  # --- Run simulation ---
  res <- simulate_event(
    event     = ev,
    swath     = sw,
    hub_info  = hub_info,
    N_crit    = N_crit,
    g_max     = g_max,
    log_gens  = do_log_gens,
    log_cells = do_log_cells
  )
  
  # --- Store result ---
  results_list[[i]] <- list(
    event_id       = res$event_id,
    route          = res$route,
    R0             = res$R0,
    k              = res$k,
    outcome        = res$outcome,
    chain_size     = res$chain_size,
    hub_infections = res$hub_infections,
    generations    = res$generations,
    peak_infected  = res$peak_infected
  )
  
  # If hub was established, save the detailed gen_log for later analysis
  if (res$outcome == "hub_established" && !is.null(res$gen_log)) {
    hub_established_logs[[as.character(i)]] <- list(
      event_id = res$event_id,
      R0       = res$R0,
      k        = res$k,
      route    = res$route,
      gen_log  = res$gen_log
    )
  }
  
  # --- Update counters ---
  n_done <- n_done + 1
  if (ev$route == "single_case")    n_single <- n_single + 1
  if (ev$route == "nonspatial_sub") n_sub    <- n_sub + 1
  if (ev$route == "spatial")        n_spat   <- n_spat + 1
  if (ev$route == "nonspatial_hub") n_hub    <- n_hub + 1
  
  # --- Progress reporting (every 50 events, or at milestones) ---
  if (i %% 50 == 0 || i == n_events) {
    elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    rate <- i / elapsed
    eta <- (n_events - i) / rate
    
    # Count outcomes so far
    outcomes_so_far <- sapply(results_list[1:i], function(r) r$outcome)
    n_extinct_sf  <- sum(outcomes_so_far == "extinct")
    n_hub_est_sf  <- sum(outcomes_so_far == "hub_established")
    n_gen_cap_sf  <- sum(outcomes_so_far == "generation_cap")
    n_single_sf   <- sum(outcomes_so_far == "single_case")
    
    cat(sprintf("  [%4d / %d]  %.1f sec  (%.0f events/sec, ETA %.0f sec)\n",
                i, n_events, elapsed, rate, eta))
    cat(sprintf("    Outcomes: extinct=%d, hub_est=%d, gen_cap=%d, single=%d\n",
                n_extinct_sf, n_hub_est_sf, n_gen_cap_sf, n_single_sf))
  }
}

t_end <- Sys.time()
elapsed_total <- as.numeric(difftime(t_end, t_start, units = "secs"))
cat(sprintf("\nAll %d events processed in %.1f seconds (%.0f events/sec)\n\n",
            n_events, elapsed_total, n_events / elapsed_total))


# =============================================================================
# BUILD RESULTS DATA FRAME
# =============================================================================
cat("Building results data frame...\n")

results_df <- data.frame(
  event_id       = sapply(results_list, function(r) r$event_id),
  route          = sapply(results_list, function(r) r$route),
  R0             = sapply(results_list, function(r) r$R0),
  k              = sapply(results_list, function(r) r$k),
  outcome        = sapply(results_list, function(r) r$outcome),
  chain_size     = sapply(results_list, function(r) r$chain_size),
  hub_infections = sapply(results_list, function(r) r$hub_infections),
  generations    = sapply(results_list, function(r) r$generations),
  peak_infected  = sapply(results_list, function(r) r$peak_infected),
  stringsAsFactors = FALSE
)

# --- Merge spatial metadata from event_table ---
results_df$at_hub       <- event_table$at_hub
results_df$spill_lat    <- event_table$spill_lat
results_df$spill_lon    <- event_table$spill_lon
results_df$hub_id       <- event_table$hub_id
results_df$hub_dist_km  <- event_table$hub_dist_km
results_df$travel_time  <- event_table$travel_time
results_df$n_swath_cells <- event_table$n_swath_cells
results_df$n_hub_cells  <- event_table$n_hub_cells

# --- Compute mortality (Stage 5 of the pipeline) ---
# For sub-pandemic events: M = chain_size * IFR
# For hub-established events: M = NA (would be determined by PanDORA)
results_df$mortality_sub_pandemic <- ifelse(
  results_df$outcome %in% c("single_case", "extinct", "generation_cap"),
  results_df$chain_size * IFR,
  NA_real_
)
results_df$pandora_handoff <- results_df$outcome == "hub_established"

cat(sprintf("  %d rows, %d columns\n\n", nrow(results_df), ncol(results_df)))


# =============================================================================
# SUMMARY STATISTICS
# =============================================================================
cat("================================================================\n")
cat("  RESULTS SUMMARY\n")
cat("================================================================\n\n")

# --- Overall outcome breakdown ---
cat("--- Overall outcome breakdown ---\n")
outcome_tab <- table(results_df$outcome)
for (o in names(outcome_tab)) {
  cat(sprintf("  %-20s  %4d  (%5.1f%%)\n", o, outcome_tab[o],
              100 * outcome_tab[o] / n_events))
}
cat("\n")

# --- Outcome by route ---
cat("--- Outcome breakdown by route ---\n")
cross_tab <- table(results_df$route, results_df$outcome)
print(cross_tab)
cat("\n")

# --- Hub establishment details ---
n_hub_est <- sum(results_df$outcome == "hub_established")
cat(sprintf("--- Hub establishment: %d events (%.2f%% of all, %.2f%% of R₀≥1) ---\n",
            n_hub_est,
            100 * n_hub_est / n_events,
            ifelse(sum(results_df$R0 >= 1) > 0,
                   100 * n_hub_est / sum(results_df$R0 >= 1), 0)))

if (n_hub_est > 0) {
  hub_est <- results_df[results_df$outcome == "hub_established", ]
  cat(sprintf("  By route:\n"))
  hub_route_tab <- table(hub_est$route)
  for (r in names(hub_route_tab)) {
    cat(sprintf("    %-18s  %d\n", r, hub_route_tab[r]))
  }
  cat(sprintf("  R₀:  median = %.2f, range = [%.2f, %.2f]\n",
              median(hub_est$R0), min(hub_est$R0), max(hub_est$R0)))
  cat(sprintf("  Chain size at handoff: median = %d, range = [%d, %d]\n",
              median(hub_est$chain_size), min(hub_est$chain_size), max(hub_est$chain_size)))
  cat(sprintf("  Generations to establishment: median = %d, range = [%d, %d]\n",
              median(hub_est$generations), min(hub_est$generations), max(hub_est$generations)))
  
  # Spatial hub establishments specifically
  spatial_hub_est <- hub_est[hub_est$route == "spatial", ]
  if (nrow(spatial_hub_est) > 0) {
    cat(sprintf("\n  Spatial hub establishments (%d events):\n", nrow(spatial_hub_est)))
    cat(sprintf("    R₀:           median = %.2f, range = [%.2f, %.2f]\n",
                median(spatial_hub_est$R0), min(spatial_hub_est$R0), max(spatial_hub_est$R0)))
    cat(sprintf("    Hub dist (km): median = %.0f, range = [%.0f, %.0f]\n",
                median(spatial_hub_est$hub_dist_km, na.rm = TRUE),
                min(spatial_hub_est$hub_dist_km, na.rm = TRUE),
                max(spatial_hub_est$hub_dist_km, na.rm = TRUE)))
    cat(sprintf("    Generations:   median = %d, range = [%d, %d]\n",
                median(spatial_hub_est$generations), min(spatial_hub_est$generations),
                max(spatial_hub_est$generations)))
  }
}
cat("\n")

# --- Chain size statistics by outcome ---
cat("--- Chain size statistics ---\n")
for (o in c("single_case", "extinct", "generation_cap", "hub_established")) {
  subset_o <- results_df[results_df$outcome == o, ]
  if (nrow(subset_o) > 0) {
    cat(sprintf("  %-20s  n=%4d  chain: median=%6d  mean=%8.1f  max=%6d\n",
                o, nrow(subset_o),
                median(subset_o$chain_size), mean(subset_o$chain_size),
                max(subset_o$chain_size)))
  }
}
cat("\n")

# --- Spatial bottleneck effectiveness ---
n_spatial_total <- sum(results_df$route == "spatial")
n_spatial_est   <- sum(results_df$route == "spatial" & results_df$outcome == "hub_established")
n_spatial_ext   <- sum(results_df$route == "spatial" & results_df$outcome == "extinct")
n_spatial_cap   <- sum(results_df$route == "spatial" & results_df$outcome == "generation_cap")

cat("--- Spatial bottleneck effectiveness ---\n")
cat(sprintf("  Spatial events (R₀ ≥ 1, not at hub): %d\n", n_spatial_total))
if (n_spatial_total > 0) {
  cat(sprintf("    → Extinct:          %4d  (%5.1f%%)\n",
              n_spatial_ext, 100 * n_spatial_ext / n_spatial_total))
  cat(sprintf("    → Hub established:  %4d  (%5.1f%%)\n",
              n_spatial_est, 100 * n_spatial_est / n_spatial_total))
  cat(sprintf("    → Generation cap:   %4d  (%5.1f%%)\n",
              n_spatial_cap, 100 * n_spatial_cap / n_spatial_total))
  cat(sprintf("  Spatial bottleneck blocks %.1f%% of supercritical rural spillovers\n",
              100 * (1 - n_spatial_est / n_spatial_total)))
}
cat("\n")

# --- Pipeline funnel (overall) ---
cat("--- Pipeline funnel ---\n")
cat(sprintf("  500 spillover events\n"))
cat(sprintf("   → %d dead-end (R₀ = 0)           [%.1f%%]\n",
            sum(results_df$R0 == 0), 100 * mean(results_df$R0 == 0)))
cat(sprintf("   → %d single case (R₀ < R_min)     [%.1f%%]\n",
            sum(results_df$outcome == "single_case"),
            100 * mean(results_df$outcome == "single_case")))
cat(sprintf("   → %d subcritical chains (R₀ < 1)  [%.1f%%]\n",
            sum(results_df$outcome == "extinct" & results_df$route == "nonspatial_sub"),
            100 * mean(results_df$outcome == "extinct" & results_df$route == "nonspatial_sub")))
cat(sprintf("   → %d supercritical (R₀ ≥ 1)       [%.1f%%]\n",
            sum(results_df$R0 >= 1), 100 * mean(results_df$R0 >= 1)))
cat(sprintf("     → %d blocked by spatial bottleneck\n",
            n_spatial_ext + n_spatial_cap))
cat(sprintf("     → %d reach hub establishment → PanDORA handoff\n", n_hub_est))
cat(sprintf("  Overall P(PanDORA handoff | spillover) = %.4f  (1 in %.0f)\n",
            n_hub_est / n_events,
            ifelse(n_hub_est > 0, n_events / n_hub_est, Inf)))
cat("\n")


# =============================================================================
# SAVE OUTPUTS
# =============================================================================
cat("--- Saving outputs to:", out_dir, "---\n")

# Results table
saveRDS(results_df, file.path(out_dir, "results_table.rds"))
cat("  results_table.rds\n")

write.csv(results_df, file.path(out_dir, "results_table.csv"), row.names = FALSE)
cat("  results_table.csv\n")

# Hub-established event details (gen_logs)
if (length(hub_established_logs) > 0) {
  saveRDS(hub_established_logs, file.path(out_dir, "hub_established_details.rds"))
  cat(sprintf("  hub_established_details.rds (%d events)\n", length(hub_established_logs)))
}

# Parameter record
step5_params <- list(
  N_crit        = N_crit,
  g_max         = g_max,
  IFR           = IFR,
  random_seed   = random_seed,
  n_events      = n_events,
  swath_file    = swath_file,
  event_file    = event_file,
  elapsed_secs  = elapsed_total,
  timestamp     = Sys.time()
)
saveRDS(step5_params, file.path(out_dir, "step5_params.rds"))
cat("  step5_params.rds\n\n")


# =============================================================================
# DIAGNOSTIC PLOTS
# =============================================================================
cat("--- Generating diagnostic plots ---\n\n")

# Colour palette
outcome_colours <- c(
  "single_case"     = "#999999",
  "extinct"         = "#56B4E9",
  "generation_cap"  = "#F0E442",
  "hub_established" = "#D55E00"
)

route_colours <- c(
  "single_case"    = "#999999",
  "nonspatial_sub" = "#56B4E9",
  "spatial"        = "#E69F00",
  "nonspatial_hub" = "#D55E00"
)

outcome_labels <- c(
  "single_case"     = "Single case",
  "extinct"         = "Extinct",
  "generation_cap"  = "Generation cap",
  "hub_established" = "Hub established (→ PanDORA)"
)

# --- Plot 1: Overall outcome breakdown ---
outcome_counts <- as.data.frame(table(results_df$outcome))
names(outcome_counts) <- c("outcome", "count")
outcome_counts$pct <- 100 * outcome_counts$count / sum(outcome_counts$count)
outcome_counts$outcome <- factor(outcome_counts$outcome,
                                 levels = c("single_case", "extinct",
                                            "generation_cap", "hub_established"))

p1 <- ggplot(outcome_counts, aes(x = outcome, y = count, fill = outcome)) +
  geom_col(width = 0.7, alpha = 0.9) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", count, pct)),
            vjust = -0.4, size = 4, fontface = "bold") +
  scale_fill_manual(values = outcome_colours, labels = outcome_labels) +
  scale_x_discrete(labels = outcome_labels) +
  labs(title = "Overall Outcome Distribution",
       subtitle = sprintf("%d spillover events | N_crit = %d | g_max = %d",
                          n_events, N_crit, g_max),
       x = NULL, y = "Number of events") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 15, hjust = 1),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, max(outcome_counts$count) * 1.15))

print(p1)

# --- Plot 2: Outcome by route (stacked bar) ---
cross_df <- as.data.frame(table(results_df$route, results_df$outcome))
names(cross_df) <- c("route", "outcome", "count")
cross_df$route <- factor(cross_df$route,
                         levels = c("single_case", "nonspatial_sub",
                                    "spatial", "nonspatial_hub"))
cross_df$outcome <- factor(cross_df$outcome,
                           levels = c("single_case", "extinct",
                                      "generation_cap", "hub_established"))

p2 <- ggplot(cross_df, aes(x = route, y = count, fill = outcome)) +
  geom_col(width = 0.7, alpha = 0.9) +
  scale_fill_manual(values = outcome_colours, labels = outcome_labels, name = "Outcome") +
  labs(title = "Outcome by Route",
       subtitle = "Each bar = one route category; colours = outcome",
       x = "Route", y = "Number of events") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 15, hjust = 1),
        legend.position = "bottom")

print(p2)

# --- Plot 3: Chain size distribution (extinct events only, log scale) ---
extinct_df <- results_df[results_df$outcome == "extinct", ]
if (nrow(extinct_df) > 0) {
  
  p3 <- ggplot(extinct_df, aes(x = chain_size, fill = route)) +
    geom_histogram(bins = 50, colour = "white", alpha = 0.85) +
    scale_x_log10() +
    scale_fill_manual(values = route_colours, name = "Route") +
    labs(title = "Chain Size Distribution (Extinct Events Only)",
         subtitle = sprintf("n = %d events | median chain = %d",
                            nrow(extinct_df), median(extinct_df$chain_size)),
         x = "Total chain size (log scale)", y = "Count") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  
  print(p3)
}

# --- Plot 4: R₀ vs outcome (for events with R₀ > 0) ---
df_r0 <- results_df[results_df$R0 > 0, ]
df_r0$outcome <- factor(df_r0$outcome,
                        levels = c("single_case", "extinct",
                                   "generation_cap", "hub_established"))

if (nrow(df_r0) > 0) {
  p4 <- ggplot(df_r0, aes(x = R0, fill = outcome)) +
    geom_histogram(bins = 50, colour = "white", alpha = 0.85) +
    scale_x_log10(breaks = c(0.01, 0.05, 0.1, 0.3, 1, 3, 10),
                  labels = c("0.01", "0.05", "0.1", "0.3", "1", "3", "10")) +
    scale_fill_manual(values = outcome_colours, labels = outcome_labels, name = "Outcome") +
    geom_vline(xintercept = 1.0, colour = "red", linewidth = 0.8, linetype = "dashed") +
    labs(title = "R₀ Distribution Coloured by Outcome (non-zero R₀ only)",
         subtitle = sprintf("Red dashed line = R₀ = 1 | %d events shown", nrow(df_r0)),
         x = expression(R[0] ~ "(log scale)"), y = "Count") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  
  print(p4)
}

# --- Plot 5: Spatial events — hub distance vs outcome ---
spatial_df <- results_df[results_df$route == "spatial", ]
if (nrow(spatial_df) > 0) {
  
  p5 <- ggplot(spatial_df, aes(x = hub_dist_km, fill = outcome)) +
    geom_histogram(bins = 30, colour = "white", alpha = 0.85) +
    scale_fill_manual(values = outcome_colours, labels = outcome_labels, name = "Outcome") +
    labs(title = "Spatial Events: Distance to Hub by Outcome",
         subtitle = sprintf("n = %d spatial events (R₀ ≥ 1, not at hub)", nrow(spatial_df)),
         x = "Distance to nearest hub (km)", y = "Count") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  
  print(p5)
}

# --- Plot 6: Spatial events — R₀ vs hub distance, coloured by outcome ---
if (nrow(spatial_df) > 0) {
  
  p6 <- ggplot(spatial_df, aes(x = R0, y = hub_dist_km, colour = outcome)) +
    geom_point(alpha = 0.7, size = 2.5) +
    scale_colour_manual(values = outcome_colours, labels = outcome_labels, name = "Outcome") +
    labs(title = "Spatial Events: R₀ vs Hub Distance",
         subtitle = "Each point is one spatial spillover event",
         x = expression(R[0]),
         y = "Distance to nearest hub (km)") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  
  print(p6)
}

# --- Plot 7: Generations to outcome (spatial + hub events) ---
supercrit_df <- results_df[results_df$R0 >= 1 & results_df$outcome != "single_case", ]
if (nrow(supercrit_df) > 0) {
  
  p7 <- ggplot(supercrit_df, aes(x = generations, fill = outcome)) +
    geom_histogram(bins = 40, colour = "white", alpha = 0.85) +
    scale_fill_manual(values = outcome_colours, labels = outcome_labels, name = "Outcome") +
    labs(title = "Generations to Outcome (Supercritical Events, R₀ ≥ 1)",
         subtitle = sprintf("n = %d events", nrow(supercrit_df)),
         x = "Number of generations", y = "Count") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  
  print(p7)
}

# --- Plot 8: Pipeline funnel (waterfall-style) ---
funnel_data <- data.frame(
  stage = c(
    "All spillovers",
    "Has R₀ > 0",
    "Has R₀ ≥ R_min",
    "Has R₀ ≥ 1",
    "Hub established"
  ),
  count = c(
    n_events,
    sum(results_df$R0 > 0),
    sum(results_df$R0 >= 0.05),
    sum(results_df$R0 >= 1),
    n_hub_est
  )
)
funnel_data$stage <- factor(funnel_data$stage, levels = funnel_data$stage)
funnel_data$pct <- 100 * funnel_data$count / n_events

p8 <- ggplot(funnel_data, aes(x = stage, y = count)) +
  geom_col(fill = "steelblue", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", count, pct)),
            vjust = -0.3, size = 3.5, fontface = "bold") +
  labs(title = "Pipeline Funnel: Spillover to Pandemic Handoff",
       subtitle = sprintf("%d events → %d PanDORA handoffs (%.2f%%)",
                          n_events, n_hub_est, 100 * n_hub_est / n_events),
       x = NULL, y = "Number of events") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 15, hjust = 1)) +
  coord_cartesian(ylim = c(0, n_events * 1.15))

print(p8)


# =============================================================================
# FINAL REPORT
# =============================================================================
cat("\n================================================================\n")
cat("  Step 5 complete\n")
cat("================================================================\n")
cat(sprintf("  Events processed:  %d\n", n_events))
cat(sprintf("  Runtime:           %.1f seconds\n", elapsed_total))
cat(sprintf("  Output directory:  %s\n", out_dir))
cat(sprintf("\n  Outcome summary:\n"))
for (o in c("single_case", "extinct", "generation_cap", "hub_established")) {
  n_o <- sum(results_df$outcome == o)
  cat(sprintf("    %-20s  %4d  (%5.1f%%)\n", o, n_o, 100 * n_o / n_events))
}
cat(sprintf("\n  PanDORA handoffs:  %d  (%.4f per spillover)\n", n_hub_est,
            n_hub_est / n_events))
if (n_hub_est > 0) {
  cat(sprintf("  → Implies: 1 handoff per %.0f spillovers\n", n_events / n_hub_est))
}
cat(sprintf("\n  Files saved:\n"))
cat(sprintf("    %s/results_table.rds\n", out_dir))
cat(sprintf("    %s/results_table.csv\n", out_dir))
if (length(hub_established_logs) > 0) {
  cat(sprintf("    %s/hub_established_details.rds\n", out_dir))
}
cat(sprintf("    %s/step5_params.rds\n", out_dir))
cat(sprintf("    + 8 diagnostic plots\n"))
cat("================================================================\n\n")