# =============================================================================
# Step1_run_stages1_2.R
#
# Runs Stages 1 (viral property draws) and 2 (route determination) for all
# 500 spillover events in the Uganda swath library.
#
# Inputs:
#   - outputs/step0_UGA_preprocessing/swath_library.rds  (from Stage 0)
#   - analyses/draw_R0.R                                  (function definitions)
#
# Outputs (saved to outputs/step1_viral_draws/):
#   - event_table.rds       : data frame with R0, k, route for each event
#   - event_table.csv       : same, as CSV for inspection
#   - stage1_params.rds     : parameter values used (for reproducibility)
#   - plot_R0_distribution.png
#   - plot_k_distribution.png
#   - plot_route_summary.png
#   - plot_R0_by_route.png
#
# Usage:
#   source("analyses/Step1_run_stages1_2.R")
#   # or: Rscript analyses/Step1_run_stages1_2.R
# =============================================================================

library(ggplot2)

# --- Source function definitions ---------------------------------------------
source("functions/draw_R0.R")

# --- Configuration -----------------------------------------------------------

# Output directory
out_dir <- "outputs/step1_viral_draws"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Swath library from Stage 0
swath_file <- "outputs/step0_UGA_preprocessing/swath_library.rds"

# R₀ distribution: Zero-Inflated LogNormal
R0_method <- "ziln"
R0_mu     <- log(0.30)   # median of active component = 0.30
R0_sigma  <- 1.2          # log-scale spread
R0_p0     <- 0.70         # 70% dead-end spillovers (R₀ = 0)

# k distribution: Truncated Gamma
k_method  <- "truncated_gamma"
k_shape   <- 3.0          # Gamma shape (mode = (shape-1)/rate = 0.67)
k_rate    <- 3.0          # Gamma rate  (mean = shape/rate = 1.0)
k_max     <- 3.0          # hard ceiling — enforces overdispersion

# Route determination
R_min     <- 0.05         # below this → single case, skip branching

# Reproducibility
random_seed <- 42

# =============================================================================
# LOAD SWATH LIBRARY
# =============================================================================
if (!file.exists(swath_file)) {
  stop("Swath library not found at: ", swath_file,
       "\n  Run Stage 0 first (Step0_UGA_preprocessing.R)")
}

cat("Loading swath library:", swath_file, "\n")
swath_library <- readRDS(swath_file)
n_events <- length(swath_library$swaths)
n_at_hub <- sum(sapply(swath_library$swaths, function(s) s$at_hub))
cat(sprintf("  %d spillover events (%d at-hub, %d spatial)\n\n",
            n_events, n_at_hub, n_events - n_at_hub))

# =============================================================================
# RUN STAGES 1a, 1b, 2
# =============================================================================
set.seed(random_seed)
cat("Drawing R₀ from ZILN(p0=%.2f, μ=%.2f, σ=%.2f)...\n", R0_p0, R0_mu, R0_sigma)
cat("Drawing k from Truncated Gamma(shape=%.1f, rate=%.1f, max=%.1f)...\n\n",
    k_shape, k_rate, k_max)

event_table <- process_viral_draws_and_routing(
  swath_library,
  R0_method = R0_method,
  R0_mu     = R0_mu,
  R0_sigma  = R0_sigma,
  R0_p0     = R0_p0,
  k_method  = k_method,
  k_shape   = k_shape,
  k_rate    = k_rate,
  k_max     = k_max,
  R_min     = R_min,
  verbose   = TRUE
)

# =============================================================================
# ENRICH WITH SPATIAL METADATA FROM SWATHS
# =============================================================================
# Attach key swath info so downstream stages have everything in one table.

cat("Enriching event table with swath metadata...\n")
event_table$spill_lat    <- sapply(swath_library$swaths, function(s) s$spill_lat)
event_table$spill_lon    <- sapply(swath_library$swaths, function(s) s$spill_lon)
event_table$hub_id       <- sapply(swath_library$swaths, function(s) s$hub_id)
event_table$hub_dist_km  <- sapply(swath_library$swaths, function(s) {
  if (!is.null(s$length_km)) s$length_km else s$hub_dist_km
})
event_table$travel_time  <- sapply(swath_library$swaths, function(s) s$travel_time)

# For spatial swaths, record swath size and hub cell count
event_table$n_swath_cells <- sapply(swath_library$swaths, function(s) {
  if (!s$at_hub && !is.null(s$n_cells)) s$n_cells else NA_integer_
})
event_table$n_hub_cells <- sapply(swath_library$swaths, function(s) {
  if (!s$at_hub && !is.null(s$hub_cell_mask)) sum(s$hub_cell_mask) else NA_integer_
})

cat("  Done.\n\n")

# =============================================================================
# SAVE OUTPUTS
# =============================================================================

cat("Saving outputs to:", out_dir, "\n")

# Event table as RDS (preserves types, fast to reload)
saveRDS(event_table, file.path(out_dir, "event_table.rds"))
cat("  event_table.rds\n")

# Event table as CSV (human-readable)
write.csv(event_table, file.path(out_dir, "event_table.csv"), row.names = FALSE)
cat("  event_table.csv\n")

# Parameter record
params <- list(
  R0_method   = R0_method,
  R0_mu       = R0_mu,
  R0_sigma    = R0_sigma,
  R0_p0       = R0_p0,
  k_method    = k_method,
  k_shape     = k_shape,
  k_rate      = k_rate,
  k_max       = k_max,
  R_min       = R_min,
  random_seed = random_seed,
  n_events    = n_events,
  timestamp   = Sys.time()
)
saveRDS(params, file.path(out_dir, "stage1_params.rds"))
cat("  stage1_params.rds\n\n")

# =============================================================================
# DIAGNOSTIC PLOTS
# =============================================================================

cat("Generating diagnostic plots...\n\n")

# --- Colour palette for routes -----------------------------------------------
route_colours <- c(
  "single_case"    = "#999999",
  "nonspatial_sub" = "#56B4E9",
  "spatial"        = "#E69F00",
  "nonspatial_hub" = "#D55E00"
)

route_labels <- c(
  "single_case"    = "Single case (R₀ < 0.05)",
  "nonspatial_sub" = "Subcritical chain (0.05 ≤ R₀ < 1)",
  "spatial"        = "Spatial branching (R₀ ≥ 1, rural)",
  "nonspatial_hub" = "Hub branching (R₀ ≥ 1, at hub)"
)

# --- Plot 1: R₀ distribution ------------------------------------------------

# Build data for the histogram + annotation of the zero mass
R0_nonzero <- event_table$R0[event_table$R0 > 0]
n_zero     <- sum(event_table$R0 == 0)
frac_zero  <- n_zero / nrow(event_table)
p_gt1      <- mean(event_table$R0 > 1)

p1 <- ggplot(data.frame(R0 = R0_nonzero), aes(x = R0)) +
  geom_histogram(aes(y = after_stat(count)),
                 bins = 60, fill = "steelblue", colour = "white", alpha = 0.85) +
  scale_x_log10(
    breaks = c(0.001, 0.01, 0.05, 0.1, 0.3, 1, 3, 10),
    labels = c("0.001", "0.01", "0.05", "0.1", "0.3", "1", "3", "10")
  ) +
  geom_vline(xintercept = R_min, colour = "orange", linewidth = 0.9, linetype = "dashed") +
  geom_vline(xintercept = 1.0,   colour = "red",    linewidth = 0.9, linetype = "dashed") +
  annotate("text", x = R_min * 0.6, y = Inf, label = sprintf("R_min = %.2f", R_min),
           vjust = 1.8, colour = "orange", size = 3.5, fontface = "bold") +
  annotate("text", x = 1.5, y = Inf, label = "R₀ = 1",
           vjust = 1.8, colour = "red", size = 3.5, fontface = "bold") +
  annotate("label", x = max(R0_nonzero) * 0.5, y = Inf,
           label = sprintf("Zero-inflation: %d / %d events (%.0f%%)\nP(R₀ > 1) = %.2f%%",
                           n_zero, nrow(event_table), frac_zero * 100, p_gt1 * 100),
           vjust = 1.5, size = 3.2, fill = "white", alpha = 0.8) +
  labs(
    title    = "Stage 1a: R₀ Distribution (non-zero draws)",
    subtitle = sprintf("ZILN(p0=%.2f, μ=log(%.2f), σ=%.2f)  |  n = %d events",
                       R0_p0, exp(R0_mu), R0_sigma, n_events),
    x = expression(R[0] ~ "(log scale)"),
    y = "Count"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

p1

# --- Plot 2: k distribution -------------------------------------------------

p2 <- ggplot(event_table, aes(x = k)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 45, fill = "darkorange", colour = "white", alpha = 0.85) +
  geom_vline(xintercept = k_max, colour = "red", linewidth = 0.9, linetype = "dashed") +
  annotate("text", x = k_max - 0.1, y = Inf,
           label = sprintf("k_max = %.0f\n(overdispersion\nguaranteed)", k_max),
           vjust = 1.5, hjust = 1, colour = "red", size = 3.2, fontface = "bold") +
  # Overlay the theoretical truncated gamma density for reference
  stat_function(
    fun = function(x) {
      raw <- dgamma(x, shape = k_shape, rate = k_rate)
      # normalise to truncated support
      norm_const <- pgamma(k_max, shape = k_shape, rate = k_rate) -
        pgamma(0.01, shape = k_shape, rate = k_rate)
      ifelse(x >= 0.01 & x <= k_max, raw / norm_const, 0)
    },
    colour = "black", linewidth = 0.8, linetype = "solid", n = 300
  ) +
  labs(
    title    = "Stage 1b: Dispersion Parameter k",
    subtitle = sprintf("Truncated Gamma(shape=%.0f, rate=%.0f) on [0.01, %.0f]  |  median = %.2f",
                       k_shape, k_rate, k_max, median(event_table$k)),
    x = "k (NegBin dispersion — smaller = more superspreading)",
    y = "Density"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

p2

# --- Plot 3: Route summary (bar chart + pie-style breakdown) -----------------

route_counts <- as.data.frame(table(event_table$route))
names(route_counts) <- c("route", "count")
route_counts$pct   <- 100 * route_counts$count / sum(route_counts$count)
route_counts$label <- route_labels[as.character(route_counts$route)]
# Order by pipeline stage
route_order <- c("single_case", "nonspatial_sub", "spatial", "nonspatial_hub")
route_counts$route <- factor(route_counts$route, levels = route_order)
route_counts <- route_counts[order(route_counts$route), ]

p3 <- ggplot(route_counts, aes(x = route, y = count, fill = route)) +
  geom_col(width = 0.7, alpha = 0.9) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", count, pct)),
            vjust = -0.4, size = 4, fontface = "bold") +
  scale_fill_manual(values = route_colours, labels = route_labels, name = "Route") +
  scale_x_discrete(labels = route_labels) +
  labs(
    title    = "Stage 2: Route Determination",
    subtitle = sprintf("%d total events  |  R_min = %.2f", n_events, R_min),
    x = NULL,
    y = "Number of events"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    axis.text.x   = element_text(angle = 15, hjust = 1, size = 10),
    legend.position = "none"
  ) +
  coord_cartesian(ylim = c(0, max(route_counts$count) * 1.15))

p3

# --- Plot 4: R₀ by route (faceted densities) --------------------------------

# Exclude zero-R₀ events for the density plot (they would need a spike)
df_plot <- event_table[event_table$R0 > 0, ]
df_plot$route <- factor(df_plot$route, levels = route_order)
df_plot$route_label <- route_labels[as.character(df_plot$route)]

p4 <- ggplot(df_plot, aes(x = R0, fill = route)) +
  geom_histogram(aes(y = after_stat(count)),
                 bins = 40, colour = "white", alpha = 0.85) +
  scale_x_log10(
    breaks = c(0.01, 0.05, 0.1, 0.3, 1, 3, 10),
    labels = c("0.01", "0.05", "0.1", "0.3", "1", "3", "10")
  ) +
  scale_fill_manual(values = route_colours) +
  geom_vline(xintercept = 1.0, colour = "red", linewidth = 0.7, linetype = "dashed") +
  facet_wrap(~ route_label, scales = "free_y", ncol = 2) +
  labs(
    title    = "R₀ Distribution by Route (non-zero draws only)",
    subtitle = sprintf("%d zero-R₀ dead-ends not shown (route: single_case)", n_zero),
    x = expression(R[0] ~ "(log scale)"),
    y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "none",
    strip.text      = element_text(face = "bold", size = 10)
  )

p4

# --- Plot 5: k vs R₀ scatter (coloured by route) ----------------------------

df_scatter <- event_table[event_table$R0 > 0, ]
df_scatter$route <- factor(df_scatter$route, levels = route_order)

p5 <- ggplot(df_scatter, aes(x = R0, y = k, colour = route)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_x_log10(
    breaks = c(0.01, 0.05, 0.1, 0.3, 1, 3, 10),
    labels = c("0.01", "0.05", "0.1", "0.3", "1", "3", "10")
  ) +
  scale_colour_manual(values = route_colours, labels = route_labels, name = "Route") +
  geom_vline(xintercept = R_min, colour = "orange", linewidth = 0.6, linetype = "dashed") +
  geom_vline(xintercept = 1.0,   colour = "red",    linewidth = 0.6, linetype = "dashed") +
  geom_hline(yintercept = k_max, colour = "red",    linewidth = 0.6, linetype = "dotted") +
  labs(
    title    = "Joint R₀ × k Space (non-zero R₀ events)",
    subtitle = "Each point is one spillover event; k and R₀ are drawn independently",
    x = expression(R[0] ~ "(log scale)"),
    y = "k (dispersion)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title    = element_text(face = "bold")
  )

p5

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  Stages 1–2 complete\n")
cat("================================================================\n")
cat(sprintf("  Events processed:  %d\n", nrow(event_table)))
cat(sprintf("  Output directory:  %s\n", out_dir))
cat(sprintf("\n  Route breakdown:\n"))
for (r in route_order) {
  n_r <- sum(event_table$route == r)
  cat(sprintf("    %-18s  %4d  (%5.1f%%)\n", r, n_r, 100 * n_r / nrow(event_table)))
}
cat(sprintf("\n  Events entering Stage 3b (spatial):      %d\n",
            sum(event_table$route == "spatial")))
cat(sprintf("  Events entering Stage 3c (hub branching): %d\n",
            sum(event_table$route == "nonspatial_hub")))
cat(sprintf("\n  Files saved:\n"))
cat(sprintf("    %s/event_table.rds   (reload with readRDS())\n", out_dir))
cat(sprintf("    %s/event_table.csv   (human-readable)\n", out_dir))
cat(sprintf("    %s/stage1_params.rds (parameter record)\n", out_dir))
cat(sprintf("    + 5 diagnostic plots\n"))
cat("\n  Next: source('analyses/Step2_branching_process.R')\n")
cat("================================================================\n\n")