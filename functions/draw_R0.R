# =============================================================================
# STAGE 1a: Draw R₀
# =============================================================================
#
# Two distribution options:
#
#   (A) LogNormal(μ, σ)
#       - log(R₀) ~ Normal(μ, σ)
#       - Median R₀ = exp(μ), so μ = log(median)
#       - σ controls tail weight: larger σ → more P(R₀ > 1)
#       - Calibration target: P(R₀ > 1) should yield ~4 pandemics per 125 yr
#         given global spillover volume
#
#   (B) Zero-Inflated LogNormal (ZILN)
#       - With probability p0, R₀ = 0  (dead-end spillover, no human transmission)
#       - With probability (1 - p0), R₀ ~ LogNormal(μ, σ)
#       - Rationale: most zoonotic viruses have no meaningful human-to-human
#         transmission at all; the lognormal component represents the subset
#         that achieve any transmissibility.
#       - p0 effectively compresses the left tail into a point mass, letting
#         the lognormal component have a higher median while preserving the
#         overall rarity of R₀ > 1 events.
#
# Computational shortcut: R₀ < R_min → record as single case, skip branching.
# =============================================================================


#' Draw R₀ values for spillover events
#'
#' @param n          Number of R₀ values to draw
#' @param method     "lognormal" or "ziln" (zero-inflated lognormal)
#' @param mu         Log-scale location parameter. Median of lognormal = exp(mu).
#' @param sigma      Log-scale spread parameter. Controls right-tail weight.
#' @param p0         Zero-inflation probability (only used if method = "ziln").
#'                   Fraction of spillovers with R₀ = 0 (no human transmission).
#' @return Numeric vector of length n with R₀ values (≥ 0)
#'
#' @details
#' Default mu = log(0.15) ≈ -1.90, giving median R₀ ≈ 0.15 (most spillovers
#' transmit poorly). Default sigma = 1.5, giving P(R₀ > 1) ≈ 0.01–0.03
#' depending on method. These are starting points for calibration.
#'
#' For the ZILN method, p0 = 0.7 means 70% of spillovers are complete dead
#' ends (R₀ = 0), and the remaining 30% draw from LogNormal(mu, sigma).
#' The lognormal component can then have a somewhat higher median (e.g., 0.3)
#' while the overall P(R₀ > 1) remains small.
draw_R0 <- function(n,
                    method = c("lognormal", "ziln"),
                    mu     = log(0.15),
                    sigma  = 1.5,
                    p0     = 0.7) {
  
  method <- match.arg(method)
  
  if (method == "lognormal") {
    # --- Standard LogNormal ---
    R0 <- rlnorm(n, meanlog = mu, sdlog = sigma)
    
  } else if (method == "ziln") {
    # --- Zero-Inflated LogNormal ---
    # Step 1: decide which events are dead-ends (R₀ = 0)
    is_dead_end <- rbinom(n, size = 1, prob = p0) == 1
    
    # Step 2: draw lognormal R₀ for the non-dead-end events
    R0 <- numeric(n)
    n_active <- sum(!is_dead_end)
    if (n_active > 0) {
      R0[!is_dead_end] <- rlnorm(n_active, meanlog = mu, sdlog = sigma)
    }
    # Dead-end events stay at R₀ = 0
  }
  
  return(R0)
}

#' Report summary statistics of an R₀ distribution
#'
#' Useful for calibration: prints median, key quantiles, P(R₀ > 1),
#' and implied pandemic frequency given a spillover rate.
#'
#' @param R0_vec                Numeric vector of R₀ draws
#' @param annual_spillovers     Expected number of spillovers per year
#'                              (for pandemic frequency calculation)
#' @param p_survive_bottleneck  Approximate probability that an R₀ > 1 event
#'                              survives the spatial bottleneck to become a
#'                              pandemic. Default 0.1 (rough placeholder).
summarise_R0_distribution <- function(R0_vec,
                                      annual_spillovers     = 50,
                                      p_survive_bottleneck  = 0.1) {
  
  cat("\n--- R₀ Distribution Summary ---\n")
  cat(sprintf("  n draws:       %d\n", length(R0_vec)))
  cat(sprintf("  Median R₀:     %.4f\n", median(R0_vec)))
  cat(sprintf("  Mean R₀:       %.4f\n", mean(R0_vec)))
  cat(sprintf("  Max R₀:        %.4f\n", max(R0_vec)))
  
  cat("  Quantiles:\n")
  probs <- c(0.50, 0.90, 0.95, 0.99, 0.999, 1.0)
  qs <- quantile(R0_vec, probs = probs)
  for (i in seq_along(qs)) {
    cat(sprintf("    %6.2f%%:  %.4f\n", probs[i] * 100, qs[i]))
  }
  
  p_gt1 <- mean(R0_vec > 1)
  p_gt2 <- mean(R0_vec > 2)
  cat(sprintf("\n  P(R₀ > 1):     %.6f  (1 in %.0f)\n", p_gt1,
              ifelse(p_gt1 > 0, 1 / p_gt1, Inf)))
  cat(sprintf("  P(R₀ > 2):     %.6f  (1 in %.0f)\n", p_gt2,
              ifelse(p_gt2 > 0, 1 / p_gt2, Inf)))
  
}


# =============================================================================
# STAGE 1b: Draw Dispersion Parameter k
# =============================================================================
#
# The offspring distribution is NegBin(mean = R_eff, size = k).
# Smaller k → more overdispersion (superspreading).
# k → ∞ recovers Poisson (no overdispersion).
#
# For influenza, literature estimates: k ≈ 0.5–2.
# We want to ENFORCE overdispersion by bounding k from above.
#
# Options:
#   (A) Fixed k = k0 (simplest; appropriate if k is well-known)
#   (B) Draw from Truncated Gamma on (0, k_max]
#       - Gamma(shape, rate) truncated at k_max
#       - shape and rate chosen so the bulk of the distribution is in [0.3, 2]
#       - k_max = 3 (hard ceiling) ensures we never enter near-Poisson territory
#
# We use rejection sampling for the truncated gamma: draw from Gamma,
# reject if > k_max. Since k_max is well into the tail, rejection rate is low.
# =============================================================================

#' Draw dispersion parameter k for NegBin offspring distribution
#'
#' @param n       Number of values to draw
#' @param method  "fixed" or "truncated_gamma"
#' @param k0      Fixed value (used if method = "fixed"). Default 1.0.
#' @param shape   Shape parameter of the Gamma distribution. Default 3.0.
#' @param rate    Rate parameter of the Gamma distribution. Default 3.0.
#'                With shape=3, rate=3: mean=1.0, mode=0.67, sd=0.58.
#' @param k_max   Upper bound (hard ceiling). Default 3.0.
#'                NegBin with k=3 still shows meaningful overdispersion:
#'                variance = mean + mean²/3 ≈ 1.33× Poisson variance at mean=1.
#' @param k_min   Lower bound (floor). Default 0.01.
#'                Prevents numerical issues with extremely small k.
#' @return Numeric vector of length n with k values in [k_min, k_max]
draw_k <- function(n,
                   method = c("fixed", "truncated_gamma"),
                   k0     = 1.0,
                   shape  = 3.0,
                   rate   = 3.0,
                   k_max  = 3.0,
                   k_min  = 0.01) {
  
  method <- match.arg(method)
  
  if (method == "fixed") {
    return(rep(k0, n))
  }
  
  # --- Truncated Gamma via rejection sampling ---
  # Draw from Gamma(shape, rate), reject values outside [k_min, k_max].
  # Pre-compute the acceptance probability to warn if it's very low.
  p_accept <- pgamma(k_max, shape = shape, rate = rate) -
    pgamma(k_min, shape = shape, rate = rate)
  
  if (p_accept < 0.5) {
    warning(sprintf(
      "Low acceptance rate (%.1f%%) for truncated gamma with shape=%.1f, rate=%.1f, k_max=%.1f. Consider adjusting parameters.",
      p_accept * 100, shape, rate, k_max
    ))
  }
  
  k_vals <- numeric(n)
  n_filled <- 0
  max_attempts <- n * 20  # safety valve
  attempts <- 0
  
  while (n_filled < n && attempts < max_attempts) {
    # Draw a batch
    batch_size <- min((n - n_filled) * 3, 10000)  # overshoot to reduce loops
    candidates <- rgamma(batch_size, shape = shape, rate = rate)
    
    # Accept those within bounds
    valid <- candidates[candidates >= k_min & candidates <= k_max]
    
    if (length(valid) > 0) {
      n_to_take <- min(length(valid), n - n_filled)
      k_vals[(n_filled + 1):(n_filled + n_to_take)] <- valid[1:n_to_take]
      n_filled <- n_filled + n_to_take
    }
    
    attempts <- attempts + batch_size
  }
  
  if (n_filled < n) {
    warning(sprintf("Only drew %d / %d k values after %d attempts.", n_filled, n, attempts))
    k_vals[(n_filled + 1):n] <- k0  # fill remainder with default
  }
  
  return(k_vals)
}


# =============================================================================
# STAGE 2: Route Determination
# =============================================================================
#
# Given R₀ and swath properties, classify each spillover event into one of
# four routes:
#
#   Route              Condition                        Next stage
#   ─────────────────  ───────────────────────────────  ──────────
#   "single_case"      R₀ < R_min (~0.05)               Stage 5 (1 death × IFR)
#   "nonspatial_sub"   R_min ≤ R₀ < 1                   Stage 3a (chain dies out)
#   "spatial"          R₀ ≥ 1 AND not at-hub             Stage 3b (spatial branching)
#   "nonspatial_hub"   R₀ ≥ 1 AND at-hub                Stage 3c (non-spatial at hub)
#
# The route fully determines which branching process variant to run.
# =============================================================================

#' Determine route for a single spillover event
#'
#' @param R0     Reproduction number drawn in Stage 1a
#' @param at_hub Logical: is the spillover location within a hub patch?
#' @param R_min  Minimum R₀ threshold below which transmission is negligible.
#'               Default 0.05. Events below this are recorded as single cases
#'               with no secondary transmission.
#' @return Character string: one of "single_case", "nonspatial_sub",
#'         "spatial", "nonspatial_hub"
determine_route <- function(R0, at_hub, R_min = 0.05) {
  
  if (R0 < R_min) {
    return("single_case")
  } else if (R0 < 1.0) {
    return("nonspatial_sub")
  } else {
    # R₀ ≥ 1: supercritical. Route depends on spatial location.
    if (at_hub) {
      return("nonspatial_hub")
    } else {
      return("spatial")
    }
  }
}

# =============================================================================
# COMBINED: Process a full batch of spillover events (Stages 1–2)
# =============================================================================

#' Run Stages 1a, 1b, and 2 for all spillover events in a swath library
#'
#' @param swath_library  List loaded from swath_library.rds (contains $swaths)
#' @param R0_method      "lognormal" or "ziln"
#' @param R0_mu          Log-scale location for R₀ distribution
#' @param R0_sigma       Log-scale spread for R₀ distribution
#' @param R0_p0          Zero-inflation probability (only for "ziln")
#' @param k_method       "fixed" or "truncated_gamma"
#' @param k_fixed        Fixed k value (if k_method = "fixed")
#' @param k_shape        Gamma shape (if k_method = "truncated_gamma")
#' @param k_rate         Gamma rate (if k_method = "truncated_gamma")
#' @param k_max          Upper bound on k
#' @param R_min          Threshold below which R₀ is treated as single case
#' @param verbose        Print summary statistics?
#'
#' @return Data frame with one row per spillover event:
#'   - event_id:  spillover index (1..n)
#'   - R0:        drawn reproduction number
#'   - k:         drawn dispersion parameter
#'   - at_hub:    logical, from swath library
#'   - route:     one of "single_case", "nonspatial_sub", "spatial", "nonspatial_hub"
process_viral_draws_and_routing <- function(
    swath_library,
    R0_method  = "lognormal",
    R0_mu      = log(0.15),
    R0_sigma   = 1.5,
    R0_p0      = 0.7,
    k_method   = "fixed",
    k_fixed    = 1.0,
    k_shape    = 3.0,
    k_rate     = 3.0,
    k_max      = 3.0,
    R_min      = 0.05,
    verbose    = TRUE
) {
  
  swaths <- swath_library$swaths
  n <- length(swaths)
  
  # --- Stage 1a: Draw R₀ for each event ---
  R0_vec <- draw_R0(n,
                    method = R0_method,
                    mu     = R0_mu,
                    sigma  = R0_sigma,
                    p0     = R0_p0)
  
  # --- Stage 1b: Draw k for each event ---
  k_vec <- draw_k(n,
                  method = k_method,
                  k0     = k_fixed,
                  shape  = k_shape,
                  rate   = k_rate,
                  k_max  = k_max)
  
  # --- Stage 2: Route determination ---
  at_hub_vec <- sapply(swaths, function(s) s$at_hub)
  route_vec  <- character(n)
  
  for (i in 1:n) {
    route_vec[i] <- determine_route(R0_vec[i], at_hub_vec[i], R_min)
  }
  
  # --- Build result table ---
  results <- data.frame(
    event_id = 1:n,
    R0       = R0_vec,
    k        = k_vec,
    at_hub   = at_hub_vec,
    route    = route_vec,
    stringsAsFactors = FALSE
  )
  
  # --- Summary ---
  if (verbose) {
    summarise_R0_distribution(R0_vec = R0_vec)
    
    cat("--- Route Determination Summary ---\n")
    route_table <- table(results$route)
    for (r in names(route_table)) {
      cat(sprintf("  %-18s  %4d  (%5.1f%%)\n",
                  r, route_table[r], 100 * route_table[r] / n))
    }
    
    cat(sprintf("\n  k summary: median = %.2f, range = [%.3f, %.3f]\n",
                median(k_vec), min(k_vec), max(k_vec)))
    cat("---\n\n")
  }
  
  return(results)
}




# =============================================================================
# DIAGNOSTIC: Visualise the R₀ and k distributions
# =============================================================================

#' Plot R₀ distribution with key thresholds marked
#'
#' @param R0_vec   Numeric vector of R₀ draws
#' @param R_min    Single-case threshold
#' @param title    Plot title
plot_R0_distribution <- function(R0_vec, R_min = 0.05, title = NULL) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    # Fallback to base R
    hist(log10(R0_vec[R0_vec > 0] + 1e-6), breaks = 80,
         main = ifelse(is.null(title), "R₀ Distribution (log10 scale)", title),
         xlab = "log10(R₀)", col = "steelblue", border = "white")
    abline(v = log10(R_min), col = "orange", lwd = 2, lty = 2)
    abline(v = 0, col = "red", lwd = 2, lty = 2)  # R₀ = 1
    legend("topright",
           legend = c(sprintf("R_min = %.2f", R_min), "R₀ = 1 (pandemic threshold)"),
           col = c("orange", "red"), lty = 2, lwd = 2)
    return(invisible(NULL))
  }
  
  library(ggplot2)
  
  df <- data.frame(R0 = R0_vec)
  
  # Separate zero and non-zero for display
  n_zero <- sum(R0_vec == 0)
  frac_zero <- n_zero / length(R0_vec)
  
  p <- ggplot(df[df$R0 > 0, , drop = FALSE], aes(x = R0)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 80, fill = "steelblue", colour = "white", alpha = 0.8) +
    scale_x_log10(
      limits = c(1e-4, max(10, max(R0_vec) * 1.2)),
      breaks = c(0.001, 0.01, 0.1, 1, 10),
      labels = c("0.001", "0.01", "0.1", "1", "10")
    ) +
    geom_vline(xintercept = R_min, colour = "orange", linewidth = 1, linetype = "dashed") +
    geom_vline(xintercept = 1.0, colour = "red", linewidth = 1, linetype = "dashed") +
    annotate("text", x = R_min, y = Inf, label = paste0("R_min = ", R_min),
             vjust = 2, hjust = -0.1, colour = "orange", size = 3.5) +
    annotate("text", x = 1.0, y = Inf, label = "R₀ = 1",
             vjust = 2, hjust = -0.1, colour = "red", size = 3.5) +
    labs(
      title = ifelse(is.null(title), "R₀ Distribution (log scale)", title),
      subtitle = sprintf(
        "P(R₀ = 0) = %.1f%% | P(R₀ > 1) = %.4f%% | median(R₀ > 0) = %.3f",
        frac_zero * 100,
        mean(R0_vec > 1) * 100,
        median(R0_vec[R0_vec > 0])
      ),
      x = expression(R[0]),
      y = "Density (non-zero draws only)"
    ) +
    theme_minimal(base_size = 13)
  
  print(p)
  return(invisible(p))
}


#' Plot k distribution
#'
#' @param k_vec  Numeric vector of k draws
#' @param k_max  Upper bound used in sampling
plot_k_distribution <- function(k_vec, k_max = 3.0) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    hist(k_vec, breaks = 40, main = "Dispersion Parameter k",
         xlab = "k", col = "darkorange", border = "white")
    abline(v = k_max, col = "red", lwd = 2, lty = 2)
    return(invisible(NULL))
  }
  
  library(ggplot2)
  
  df <- data.frame(k = k_vec)
  
  p <- ggplot(df, aes(x = k)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 50, fill = "darkorange", colour = "white", alpha = 0.8) +
    geom_vline(xintercept = k_max, colour = "red", linewidth = 1, linetype = "dashed") +
    annotate("text", x = k_max, y = Inf,
             label = paste0("k_max = ", k_max, " (overdispersion floor)"),
             vjust = 2, hjust = 1.1, colour = "red", size = 3.5) +
    labs(
      title = "Dispersion Parameter k",
      subtitle = sprintf("median = %.2f | range = [%.3f, %.3f] | smaller k → more superspreading",
                         median(k_vec), min(k_vec), max(k_vec)),
      x = "k (NegBin dispersion)",
      y = "Density"
    ) +
    theme_minimal(base_size = 13)
  
  print(p)
  return(invisible(p))
}

