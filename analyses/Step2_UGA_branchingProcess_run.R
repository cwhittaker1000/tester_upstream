# =============================================================================
# branching_process.R
#
# Core simulation engine for Stages 3a, 3b, 3c of the spillover-pandemic
# pipeline. One generalised branching process handles all three modes:
#
#   Stage 3a (nonspatial_sub):  1-cell swath, no hub, R₀ < 1 → always dies out
#   Stage 3b (spatial):         Multi-cell swath with gravity kernel
#   Stage 3c (nonspatial_hub):  1-cell swath flagged as hub, R₀ ≥ 1
#
# The route classification from Stage 2 determines what inputs we construct,
# not a different algorithm. The same generation loop runs in all cases.
#
# Functions:
#   run_branching_process()  — the core generation-based SIR branching engine
#   make_nonspatial_swath()  — constructs trivial 1-cell input for 3a / 3c
#   simulate_event()         — dispatcher: takes event row + swath, runs the
#                              appropriate branching process variant
#
# Dependencies: none (base R only)
# =============================================================================
source("functions/branching_process.R")

# =============================================================================
# TESTING / DEMO
# =============================================================================
set.seed(42)

# --- Test 1: Single cell, subcritical (Stage 3a) ---------------------------
cat("--- Test 1: Subcritical (R₀ = 0.5, k = 1) ---\n")
ns <- make_nonspatial_swath(pop = 10000, is_hub = FALSE)
res <- run_branching_process(
  ns$N_pop, ns$hub_mask, ns$kernel, ns$spill_idx,
  R0 = 0.5, k = 1.0, log_gens = TRUE
)
cat(sprintf("  Outcome: %s | Chain: %d | Generations: %d\n", res$outcome, res$chain_size, res$generations))
stopifnot(res$outcome == "extinct")

# Run 1000 subcritical chains to check mean chain size ≈ 1/(1-R0) = 2
cat("  Running 1000 subcritical chains (R₀=0.5)...\n")
chains <- replicate(1000, {
  r <- run_branching_process(ns$N_pop, ns$hub_mask, ns$kernel, ns$spill_idx,
                             R0 = 0.5, k = 1.0)
  r$chain_size
})
cat(sprintf("  Mean chain size: %.2f (theory ≈ %.2f)\n", mean(chains), 1 / (1 - 0.5)))
cat(sprintf("  Max chain size:  %d\n\n", max(chains)))

# --- Test 2: Single cell, supercritical at hub (Stage 3c) ------------------
cat("--- Test 2: Supercritical at hub (R₀ = 1.5, k = 1) ---\n")
ns_hub <- make_nonspatial_swath(pop = 500000, is_hub = TRUE)

outcomes <- replicate(500, {
  r <- run_branching_process(
    ns_hub$N_pop, ns_hub$hub_mask, ns_hub$kernel, ns_hub$spill_idx,
    R0 = 1.5, k = 1.0, N_crit = 500
  )
  r$outcome
})

p_extinct <- mean(outcomes == "extinct")
p_hub     <- mean(outcomes == "hub_established")
cat(sprintf("  P(extinct) = %.3f (theory ≈ 0.667)\n", p_extinct))
cat(sprintf("  P(hub_est) = %.3f (theory ≈ 0.333)\n", p_hub))

# --- Test 3: Multi-cell spatial swath (Stage 3b) ---------------------------
cat("--- Test 3: Spatial branching (5-cell corridor, R₀ = 1.5) ---\n")

# Build a simple 5-cell linear corridor: [spill] -- [rural] -- [rural] -- [rural] -- [hub]
# Cells are 10km apart
n_test <- 5
pops_test <- c(500, 200, 300, 200, 50000)   # hub is cell 5
hub_test  <- c(FALSE, FALSE, FALSE, FALSE, TRUE)

# Distance matrix (linear, 10km spacing)
coords_test <- cbind(seq(0, 40, by = 10), rep(0, 5))
D_test <- as.matrix(dist(coords_test))

# Gravity kernel
alpha <- 1.0; beta <- 0.1; w_self <- 100
pop_attract_test <- pops_test^alpha
W_test <- sweep(exp(-beta * D_test), 2, pop_attract_test, FUN = "*")
diag(W_test) <- diag(W_test) + w_self
K_test <- W_test / rowSums(W_test)

# Capture fractions: assume full capture for this simple test
cf_test <- rep(1.0, n_test)

# Run one event with BOTH gen_log and cell_log
res_spatial <- run_branching_process(
  N_pop = pops_test, hub_mask = hub_test, kernel = K_test, spill_idx = 1L,
  R0 = 1.5, k = 1.0, N_crit = 100, capture_frac = cf_test, 
  log_gens = TRUE, log_cells = TRUE
)

cat(sprintf("  Outcome: %s | Chain: %d | Hub inf: %d | Gens: %d\n",
            res_spatial$outcome, res_spatial$chain_size,
            res_spatial$hub_infections, res_spatial$generations))
cat("  Generation log:\n")
print(res_spatial$gen_log)
cat("\n  Per-cell infection log (new infections per cell per generation):\n")
cat("  Columns: [spill] [rural] [rural] [rural] [HUB]\n")
print(res_spatial$cell_log)
cat("\n  Note: infections can flow backwards (cell 2 → cell 1, etc.)\n")
cat("  This is correct: the gravity kernel is isotropic (distance-based).\n\n")

# Run 500 spatial events to get outcome distribution
cat("  Running 500 spatial events...\n")
spatial_outcomes <- replicate(500, {
  r <- run_branching_process(
    N_pop = pops_test, hub_mask = hub_test, kernel = K_test, spill_idx = 1L,
    R0 = 1.5, k = 1.0, N_crit = 100, capture_frac = cf_test
  )
  c(outcome = r$outcome, chain = r$chain_size, hub = r$hub_infections)
})

outcome_vec <- spatial_outcomes["outcome", ]
chain_vec   <- as.integer(spatial_outcomes["chain", ])

cat(sprintf("  P(extinct):         %.3f\n", mean(outcome_vec == "extinct")))
cat(sprintf("  P(hub_established): %.3f\n", mean(outcome_vec == "hub_established")))
cat(sprintf("  P(generation_cap):  %.3f\n", mean(outcome_vec == "generation_cap")))
cat(sprintf("  Mean chain (extinct only): %.1f\n", mean(chain_vec[outcome_vec == "extinct"])))

# --- Test 4: simulate_event() dispatcher -----------------------------------
cat("--- Test 4: simulate_event() dispatcher ---\n")

# Mock event table rows
event_sc <- list(event_id = 1, R0 = 0.01, k = 1.0, route = "single_case")
event_ns <- list(event_id = 2, R0 = 0.5,  k = 1.0, route = "nonspatial_sub")
event_nh <- list(event_id = 3, R0 = 1.5,  k = 1.0, route = "nonspatial_hub")
event_sp <- list(event_id = 4, R0 = 1.5,  k = 1.0, route = "spatial")

# Mock swath for spatial event
mock_swath_spatial <- list(
  at_hub           = FALSE,
  cell_pop_total   = pops_test,
  hub_cell_mask    = hub_test,
  kernel_matrix    = K_test,
  spill_cell_idx   = 1L,
  capture_fraction = cf_test,
  hub_id           = 1
)

# Mock swath for at-hub event
mock_swath_hub <- list(
  at_hub = TRUE,
  hub_id = 1
)

# Mock hub_info
mock_hub_info <- data.frame(hub_id = 1, pop_k = 500, n_cells = 10,
                            lat = 0.3, lon = 32.6)

r1 <- simulate_event(event_sc, mock_swath_spatial, mock_hub_info)
cat(sprintf("  Single case:    outcome=%s, chain=%d\n", r1$outcome, r1$chain_size))
stopifnot(r1$outcome == "single_case", r1$chain_size == 1)

r2 <- simulate_event(event_ns, mock_swath_spatial, mock_hub_info)
cat(sprintf("  Nonspatial sub: outcome=%s, chain=%d\n", r2$outcome, r2$chain_size))
stopifnot(r2$outcome == "extinct")

r3 <- simulate_event(event_nh, mock_swath_hub, mock_hub_info)
cat(sprintf("  Nonspatial hub: outcome=%s, chain=%d, hub_inf=%d\n",
            r3$outcome, r3$chain_size, r3$hub_infections))

r4 <- simulate_event(event_sp, mock_swath_spatial, mock_hub_info,
                     log_gens = TRUE, log_cells = TRUE)
cat(sprintf("  Spatial:        outcome=%s, chain=%d, hub_inf=%d\n",
            r4$outcome, r4$chain_size, r4$hub_infections))

cat("  PASSED\n\n")
