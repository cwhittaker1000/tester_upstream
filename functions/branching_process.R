# =============================================================================
# run_branching_process()
# =============================================================================
#
# Runs a discrete-generation branching process on a swath (1 to N cells).
#
# Each generation g:
#   1. For each cell j with I_j(g) > 0:
#      a) Compute R_eff = R0 * S_j(g) / N_j
#      b) Draw total offspring:  NegBin(I_j(g) * R_eff,  I_j(g) * k)
#         (vectorised sum of I_j(g) iid NegBin(R_eff, k) draws)
#      c) Thin by capture fraction:  Binomial(N_offspring, capture_frac[j])
#      d) Distribute captured offspring across cells: Multinomial(N_captured, P[j,])
#   2. For each destination cell ℓ with assigned offspring:
#      - Cap at S_ℓ (can't infect more than remaining susceptibles)
#      - Decrement S_ℓ, increment I_ℓ(g+1), accumulate C_hub if ℓ is hub
#   3. Zero out I_j(g) for all cells (recovered / removed)
#   4. Check termination conditions
#
# @param N_pop        Numeric vector length n_cells: total population per cell
# @param hub_mask     Logical vector length n_cells: TRUE for hub cells
# @param kernel       n_cells × n_cells row-stochastic matrix P[j,ℓ].
#                     For 1-cell swaths this is matrix(1, 1, 1).
# @param spill_idx    Integer: index of the spillover cell (1-based)
# @param R0           Numeric: basic reproduction number for this event
# @param k            Numeric: NegBin dispersion parameter
# @param N_crit       Integer: cumulative hub infections to trigger handoff.
#                     Default 500.
# @param g_max        Integer: maximum generations before forced termination.
#                     Default 100.
# @param capture_frac Numeric vector length n_cells: fraction of offspring 
#                     retained in swath (1.0 = no leakage). If NULL, no
#                     thinning is applied (appropriate for 1-cell swaths).
# @param log_gens     Logical: if TRUE, record per-generation diagnostics.
#                     Default FALSE for speed in production runs.
# @param log_cells    Logical: if TRUE, record infections per cell per 
#                     generation as a matrix. Useful for visualising wavefront
#                     propagation through the corridor. Default FALSE.
# @param offspring_cap Integer: maximum total offspring per cell per generation.
#                     Safety valve against extreme NegBin draws. Default 10000.
#
# @return List with:
#   $outcome         Character: "extinct", "hub_established", or "generation_cap"
#   $chain_size      Integer: total cumulative infections across all cells
#   $hub_infections  Integer: cumulative infections in hub cells only
#   $generations     Integer: number of generations elapsed
#   $peak_infected   Integer: maximum simultaneously infected across any generation
#   $gen_log         Data frame (only if log_gens = TRUE): per-generation stats
#   $cell_log        Matrix (only if log_cells = TRUE): rows = generations,
#                    cols = cells. Entry [g, j] = new infections in cell j at gen g.
# =============================================================================

run_branching_process <- function(N_pop,
                                  hub_mask,
                                  kernel,
                                  spill_idx,
                                  R0,
                                  k,
                                  N_crit       = 500,
                                  g_max        = 100,
                                  capture_frac = NULL,
                                  log_gens     = FALSE,
                                  log_cells    = FALSE,
                                  offspring_cap = 10000) {
  
  n_cells <- length(N_pop)
  
  # --- Validate inputs -------------------------------------------------------
  stopifnot(length(hub_mask) == n_cells)
  stopifnot(nrow(kernel) == n_cells && ncol(kernel) == n_cells)
  stopifnot(spill_idx >= 1 && spill_idx <= n_cells)
  stopifnot(R0 >= 0, k > 0)
  
  if (is.null(capture_frac)) {
    capture_frac <- rep(1.0, n_cells)  # no thinning
  }
  stopifnot(length(capture_frac) == n_cells)
  
  # --- Initialise state vectors ----------------------------------------------
  S <- N_pop                        # susceptibles per cell
  I <- integer(n_cells)             # actively infected this generation
  I[spill_idx] <- 1L                # seed: 1 infected in spillover cell
  S[spill_idx] <- S[spill_idx] - 1  # that individual came from the susceptibles
  
  C_total <- 1L                     # cumulative total infections
  C_hub   <- 0L                     # cumulative hub infections
  
  # Count the initial case toward hub infections if spillover is in a hub cell
  if (hub_mask[spill_idx]) {
    C_hub <- 1L
  }
  
  peak_infected <- 1L
  
  # --- Optional generation log -----------------------------------------------
  if (log_gens) {
    gen_log <- data.frame(
      generation     = 0L,
      total_infected = 1L,
      new_infections = 1L,
      hub_infections = C_hub,
      active_cells   = 1L,
      cum_total      = C_total,
      cum_hub        = C_hub
    )
  }
  
  # --- Optional per-cell log -------------------------------------------------
  # Matrix where row g = generation, col j = cell. Stores new infections.
  # Pre-allocate with g_max+1 rows; we trim at the end.
  if (log_cells) {
    cell_log <- matrix(0L, nrow = g_max + 1, ncol = n_cells)
    cell_log[1, spill_idx] <- 1L  # generation 0: initial case
  }
  
  # --- Generation loop -------------------------------------------------------
  for (g in 1:g_max) {
    
    # Which cells have active infections this generation?
    active <- which(I > 0)
    
    if (length(active) == 0) {
      # Extinction — no active infections remain
      outcome <- "extinct"
      break
    }
    
    # Next generation's infections (accumulator)
    I_next <- integer(n_cells)
    new_infections_this_gen <- 0L
    
    # --- Process each active cell --------------------------------------------
    for (j in active) {
      
      n_infected_j <- I[j]
      
      # Effective R in this cell: R0 * S_j / N_j
      # If N_pop is 0 (shouldn't happen) or S is 0, R_eff = 0
      if (N_pop[j] <= 0 || S[j] <= 0) next
      
      R_eff <- R0 * S[j] / N_pop[j]
      
      if (R_eff <= 0) next
      
      # --- Draw total offspring from this cell --------------------------------
      # Sum of n_infected_j iid NegBin(R_eff, k) = NegBin(n * R_eff, n * k)
      # R's rnbinom: size = dispersion, mu = mean
      total_offspring <- rnbinom(1, 
                                 size = n_infected_j * k, 
                                 mu   = n_infected_j * R_eff)
      
      if (total_offspring == 0) next
      
      # Safety cap
      total_offspring <- min(total_offspring, offspring_cap)
      
      # --- Thin by capture fraction ------------------------------------------
      # Accounts for offspring that would land outside the swath
      if (capture_frac[j] < 1.0) {
        total_offspring <- rbinom(1, size = total_offspring, prob = capture_frac[j])
        if (total_offspring == 0) next
      }
      
      # --- Distribute offspring across destination cells ----------------------
      # Multinomial draw using pre-computed kernel row for cell j
      if (n_cells == 1) {
        # Trivial case: all offspring stay in the single cell
        dest_counts <- total_offspring
      } else {
        dest_counts <- as.integer(rmultinom(1, size = total_offspring, 
                                            prob = kernel[j, ]))
      }
      
      # --- Apply to destination cells: check susceptibles --------------------
      for (ell in which(dest_counts > 0)) {
        
        assigned <- dest_counts[ell]
        
        # Can't infect more than remaining susceptibles
        actual <- min(assigned, S[ell])
        
        if (actual > 0) {
          I_next[ell]  <- I_next[ell] + actual
          S[ell]       <- S[ell] - actual
          new_infections_this_gen <- new_infections_this_gen + actual
          
          # Track hub infections
          if (hub_mask[ell]) {
            C_hub <- C_hub + actual
          }
        }
      }
    }
    
    # --- Update state for next generation ------------------------------------
    I <- I_next
    C_total <- C_total + new_infections_this_gen
    
    total_active <- sum(I)
    peak_infected <- max(peak_infected, total_active)
    
    # --- Log per-cell infections this generation -----------------------------
    if (log_cells) {
      cell_log[g + 1, ] <- I_next  # row g+1 because row 1 = generation 0
    }
    
    # --- Log this generation -------------------------------------------------
    if (log_gens) {
      gen_log <- rbind(gen_log, data.frame(
        generation     = g,
        total_infected = total_active,
        new_infections = new_infections_this_gen,
        hub_infections = sum(I[hub_mask]),
        active_cells   = sum(I > 0),
        cum_total      = C_total,
        cum_hub        = C_hub
      ))
    }
    
    # --- Check termination conditions ----------------------------------------
    
    # Extinction
    if (total_active == 0) {
      outcome <- "extinct"
      break
    }
    
    # Hub establishment
    if (C_hub >= N_crit) {
      outcome <- "hub_established"
      break
    }
    
    # Generation cap (only triggers on last iteration of the loop)
    if (g == g_max) {
      outcome <- "generation_cap"
    }
  }
  
  # Handle case where loop completes without breaking (g == g_max and still active)
  if (!exists("outcome")) {
    outcome <- "generation_cap"
  }
  
  # --- Build return value ----------------------------------------------------
  result <- list(
    outcome        = outcome,
    chain_size     = C_total,
    hub_infections = C_hub,
    generations    = g,
    peak_infected  = peak_infected
  )
  
  if (log_gens) {
    result$gen_log <- gen_log
  }
  
  if (log_cells) {
    # Trim to actual generations elapsed (g+1 rows: gen 0 through gen g)
    result$cell_log <- cell_log[1:(g + 1), , drop = FALSE]
    rownames(result$cell_log) <- paste0("g", 0:g)
  }
  
  return(result)
}


# =============================================================================
# make_nonspatial_swath()
# =============================================================================
#
# Constructs a trivial 1-cell "swath" for non-spatial branching processes:
#   - Stage 3a (subcritical, R₀ < 1): large population, no hub → always extinct
#   - Stage 3c (at-hub, R₀ ≥ 1): large hub population, hub flagged → can 
#     trigger establishment
#
# @param pop       Total population of the single cell. For 3a, a moderate
#                  value (10000) suffices since depletion is negligible. For
#                  3c, use the actual hub population from swath metadata.
# @param is_hub    Logical: should this cell be flagged as hub?
#                  FALSE for Stage 3a, TRUE for Stage 3c.
#
# @return List with the same structure expected by run_branching_process():
#   $N_pop, $hub_mask, $kernel, $spill_idx, $capture_frac
# =============================================================================

make_nonspatial_swath <- function(pop = 10000, is_hub = FALSE) {
  
  list(
    N_pop        = pop,
    hub_mask     = is_hub,
    kernel       = matrix(1.0, nrow = 1, ncol = 1),
    spill_idx    = 1L,
    capture_frac = 1.0  # no leakage in a 1-cell swath
  )
}

# =============================================================================
# simulate_event()
# =============================================================================
#
# Dispatcher: given one event row (from the Stage 1-2 event table) and the
# corresponding swath (from swath_library), run the appropriate simulation
# and return a result.
#
# Routes:
#   "single_case"    → no branching process. Return chain_size = 1.
#   "nonspatial_sub" → build 1-cell non-hub swath, run branching process.
#   "spatial"        → use the full multi-cell swath, run branching process.
#   "nonspatial_hub" → build 1-cell hub swath, run branching process.
#
# @param event       One-row data frame (or named list) with fields:
#                    R0, k, route, event_id
# @param swath       List element from swath_library$swaths[[i]]. Must contain
#                    cell_pop_total, hub_cell_mask, kernel_matrix, spill_cell_idx,
#                    capture_fraction (for spatial swaths).
# @param hub_info    Data frame of hub metadata from swath_library$hub_info,
#                    used to look up hub population for at-hub events.
# @param N_crit      Hub establishment threshold (cumulative infections).
# @param g_max       Maximum generations.
# @param log_gens    Whether to record per-generation diagnostics.
# @param nonspatial_pop  Default population for non-spatial subcritical swaths.
#
# @return List with:
#   $event_id, $route, $R0, $k,
#   $outcome, $chain_size, $hub_infections, $generations, $peak_infected,
#   $gen_log (if log_gens = TRUE)
# =============================================================================

simulate_event <- function(event,
                           swath,
                           hub_info        = NULL,
                           N_crit          = 500,
                           g_max           = 100,
                           log_gens        = FALSE,
                           log_cells       = FALSE,
                           nonspatial_pop  = 10000) {
  
  route <- event$route
  R0    <- event$R0
  k_val <- event$k
  
  # --- Single case: no branching needed --------------------------------------
  if (route == "single_case") {
    result <- list(
      event_id       = event$event_id,
      route          = route,
      R0             = R0,
      k              = k_val,
      outcome        = "single_case",
      chain_size     = 1L,
      hub_infections = 0L,
      generations    = 0L,
      peak_infected  = 1L
    )
    if (log_gens) {
      result$gen_log <- data.frame(
        generation = 0L, total_infected = 1L, new_infections = 1L,
        hub_infections = 0L, active_cells = 1L, cum_total = 1L, cum_hub = 0L
      )
    }
    if (log_cells) {
      result$cell_log <- matrix(1L, nrow = 1, ncol = 1,
                                dimnames = list("g0", "cell1"))
    }
    return(result)
  }
  
  # --- Non-spatial subcritical (Stage 3a) ------------------------------------
  if (route == "nonspatial_sub") {
    ns <- make_nonspatial_swath(pop = nonspatial_pop, is_hub = FALSE)
    
    bp_result <- run_branching_process(
      N_pop        = ns$N_pop,
      hub_mask     = ns$hub_mask,
      kernel       = ns$kernel,
      spill_idx    = ns$spill_idx,
      R0           = R0,
      k            = k_val,
      N_crit       = N_crit,
      g_max        = g_max,
      capture_frac = ns$capture_frac,
      log_gens     = log_gens,
      log_cells    = log_cells
    )
  }
  
  # --- Spatial branching process (Stage 3b) ----------------------------------
  if (route == "spatial") {
    
    # Extract swath components
    bp_result <- run_branching_process(
      N_pop        = swath$cell_pop_total,
      hub_mask     = swath$hub_cell_mask,
      kernel       = swath$kernel_matrix,
      spill_idx    = swath$spill_cell_idx,
      R0           = R0,
      k            = k_val,
      N_crit       = N_crit,
      g_max        = g_max,
      capture_frac = swath$capture_fraction,
      log_gens     = log_gens,
      log_cells    = log_cells
    )
  }
  
  # --- Non-spatial at hub (Stage 3c) -----------------------------------------
  if (route == "nonspatial_hub") {
    
    # Determine hub population: use swath metadata if available
    hub_pop <- nonspatial_pop  # fallback
    
    if (!is.null(hub_info) && !is.null(swath$hub_id)) {
      hub_row <- hub_info[hub_info$hub_id == swath$hub_id, ]
      if (nrow(hub_row) > 0) {
        # hub_info$pop_k is in thousands
        hub_pop <- hub_row$pop_k[1] * 1000
      }
    }
    
    # Use a generous population: hub populations are large enough that 
    # susceptible depletion is negligible in the early phase
    hub_pop <- max(hub_pop, 50000)
    
    ns <- make_nonspatial_swath(pop = hub_pop, is_hub = TRUE)
    
    bp_result <- run_branching_process(
      N_pop        = ns$N_pop,
      hub_mask     = ns$hub_mask,
      kernel       = ns$kernel,
      spill_idx    = ns$spill_idx,
      R0           = R0,
      k            = k_val,
      N_crit       = N_crit,
      g_max        = g_max,
      capture_frac = ns$capture_frac,
      log_gens     = log_gens,
      log_cells    = log_cells
    )
  }
  
  # --- Package result --------------------------------------------------------
  result <- list(
    event_id       = event$event_id,
    route          = route,
    R0             = R0,
    k              = k_val,
    outcome        = bp_result$outcome,
    chain_size     = bp_result$chain_size,
    hub_infections = bp_result$hub_infections,
    generations    = bp_result$generations,
    peak_infected  = bp_result$peak_infected
  )
  
  if (log_gens && !is.null(bp_result$gen_log)) {
    result$gen_log <- bp_result$gen_log
  }
  
  if (log_cells && !is.null(bp_result$cell_log)) {
    result$cell_log <- bp_result$cell_log
  }
  
  return(result)
}
