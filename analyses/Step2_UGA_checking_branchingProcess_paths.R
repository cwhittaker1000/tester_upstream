# =============================================================================
# trace_hub_paths.R
#
# Runs 1000 spatial branching process simulations on the 5-cell test corridor
# and classifies the transmission path for each hub_established outcome.
#
# Path classification:
#   "direct"    — hub (cell 5) was infected before any intermediate cell (2,3,4)
#   "relay"     — at least one intermediate cell was infected before the hub
#   "simultaneous" — hub and an intermediate cell infected in the same generation
#
# For relay paths, we also record the full sequence of first-infection generations.
#
# Dependencies: source("functions/branching_process.R")
# =============================================================================

source("functions/branching_process.R")

# --- Build the same 5-cell test corridor as in Step2 -------------------------
n_test <- 5
pops_test <- c(500, 200, 300, 200, 50000)   # cell 1 = spillover, cell 5 = hub
hub_test  <- c(FALSE, FALSE, FALSE, FALSE, TRUE)

# Distance matrix (linear, 10km spacing)
coords_test <- cbind(seq(0, 40, by = 10), rep(0, 5))
D_test <- as.matrix(dist(coords_test))

# Gravity kernel (same parameters as Step2)
alpha <- 1.0; beta <- 0.1; w_self <- 100
pop_attract_test <- pops_test^alpha
W_test <- sweep(exp(-beta * D_test), 2, pop_attract_test, FUN = "*")
diag(W_test) <- diag(W_test) + w_self
K_test <- W_test / rowSums(W_test)

cf_test <- rep(1.0, n_test)

# --- Print the kernel to understand expected behaviour -----------------------
cat("=== Gravity kernel (row = source, col = destination) ===\n")
cat("Columns: [Spill(500)] [Rural(200)] [Rural(300)] [Rural(200)] [Hub(50k)]\n\n")
print(round(K_test, 4))

cat("\nKey insight: From cell 1, P(offspring → cell 5) =", 
    round(K_test[1, 5], 4), "\n")
cat("  vs P(offspring → cell 2) =", round(K_test[1, 2], 4), "\n")
cat("  vs P(offspring → cell 3) =", round(K_test[1, 3], 4), "\n")
cat("  vs P(offspring → cell 4) =", round(K_test[1, 4], 4), "\n")
cat("  vs P(offspring → cell 1) =", round(K_test[1, 1], 4), "(self)\n\n")

# =============================================================================
# classify_path(): Examine a cell_log matrix and classify the route to hub
# =============================================================================
#
# For each cell, find the first generation it received any infection.
# Then compare the hub's first-infection generation to intermediates.
#
# Returns a list with:
#   $path_type     "direct", "relay", or "simultaneous"
#   $first_gen     Named vector: generation of first infection per cell (NA if never)
#   $hub_first_gen Generation when hub was first infected
#   $first_intermediate  Which intermediate cell (2,3,4) was infected first (NA if direct)
#   $sequence      Character string describing the full ordering, e.g. "1→5" or "1→3→5"

classify_path <- function(cell_log) {
  
  n_cells <- ncol(cell_log)
  n_gens  <- nrow(cell_log)
  
  # Find the first generation each cell received infections (NA if never)
  first_gen <- rep(NA_integer_, n_cells)
  for (j in 1:n_cells) {
    infected_gens <- which(cell_log[, j] > 0)
    if (length(infected_gens) > 0) {
      first_gen[j] <- infected_gens[1] - 1L  # convert row index to generation number
      # (row 1 = generation 0, row 2 = generation 1, etc.)
    }
  }
  names(first_gen) <- paste0("cell", 1:n_cells)
  
  hub_first <- first_gen[5]         # hub = cell 5
  intermediate_firsts <- first_gen[2:4]  # cells 2, 3, 4
  
  # If hub was never infected, this shouldn't be a hub_established outcome
  # but handle gracefully
  if (is.na(hub_first)) {
    return(list(
      path_type          = "no_hub_infection",
      first_gen          = first_gen,
      hub_first_gen      = NA,
      first_intermediate = NA,
      sequence           = "1→extinct"
    ))
  }
  
  # Were any intermediates infected before or at the same time as the hub?
  intermediates_before <- which(!is.na(intermediate_firsts) & 
                                  intermediate_firsts < hub_first)
  intermediates_same   <- which(!is.na(intermediate_firsts) & 
                                  intermediate_firsts == hub_first)
  
  # --- Build the full sequence string ----------------------------------------
  # Collect all cells that were ever infected, ordered by first_gen
  infected_cells <- which(!is.na(first_gen))
  infection_order <- infected_cells[order(first_gen[infected_cells])]
  
  # Group by generation for the sequence string
  gens_present <- sort(unique(first_gen[infected_cells]))
  seq_parts <- c()
  for (gen in gens_present) {
    cells_this_gen <- infected_cells[first_gen[infected_cells] == gen]
    # Label hub cell specially
    labels <- ifelse(cells_this_gen == 5, "5(hub)", as.character(cells_this_gen))
    if (length(labels) > 1) {
      seq_parts <- c(seq_parts, paste0("[", paste(labels, collapse = "+"), "]"))
    } else {
      seq_parts <- c(seq_parts, labels)
    }
  }
  sequence_str <- paste(seq_parts, collapse = "→")
  
  # --- Classify --------------------------------------------------------------
  if (length(intermediates_before) > 0) {
    # At least one intermediate was infected strictly before the hub
    first_inter_cell <- intermediates_before[which.min(
      intermediate_firsts[intermediates_before])] + 1  # +1 to get cell index (2,3,4)
    path_type <- "relay"
    
  } else if (length(intermediates_same) > 0) {
    # An intermediate was infected in the same generation as the hub
    first_inter_cell <- intermediates_same[1] + 1
    path_type <- "simultaneous"
    
  } else {
    # Hub was infected before any intermediate (or no intermediates infected)
    first_inter_cell <- NA
    path_type <- "direct"
  }
  
  return(list(
    path_type          = path_type,
    first_gen          = first_gen,
    hub_first_gen      = hub_first,
    first_intermediate = first_inter_cell,
    sequence           = sequence_str
  ))
}

# =============================================================================
# Run 1000 simulations with cell logging
# =============================================================================

set.seed(123)
n_sims <- 1000

cat(sprintf("Running %d spatial simulations (R0=1.5, k=1.0, N_crit=100)...\n\n", n_sims))

results <- vector("list", n_sims)

for (i in 1:n_sims) {
  res <- run_branching_process(
    N_pop = pops_test, hub_mask = hub_test, kernel = K_test, spill_idx = 1L,
    R0 = 1.5, k = 1.0, N_crit = 100, capture_frac = cf_test,
    log_cells = TRUE
  )
  
  # Classify path if hub was established
  path_info <- classify_path(res$cell_log)
  
  results[[i]] <- list(
    outcome    = res$outcome,
    chain_size = res$chain_size,
    hub_inf    = res$hub_infections,
    gens       = res$generations,
    path_type  = path_info$path_type,
    first_gen  = path_info$first_gen,
    hub_first  = path_info$hub_first_gen,
    sequence   = path_info$sequence
  )
  
  if (i %% 200 == 0) cat(sprintf("  %d / %d done\n", i, n_sims))
}

# =============================================================================
# Analyse results
# =============================================================================

outcomes  <- sapply(results, `[[`, "outcome")
paths     <- sapply(results, `[[`, "path_type")
sequences <- sapply(results, `[[`, "sequence")

cat("\n=== Overall Outcomes ===\n")
print(table(outcomes))

cat(sprintf("\n  P(extinct):         %.3f\n", mean(outcomes == "extinct")))
cat(sprintf("  P(hub_established): %.3f\n", mean(outcomes == "hub_established")))
cat(sprintf("  P(generation_cap):  %.3f\n", mean(outcomes == "generation_cap")))

# --- Path analysis for hub_established events --------------------------------
hub_est <- which(outcomes == "hub_established")

cat(sprintf("\n=== Path Classification (%d hub-established events) ===\n", length(hub_est)))
path_table <- table(paths[hub_est])
print(path_table)
cat("\n")
for (p in names(path_table)) {
  cat(sprintf("  %-15s  %4d  (%5.1f%% of established)\n",
              p, path_table[p], 100 * path_table[p] / length(hub_est)))
}

# --- Most common transmission sequences --------------------------------------
cat("\n=== Most Common Transmission Sequences (hub-established) ===\n")
seq_table <- sort(table(sequences[hub_est]), decreasing = TRUE)
n_show <- min(15, length(seq_table))
for (i in 1:n_show) {
  cat(sprintf("  %4d  (%5.1f%%)  %s\n",
              seq_table[i],
              100 * seq_table[i] / length(hub_est),
              names(seq_table)[i]))
}
if (length(seq_table) > n_show) {
  cat(sprintf("  ... and %d more unique sequences\n", length(seq_table) - n_show))
}

# --- Generation of first hub infection ---------------------------------------
hub_first_gens <- sapply(results[hub_est], `[[`, "hub_first")

cat("\n=== Generation of First Hub Infection (hub-established events) ===\n")
cat(sprintf("  Median: generation %d\n", median(hub_first_gens)))
cat(sprintf("  Range:  [%d, %d]\n", min(hub_first_gens), max(hub_first_gens)))
cat("  Distribution:\n")
gen_table <- table(hub_first_gens)
for (g in names(gen_table)) {
  cat(sprintf("    Gen %2s: %4d  (%5.1f%%)\n",
              g, gen_table[g], 100 * gen_table[g] / length(hub_est)))
}

# --- Path analysis for EXTINCT events (what happened before extinction) ------
extinct_idx <- which(outcomes == "extinct")

cat(sprintf("\n=== Path info for extinct events (%d total) ===\n", length(extinct_idx)))
cat("Did the hub ever get infected before the chain went extinct?\n")

hub_ever_infected <- sapply(results[extinct_idx], function(r) !is.na(r$hub_first))
cat(sprintf("  Hub reached but chain extinct: %d / %d (%.1f%%)\n",
            sum(hub_ever_infected), length(extinct_idx),
            100 * mean(hub_ever_infected)))
cat(sprintf("  Hub never reached:             %d / %d (%.1f%%)\n",
            sum(!hub_ever_infected), length(extinct_idx),
            100 * mean(!hub_ever_infected)))

# --- Example traces ----------------------------------------------------------
cat("\n=== Example Traces ===\n")

# Show a few direct paths
direct_idx <- hub_est[paths[hub_est] == "direct"]
if (length(direct_idx) > 0) {
  cat("\n--- Direct paths (1→5 first): ---\n")
  for (i in head(direct_idx, 3)) {
    r <- results[[i]]
    cat(sprintf("  Sim %d: %s  |  hub at gen %d  |  chain=%d  |  gens=%d\n",
                i, r$sequence, r$hub_first, r$chain_size, r$gens))
  }
}

# Show a few relay paths
relay_idx <- hub_est[paths[hub_est] == "relay"]
if (length(relay_idx) > 0) {
  cat("\n--- Relay paths (intermediate first): ---\n")
  for (i in head(relay_idx, 5)) {
    r <- results[[i]]
    cat(sprintf("  Sim %d: %s  |  hub at gen %d  |  chain=%d  |  gens=%d\n",
                i, r$sequence, r$hub_first, r$chain_size, r$gens))
  }
}

# Show a few simultaneous paths
simul_idx <- hub_est[paths[hub_est] == "simultaneous"]
if (length(simul_idx) > 0) {
  cat("\n--- Simultaneous paths: ---\n")
  for (i in head(simul_idx, 3)) {
    r <- results[[i]]
    cat(sprintf("  Sim %d: %s  |  hub at gen %d  |  chain=%d  |  gens=%d\n",
                i, r$sequence, r$hub_first, r$chain_size, r$gens))
  }
}

cat("\n=== Analysis complete ===\n")