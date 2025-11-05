# ============================================================================
# toy_data.R: Toy Dataset Initialization
# ============================================================================
# PURPOSE: Provide pre-defined toy sequences for testing and learning.
# WHAT: Define example DNA sequences to use throughout the pipeline.
# WHY: These sequences are small enough to trace by hand, yet complex enough
#      to demonstrate the full De Bruijn graph and Eulerian path algorithms.
# ============================================================================

# ============================================================================
# OPTION 1: Simple cyclic sequence (good for visualization)
# ============================================================================
# WHAT: A short, repeating sequence that forms an Eulerian circuit.
# WHY: Cycles are easier to understand; Eulerian circuit (not path) is a
#      special case worth studying.

toy_cyclic <- "ACGTAC"

# ============================================================================
# OPTION 2: Linear sequence (good for path understanding)
# ============================================================================
# WHAT: A sequence that will form an Eulerian path (start ≠ end).
# WHY: Linear paths are the most common case in real genome assembly.

toy_linear <- "ACGTACGTAC"

# ============================================================================
# OPTION 3: Longer sequence (more complex graph)
# ============================================================================
# WHAT: A longer sequence with more k-mers and higher graph complexity.
# WHY: Tests scalability and demonstrates how graphs grow with sequence length.

toy_complex <- "ACGTACGTACGTACGT"

# ============================================================================
# Select active toy sequence (modify this to switch)
# ============================================================================
# WHAT: Choose which toy sequence to use by default.
# WHY: Allows easy experimentation without editing all downstream scripts.

ACTIVE_TOY <- toy_cyclic

cat("Toy Dataset Loaded:\n")
cat("Available sequences:\n")
cat("  - toy_cyclic (length ", nchar(toy_cyclic), "): ", toy_cyclic, "\n", sep = "")
cat("  - toy_linear (length ", nchar(toy_linear), "): ", toy_linear, "\n", sep = "")
cat("  - toy_complex (length ", nchar(toy_complex), "): ", toy_complex, "\n", sep = "")
cat("\nActive toy sequence: toy_cyclic\n")
cat("To switch, edit ACTIVE_TOY variable or modify 01_kmers.R\n")
