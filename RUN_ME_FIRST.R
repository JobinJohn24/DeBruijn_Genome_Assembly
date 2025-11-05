# RUN_ME_FIRST.R: Complete De Bruijn Assembly Workflow
# PURPOSE: Execute the entire De Bruijn graph assembly pipeline in sequence.
# WHAT: Source all scripts in the correct order to go from sequence → metrics.
# WHY: This file demonstrates the complete workflow and validates that all
#      components work together correctly.
# 

cat("╔n")
cat("║De Bruijn Graph Genome Assembly - Complete Workflow            ║\n")
cat("╚═\n\n")


# Check for required packages
# WHAT: Verify that igraph, stringr, stringdist are installed.
# WHY: These packages are required by the scripts; missing packages cause errors.

required_packages <- c("igraph", "stringr", "stringdist")
missing_packages <- required_packages[!sapply(required_packages, 
                                               function(pkg) require(pkg, quietly = TRUE))]

if (length(missing_packages) > 0) {
  cat("Installing missing packages...\n")
  install.packages(missing_packages)
}

cat("✓ All packages loaded successfully.\n\n")

# ============================================================================
# PHASE 1: Load Toy Data
# ============================================================================
# WHAT: Initialize the toy dataset that will be used throughout.
# WHY: All subsequent scripts depend on having a sequence to process.

cat("────────────────────────────────────────────────────────────────────\n")
cat("PHASE 1: Loading toy data...\n")
cat("────────────────────────────────────────────────────────────────────\n")

source("toy_data.R")

cat("\n")
# PHASE 2: Extract K-mers
# WHAT: Break the sequence into overlapping k-mers.
# WHY: K-mers are the input to the De Bruijn graph; this generates them.

cat("────────────────────────────────────────────────────────────────────\n")
cat("PHASE 2: Extracting k-mers from sequence...\n")
cat("────────────────────────────────────────────────────────────────────\n")

source("01_kmers.R")

cat("\n")


# PHASE 3: Build De Bruijn Graph
# WHAT: Construct a directed multigraph from the k-mers.
# WHY: The graph encodes all overlapping sequence information.

cat("────────────────────────────────────────────────────────────────────\n")
cat("PHASE 3: Building De Bruijn graph...\n")
cat("────────────────────────────────────────────────────────────────────\n")

source("02_debruijn_graph.R")

cat("\n")

# PHASE 4: Compute Eulerian Path
# WHAT: Use Hierholzer's algorithm to find an Eulerian path.
# WHY: An Eulerian path traverses every edge exactly once; its sequence
#      order yields the assembled sequence.

cat("────────────────────────────────────────────────────────────────────\n")
cat("PHASE 4: Finding Eulerian path (Hierholzer algorithm)...\n")
cat("────────────────────────────────────────────────────────────────────\n")

source("03_eulerian_path.R")

cat("\n")


# PHASE 5: Reconstruct Sequence
# WHAT: Convert the Eulerian node path into a DNA sequence.
# WHY: The path order encodes the sequence; walking it reconstructs it.

cat("────────────────────────────────────────────────────────────────────\n")
cat("PHASE 5: Reconstructing sequence from Eulerian path...\n")
cat("────────────────────────────────────────────────────────────────────\n")

source("04_reconstruct_sequence.R")

cat("\n")

# PHASE 6: Compute Metrics
# WHAT: Evaluate assembly quality using standard bioinformatics metrics.
# WHY: Metrics quantify how well the assembly succeeded.

cat("\n")
cat("PHASE 6: Computing assembly metrics...\n")
cat("\n")

source("05_metrics.R")

cat("\n")

# 
# WORKFLOW COMPLETE
# 
# WHAT: Summarize what was accomplished.
# WHY: Confirms successful execution and summarizes key results.

cat("╔════════════════════════════════════════════════════════════════════╗\n")
cat("║                 WORKFLOW COMPLETE ✓                                ║\n")
cat("╚════════════════════════════════════════════════════════════════════╝\n\n")

cat("Summary of Results:\n")
cat("  Original Sequence:     ", recon_data$original_sequence, "\n")
cat("  Reconstructed Seq:     ", recon_data$reconstructed_sequence, "\n")
cat("  Match:                 ", recon_data$sequences_match, "\n")
cat("  Percent Identity:      ", 
    round(metrics$accuracy$pct_identity, 2), "%\n", sep = "")
cat("  Graph Nodes/Edges:     ", metrics$graph$num_nodes, "/", 
    metrics$graph$num_edges, "\n", sep = "")
cat("  Eulerian Path Type:    ", metrics$path$type, "\n")
cat("\nAll data objects are now in global environment:\n")
cat("  - kmer_data        (k-mer extraction results)\n")
cat("  - dbg_data         (De Bruijn graph and degree info)\n")
cat("  - euler_data       (Eulerian path computed)\n")
cat("  - recon_data       (reconstructed sequence)\n")
cat("  - metrics          (quality metrics)\n")
cat("\nNext Steps:\n")
cat("  1. Experiment with different k values (in 01_kmers.R)\n")
cat("  2. Try other toy sequences (in toy_data.R)\n")
cat("  3. Add noise and test error handling\n")
cat("  4. Visualize the graph using igraph plotting functions\n")
cat("\nNext Steps:\n")
cat("  1. Experiment with different k values (in 01_kmers.R)\n")
cat("  2. Try other toy sequences (in toy_data.R)\n")
cat("  3. Add noise and test error handling\n")
cat("  4. Visualize the graph using igraph plotting functions\n")

# Visualize the Graph
cat("\n────────────────────────────────────────────────────────────────────\n")
cat("Generating graph visualization...\n")
cat("────────────────────────────────────────────────────────────────────\n\n")

source("visualize_graph.R")

cat("\n")

# 
# END OF WORKFLOW
# 

cat("╔════════════════════════════════════════════════════════════════════╗\n")
cat("║                 ALL STEPS COMPLETE! ✓✓✓                            ║\n")
cat("╚════════════════════════════════════════════════════════════════════╝\n\n")

cat("Happy exploring! 🧬\n")

