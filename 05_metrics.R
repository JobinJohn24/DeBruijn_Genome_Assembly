
# 05_metrics.R: Assembly and Graph Metrics
# 
# PURPOSE: Evaluate assembly quality using standard bioinformatics metrics.
# WHAT: Calculate N50, assembly length, identity, and graph statistics.
# WHY: Metrics quantify assembly quality and help identify problems.
#      They're important for comparing different assembly parameters.


library(stringr)
library(stringdist)
library(igraph)

# 
# STEP 1: Verify all data is loaded
# 
# WHAT: Check that all previous scripts were sourced successfully.
# WHY: We need reconstruction results, De Bruijn graph, and k-mer data.

if (!exists("recon_data") || !exists("dbg_data")) {
  stop("Missing data. Run all scripts 01-04 first.")
}

original <- recon_data$original_sequence
reconstructed <- recon_data$reconstructed_sequence
dbg <- dbg_data$graph
k <- recon_data$k

cat("=ASSEMBLY QUALITY METRICS =\n")


# STEP 2: Basic Assembly Metrics
# WHAT: Calculate total length, contig count, and N50.
# WHY: These quantify completeness and contiguity. N50 is standard in genomics:
#      50% of bases are in contigs >= N50 length (higher is better).

assembly_length <- nchar(reconstructed)
num_contigs <- 1  # For perfect case, one contig

cat("\nAssembly Statistics:\n")
cat("  Assembly Length:", assembly_length, "bp\n")
cat("  Number of Contigs:", num_contigs, "\n")

# WHAT: N50 for single contig is simply the contig length.
# WHY: In multi-contig assemblies, N50 is computed differently; here it's trivial.
n50 <- assembly_length

cat("  N50:", n50, "bp\n")

# STEP 3: Sequence Identity & Error Rate
# WHAT: Compare reconstructed vs. original base-by-base.
# WHY: Identity measures accuracy. Perfect assembly = 100% identity.

if (recon_data$sequences_match) {
  # WHAT: If sequences match exactly, identity is 100%.
  # WHY: No mismatches, insertions, or deletions.
  pct_identity <- 100.0
  num_mismatches <- 0
} else {
  # WHAT: Compute base-by-base comparison; count matches.
  # WHY: Even with mismatches, we can measure partial identity.
  min_len <- min(nchar(original), nchar(reconstructed))
  matches <- 0
  
  for (i in 1:min_len) {
    if (str_sub(original, i, i) == str_sub(reconstructed, i, i)) {
      matches <- matches + 1
    }
  }
  
  # WHAT: Identity = (matches / reference length) * 100.
  # WHY: Standard metric; always use reference (original) as denominator.
  pct_identity <- (matches / nchar(original)) * 100
  num_mismatches <- nchar(original) - matches
}

cat("\nSequence Accuracy:\n")
cat("  Percent Identity:", round(pct_identity, 2), "%\n")
cat("  Mismatches:", num_mismatches, "\n")

# STEP 4: Edit Distance (Levenshtein)
# WHAT: Calculate minimum edits (insertions, deletions, substitutions) needed
#       to transform reconstructed into original.
# WHY: Edit distance quantifies overall difference. Lower is better.

edit_distance <- stringdist(original, reconstructed, method = "lv")

cat("  Levenshtein Distance:", edit_distance, "\n")

# STEP 5: De Bruijn Graph Metrics
# WHAT: Calculate structural properties of the De Bruijn graph.
# WHY: Graph metrics reveal assembly complexity, errors, and redundancy.

num_nodes <- vcount(dbg)
num_edges <- ecount(dbg)
avg_out_degree <- mean(degree(dbg, mode = "out"))
avg_in_degree <- mean(degree(dbg, mode = "in"))

cat("\nDe Bruijn Graph Metrics:\n")
cat("  Number of Nodes:", num_nodes, "\n")
cat("  Number of Edges:", num_edges, "\n")
cat("  Average Out-Degree:", round(avg_out_degree, 2), "\n")
cat("  Average In-Degree:", round(avg_in_degree, 2), "\n")

# WHAT: Count "tips" (dead-end nodes with in or out degree = 1, but not both).
# WHY: Tips suggest errors or incomplete data; they're common in real assemblies.
in_deg <- degree(dbg, mode = "in")
out_deg <- degree(dbg, mode = "out")
tips <- sum((in_deg == 0 & out_deg > 0) | (in_deg > 0 & out_deg == 0))

cat("  Number of Tips (dead ends):", tips, "\n")

# STEP 6: Eulerian Path Properties
# WHAT: Summarize the Eulerian path that was computed.
# WHY: Confirms the path type and validates assembly feasibility.

cat("\nEulerian Path Properties:\n")
cat("  Path Type:", recon_data$path_type, "\n")
cat("  Start Node:", recon_data$start_node, "\n")
cat("  End Node:", recon_data$end_node, "\n")
cat("  Path Length (edges):", length(recon_data$path) - 1, "\n")

# STEP 7: Final Summary Report
# WHAT: Produce a concise summary of all metrics.
# WHY: Provides quick overview of assembly quality; useful for diagnostics.

cat("\n========== FINAL ASSEMBLY REPORT ==========\n")
cat("Original Sequence Length:", nchar(original), "bp\n")
cat("Assembled Sequence Length:", nchar(reconstructed), "bp\n")
cat("Percent Identity: ", round(pct_identity, 2), "%\n", sep = "")
cat("N50 (contiguity):", n50, "bp\n")
cat("Graph Complexity (nodes/edges):", num_nodes, "/", num_edges, "\n")

if (pct_identity == 100) {
  cat("\n✓ ASSEMBLY SUCCESSFUL: Perfect reconstruction achieved!\n")
} else {
  cat("\n⚠ Imperfect assembly. Consider:\n")
  cat("  - Increasing k-mer length for accuracy\n")
  cat("  - Adjusting error filtering threshold\n")
  cat("  - Checking for incomplete coverage\n")
}

# STEP 8: Save results for reference
#structured list for later access/export.
# WHY: Facilitates downstream analysis, visualization, or report generation.

metrics <- list(
  assembly = list(
    length = assembly_length,
    n50 = n50,
    num_contigs = num_contigs
  ),
  accuracy = list(
    pct_identity = pct_identity,
    num_mismatches = num_mismatches,
    edit_distance = edit_distance
  ),
  graph = list(
    num_nodes = num_nodes,
    num_edges = num_edges,
    avg_out_degree = avg_out_degree,
    avg_in_degree = avg_in_degree,
    tips = tips
  ),
  path = list(
    type = recon_data$path_type,
    start = recon_data$start_node,
    end = recon_data$end_node,
    length = length(recon_data$path) - 1
  )
)

cat("\n========== METRICS CALCULATION COMPLETE ==========\n")
cat("All assembly metrics computed and available in 'metrics' object.\n")
