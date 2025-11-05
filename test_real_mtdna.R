# ============================================================================
# test_real_mtdna.R: Testing De Bruijn Assembly on Real Mitochondrial DNA
# ============================================================================
# PURPOSE: Test the complete assembly pipeline on real DNA sequences
# WHAT: Load mtDNA, test multiple k values, record metrics
# WHY: Understand algorithm performance on realistic data

library(igraph)

# ============================================================================
# Load Real mtDNA Sequence
# ============================================================================

cat("\n========== LOADING REAL MTDNA SEQUENCE ==========\n\n")

fasta_file <- "data/mtDNA_fragment.fasta"

# Read FASTA file
fasta_lines <- readLines(fasta_file)
header_idx <- grep("^>", fasta_lines)
sequence_lines <- fasta_lines[-header_idx]
mtdna_sequence <- paste(sequence_lines, collapse = "")
mtdna_sequence <- toupper(mtdna_sequence)

cat("Sequence loaded successfully!\n")
cat("Sequence length:", nchar(mtdna_sequence), "bp\n")
cat("First 100 bp:", substr(mtdna_sequence, 1, 100), "\n\n")

# ============================================================================
# Test Multiple K Values on Real mtDNA
# ============================================================================

k_values <- c(5, 7, 10, 15, 20, 25)
results <- data.frame(
  k = integer(),
  num_kmers = integer(),
  num_nodes = integer(),
  num_edges = integer(),
  avg_degree = numeric(),
  balanced_nodes = integer(),
  graph_connected = logical(),
  stringsAsFactors = FALSE
)

cat("========== TESTING K VALUES ON MTDNA ==========\n\n")

for (k in k_values) {
  cat(sprintf("Testing K = %d...\n", k))
  
  seq_len <- nchar(mtdna_sequence)
  num_kmers <- seq_len - k + 1
  
  if (num_kmers <= 0) {
    cat(sprintf("  ✗ K is too large (need at least %d)\n\n", k))
    next
  }
  
  # Extract k-mers
  kmers <- character(num_kmers)
  for (i in 1:num_kmers) {
    kmers[i] <- substr(mtdna_sequence, i, i + k - 1)
  }
  
  # Remove duplicates for graph building (unique k-mers)
  unique_kmers <- unique(kmers)
  
  # Build edge list
  edge_list <- data.frame(
    from = substr(unique_kmers, 1, k-1),
    to = substr(unique_kmers, 2, k),
    stringsAsFactors = FALSE
  )
  
  # Get unique nodes
  all_nodes <- unique(c(edge_list$from, edge_list$to))
  num_nodes <- length(all_nodes)
  num_edges <- nrow(edge_list)
  
  # Build graph with igraph
  g <- graph_from_data_frame(edge_list, directed = TRUE)
  
  # Calculate degrees
  in_deg <- degree(g, mode = "in")
  out_deg <- degree(g, mode = "out")
  
  # Count balanced nodes
  degree_balance <- out_deg - in_deg
  balanced_nodes <- sum(degree_balance == 0)
  
  # Check connectivity
  is_connected <- is_connected(g, mode = "weak")
  num_components <- components(g, mode = "weak")$no
  
  # Average degree
  avg_degree <- mean(c(in_deg, out_deg))
  
  # Print results
  cat(sprintf("  K-mers extracted: %d (unique: %d)\n", num_kmers, length(unique_kmers)))
  cat(sprintf("  Nodes: %d\n", num_nodes))
  cat(sprintf("  Edges: %d\n", num_edges))
  cat(sprintf("  Avg degree: %.2f\n", avg_degree))
  cat(sprintf("  Balanced nodes: %d/%d\n", balanced_nodes, num_nodes))
  cat(sprintf("  Connected: %s (%d component%s)\n", 
              if(is_connected) "YES" else "NO", 
              num_components,
              if(num_components > 1) "s" else ""))
  
  # Check if Eulerian path exists
  imbalanced <- which(degree_balance != 0)
  sources <- which(out_deg - in_deg == 1)
  sinks <- which(in_deg - out_deg == 1)
  eulerian_path_exists <- (length(sources) == 0 && length(sinks) == 0) ||
                          (length(sources) == 1 && length(sinks) == 1)
  
  cat(sprintf("  Eulerian path possible: %s\n", if(eulerian_path_exists) "YES ✓" else "NO ✗"))
  cat("\n")
  
  # Store results
  results <- rbind(results, data.frame(
    k = k,
    num_kmers = num_kmers,
    num_nodes = num_nodes,
    num_edges = num_edges,
    avg_degree = avg_degree,
    balanced_nodes = balanced_nodes,
    graph_connected = is_connected
  ))
}

# ============================================================================
# Display Summary Table
# ============================================================================

cat("\n========== SUMMARY TABLE ==========\n\n")
print(results)

# ============================================================================
# Analysis and Recommendations
# ============================================================================

cat("\n========== ANALYSIS & RECOMMENDATIONS ==========\n\n")

cat("Observations:\n")
cat("1. Smaller k values (5-10):\n")
cat("   - More k-mers → larger graphs\n")
cat("   - More likely to be connected\n")
cat("   - Risk of false overlaps\n\n")

cat("2. Larger k values (15+):\n")
cat("   - Fewer k-mers → smaller graphs\n")
cat("   - Risk of disconnected graph\n")
cat("   - Better specificity\n\n")

cat("3. For this mtDNA sequence:\n")

# Find best k (most balanced)
best_balance_idx <- which.max(results$balanced_nodes)
if(length(best_balance_idx) > 0) {
  best_k <- results$k[best_balance_idx]
  cat("   - Best k value:", best_k, "\n")
  cat("   - Balanced nodes:", results$balanced_nodes[best_balance_idx], "/", 
      results$num_nodes[best_balance_idx], "\n")
}

# Find connected graphs
connected_results <- results[results$graph_connected, ]
if(nrow(connected_results) > 0) {
  cat("   - K values with connected graphs:", paste(connected_results$k, collapse=", "), "\n")
} else {
  cat("   - No connected graphs at tested k values\n")
}

cat("\nTEST COMPLETE \n")