rm(list = ls())

k_values <- c(2, 3, 4, 5)

for (k in k_values) {
  cat("\n\n========== TESTING K =", k, "==========\n\n")
  
  sequence <- "ACGTAC"
  
  # Extract k-mers
  num_kmers <- nchar(sequence) - k + 1
  
  if (num_kmers <= 0) {
    cat("K is too large for this sequence! Skipping.\n")
    next
  }
  
  kmers <- character(num_kmers)
  for (i in 1:num_kmers) {
    kmers[i] <- substr(sequence, i, i + k - 1)
  }
  
  cat("K-mers extracted:", num_kmers, "\n")
  cat("K-mers:", paste(kmers, collapse=", "), "\n")
  
  # Build graph
  edge_list <- data.frame(
    from = substr(kmers, 1, k-1),
    to = substr(kmers, 2, k),
    stringsAsFactors = FALSE
  )
  
  cat("Nodes:", length(unique(c(edge_list$from, edge_list$to))), "\n")
  cat("Edges:", nrow(edge_list), "\n")
}