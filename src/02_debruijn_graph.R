

library(igraph)
# STEP 1: Verify k-mer data is loaded

if (!exists("kmer_data")) {
  stop("K-mer data not found. Run 01_kmers.R first.")
}

kmers <- kmer_data$kmers
prefixes <- kmer_data$prefixes
suffixes <- kmer_data$suffixes
k <- kmer_data$k

cat("Building De Bruijn graph with k =", k, "\n")

# 
# STEP 2: Create edge list

edge_list <- data.frame(
  from = prefixes,    # Source node: (k-1)-mer prefix
  to = suffixes,      # Target node: (k-1)-mer suffix
  stringsAsFactors = FALSE
)

cat("\nEdge List (first 10 rows):\n")
print(head(edge_list, 10))

# 
# STEP 3: Build igraph object

dbg <- graph_from_data_frame(
  d = edge_list,
  directed = TRUE,
  vertices = NULL  # Auto-detect nodes from edge list
)

cat("\nDe Bruijn Graph Created:\n")
cat("Number of Nodes:", vcount(dbg), "\n")
cat("Number of Edges:", ecount(dbg), "\n")

# 
# STEP 4: Compute in-degree and out-degree for each node

in_degrees <- degree(dbg, mode = "in")
out_degrees <- degree(dbg, mode = "out")

degree_table <- data.frame(
  node = names(in_degrees),
  in_degree = as.numeric(in_degrees),
  out_degree = as.numeric(out_degrees),
  balance = as.numeric(out_degrees) - as.numeric(in_degrees),
  stringsAsFactors = FALSE
)

cat("\nDegree Analysis:\n")
print(degree_table)

# STEP 5: Check graph connectivity

weakly_conn <- is_connected(dbg, mode = "weak")
num_components <- components(dbg, mode = "weak")$no

cat("\nConnectivity Analysis:\n")
cat("Weakly Connected:", weakly_conn, "\n")
cat("Number of Weak Components:", num_components, "\n")

# STEP 6: Prepare output for next script
dbg_data <- list(
  graph = dbg,
  edge_list = edge_list,
  in_degrees = in_degrees,
  out_degrees = out_degrees,
  degree_table = degree_table,
  weakly_connected = weakly_conn,
  num_components = num_components,
  kmers = kmers,
  k = k
)

cat("\n= DE BRUIJN GRAPH COMPLETE =\n")
cat("Graph ready for Eulerian path computation.\n")
