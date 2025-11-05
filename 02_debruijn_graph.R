# 
# 02_debruijn_graph.R: De Bruijn Graph Construction
# 
# PURPOSE: Build a directed multigraph where nodes are (k-1)-mers and edges
#          are k-mers, representing overlaps.
# WHAT: Create nodes and edges; use igraph to represent the graph.
# WHY: The De Bruijn graph encodes all overlapping information. An Eulerian
#      path through this graph reconstructs the original sequence.
#

library(igraph)
# 
# STEP 1: Verify k-mer data is loaded
# 
# WHAT: Check that the previous script (01_kmers.R) was sourced successfully.
# WHY: We need kmers, prefixes, and suffixes to build the graph.

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
# 
# WHAT: For each k-mer, create a directed edge from its prefix to its suffix.
# WHY: Each edge represents a k-mer, and the edge direction represents the
#      overlap information. Prefix -> Suffix encodes "these (k-1)-mers overlap
#      in this k-mer."

edge_list <- data.frame(
  from = prefixes,    # Source node: (k-1)-mer prefix
  to = suffixes,      # Target node: (k-1)-mer suffix
  stringsAsFactors = FALSE
)

cat("\nEdge List (first 10 rows):\n")
print(head(edge_list, 10))

# 
# STEP 3: Build igraph object
# 
# WHAT: Convert the edge list into an igraph directed multigraph.
# WHY: igraph provides efficient graph operations, degree calculations,
#      and connectivity checks needed for Eulerian path computation.

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
# 
# WHAT: For each node, count incoming and outgoing edges.
# WHY: Eulerian path/circuit existence depends on degree balance.
#      - Circuit: all nodes have in-degree = out-degree
#      - Path: exactly 0 or 2 nodes with imbalanced degrees

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
# 
# WHAT: Determine if the graph is weakly connected (all nodes reachable via
#       undirected paths).
# WHY: For an Eulerian path to exist, all nodes with edges must be in one
#      connected component. Multiple components mean multiple separate contigs.

weakly_conn <- is_connected(dbg, mode = "weak")
num_components <- components(dbg, mode = "weak")$no

cat("\nConnectivity Analysis:\n")
cat("Weakly Connected:", weakly_conn, "\n")
cat("Number of Weak Components:", num_components, "\n")

# STEP 6: Prepare output for next script
# WHAT: Store the graph and degree information for Eulerian path computation.
# WHY: The next script (03_eulerian_path.R) needs the graph structure and
#      degree information to run Hierholzer's algorithm.

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
