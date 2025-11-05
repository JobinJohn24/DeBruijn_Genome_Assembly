
# 03_eulerian_path.R: Hierholzer's Algorithm for Eulerian Path
# 
# PURPOSE: Find an Eulerian path through the De Bruijn graph that traverses
#          every edge exactly once.
# WHAT: Implement Hierholzer's algorithm using a depth-first search approach.
# WHY: This algorithm guarantees finding an Eulerian path in O(V + E) time if
#      one exists. We use it to determine the order of k-mers, which yields
#      the assembled sequence.

library(igraph)
# STEP 1: Verify graph data is loaded
# WHAT: Check that 02_debruijn_graph.R was sourced successfully.
# WHY: We need the graph, degrees, and connectivity info.

if (!exists("dbg_data")) {
  stop("De Bruijn graph not found. Run 02_debruijn_graph.R first.")
}

dbg <- dbg_data$graph
kmers <- dbg_data$kmers
k <- dbg_data$k

cat("Checking Eulerian Path Conditions...\n\n")

# 
# STEP 2: Extract node and degree information
# 
# WHAT: Get all nodes and their in/out degrees from the graph.
# WHY: We need this to check if an Eulerian path exists.

# Extract node names correctly
all_nodes <- names(V(dbg))

# Calculate degrees
in_deg <- degree(dbg, mode = "in")
out_deg <- degree(dbg, mode = "out")

cat("All nodes:", paste(all_nodes, collapse=", "), "\n")
cat("In-degrees:", paste(in_deg, collapse=", "), "\n")
cat("Out-degrees:", paste(out_deg, collapse=", "), "\n\n")

# 
# STEP 3: Check Eulerian path/circuit existence conditions
# 
# WHAT: Count nodes with degree imbalance; determine path type.
# WHY: - Eulerian circuit: all nodes balanced (in-deg = out-deg)
#      - Eulerian path: exactly 1 source and 1 sink
#      If conditions not met, no Eulerian path exists.

# Find imbalanced nodes
imbalanced_idx <- which(out_deg != in_deg)
source_idx <- which(out_deg - in_deg == 1)   # out-deg > in-deg
sink_idx <- which(in_deg - out_deg == 1)     # in-deg > out-deg

sources <- all_nodes[source_idx]
sinks <- all_nodes[sink_idx]

cat("Imbalanced nodes:", length(imbalanced_idx), "\n")
cat("Source nodes (out > in):", if(length(sources) > 0) paste(sources, collapse=", ") else "none", "\n")
cat("Sink nodes (in > out):", if(length(sinks) > 0) paste(sinks, collapse=", ") else "none", "\n\n")

# Check validity
path_exists <- (length(source_idx) == 0 && length(sink_idx) == 0) ||
               (length(source_idx) == 1 && length(sink_idx) == 1)

cat("Eulerian Path Exists:", path_exists, "\n\n")

if (!path_exists) {
  stop("Eulerian path conditions not satisfied. Cannot reconstruct sequence.")
}

# 
# STEP 4: Determine start node
# 
# WHAT: Choose where to begin traversal based on degree balance.
# WHY: - If there's a source (out-deg > in-deg), start there.
#      - If balanced (circuit), start at any node with edges.

start_node <- if (length(sources) > 0) sources[1] else all_nodes[1]
path_type <- if (length(source_idx) == 0) "Eulerian Circuit" else "Eulerian Path"

cat("Start Node:", start_node, "\n")
cat("Path Type:", path_type, "\n\n")

# STEP 5: Implement Hierholzer's Algorithm
# 
# WHAT: Use a stack-based DFS to traverse edges, ensuring each is visited once.
# WHY: Hierholzer ensures every edge is traversed exactly once, producing the
#      canonical Eulerian walk. The algorithm is efficient and handles cycles.

# Create edge list from graph
edge_list <- as_data_frame(dbg, what = "edges")
colnames(edge_list) <- c("from", "to")

cat("Edge List:\n")
print(edge_list)
cat("\n")

# Track which edges have been used (by row index)
edge_used <- rep(FALSE, nrow(edge_list))

# Stack for DFS traversal (LIFO - Last In First Out)
stack <- c(start_node)
path <- character(0)

cat("Executing Hierholzer's Algorithm...\n")

iteration <- 0
while (length(stack) > 0) {
  iteration <- iteration + 1
  
  # WHAT: Get current node from top of stack.
  # WHY: We process nodes in LIFO order, backtracking when stuck.
  v <- stack[length(stack)]
  
  # WHAT: Find an unused outgoing edge from current node.
  # WHY: If found, we follow it; if not, we backtrack.
  unused_edges <- which(!edge_used & edge_list$from == v)
  
  if (length(unused_edges) > 0) {
    # WHAT: Take the first unused edge from v to next node.
    # WHY: Mark it as used and push destination onto stack to continue DFS.
    edge_idx <- unused_edges[1]
    edge_used[edge_idx] <- TRUE
    next_node <- edge_list$to[edge_idx]
    stack <- c(stack, next_node)
    
    if (iteration <= 10) {
      cat("Iteration", iteration, ": From", v, "→", next_node, 
          "(Edge", edge_idx, ") - Stack size:", length(stack), "\n")
    }
  } else {
    # WHAT: If no unused edges, pop node and add to path (in reverse).
    # WHY: Node is finished; recording in path ensures correct traversal order.
    path <- c(v, path)
    stack <- stack[-length(stack)]
    
    if (iteration <= 10) {
      cat("Iteration", iteration, ": No edges from", v, "- backtrack. Path size:", length(path), "\n")
    }
  }
}

cat("\n✓ Eulerian Path Computed:\n")
cat("Path Length (nodes):", length(path), "\n")
cat("First 5 nodes:", paste(head(path, 5), collapse = " -> "), "\n")
cat("Last 5 nodes:", paste(tail(path, 5), collapse = " -> "), "\n\n")

# STEP 6: Prepare output for sequence reconstruction
# WHAT: Store the node path and metadata for the next script.
# WHY: The next script converts this node path into the assembled DNA sequence.

euler_data <- list(
  path = path,
  start_node = start_node,
  end_node = path[length(path)],
  path_type = path_type,
  path_length = length(path),
  kmers = kmers,
  k = k
)

cat("========== EULERIAN PATH COMPUTED ==========\n")
cat("Ready for sequence reconstruction.\n")

# 
# END OF SCRIPT
#