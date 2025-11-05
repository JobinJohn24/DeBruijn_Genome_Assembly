# Visualize_graph.R: Visualize the De Bruijn Graph 
# This script creates a visual representation of the De Bruijn graph
# Run this AFTER running the main pipeline (RUN_ME_FIRST.R)

library(igraph)

cat("\n=VISUALIZING DE BRUIJN GRAPH \n\n")

# Check if graph data exists
if (!exists("dbg_data")) {
  cat("Error: Graph data not found. Run RUN_ME_FIRST.R first!\n")
  stop()
}

# Get graph
g <- dbg_data$graph
kmers <- dbg_data$kmers
k <- dbg_data$k

cat("Graph Summary:\n")
cat("  Nodes:", vcount(g), "\n")
cat("  Edges:", ecount(g), "\n")
cat("  K-value:", k, "\n\n")

# Create High-Quality Visualization
# Set up the plot device
png("example_output/debruijn_graph_visualization.png", width = 1200, height = 1000)

# Calculate layout (Fruchterman-Reingold for nice spacing)
layout <- layout_with_fr(g)

# Plot the graph
plot(g,
     layout = layout,
     vertex.size = 8,
     vertex.color = "lightblue",
     vertex.frame.color = "darkblue",
     vertex.label.cex = 0.7,
     edge.arrow.size = 0.5,
     edge.arrow.width = 1.5,
     edge.color = "gray60",
     main = sprintf("De Bruijn Graph (K=%d, Nodes=%d, Edges=%d)", 
                    k, vcount(g), ecount(g)),
     sub = "Nodes = (k-1)mers, Edges = kmers")

dev.off()

cat("✓ Graph visualization saved to: example_output/debruijn_graph_visualization.png\n\n")

# 
# Create a Simpler Version for Small Graphs (Toy Data)

if (vcount(g) <= 50) {
  cat("Graph is small enough for detailed labeling...\n\n")
  
  png("example_output/debruijn_graph_labeled.png", width = 1400, height = 1000)
  
  # Use circular layout for small graphs (clearer)
  layout_circular <- layout_in_circle(g)
  
  plot(g,
       layout = layout_circular,
       vertex.size = 20,
       vertex.color = "lightgreen",
       vertex.frame.color = "darkgreen",
       vertex.label.cex = 1.2,
       vertex.label.dist = 2,
       edge.arrow.size = 0.7,
       edge.arrow.width = 2,
       edge.color = "darkgray",
       edge.width = 2,
       main = sprintf("De Bruijn Graph - Circular Layout (K=%d)", k),
       sub = "Eulerian path finds a cycle through all edges")
  
  dev.off()
  
  cat("✓ Labeled circular graph saved to: example_output/debruijn_graph_labeled.png\n\n")
}

# Create ASCII Diagram (for markdown documentation)

cat(" ASCII GRAPH REPRESENTATION \n\n")

if (vcount(g) <= 6) {
  cat("Small graph detected! Here's a text representation:\n\n")
  
  # Get edges
  edges <- as_data_frame(g, what = "edges")
  
  cat("Nodes:\n")
  for (i in 1:vcount(g)) {
    cat(sprintf("  %s\n", V(g)$name[i]))
  }
  
  cat("\nEdges (k-mers):\n")
  for (i in 1:nrow(edges)) {
    cat(sprintf("  %s → %s\n", edges$from[i], edges$to[i]))
  }
  
  cat("\nEulerian Path (if computed):\n")
  if (exists("euler_data")) {
    path_str <- paste(euler_data$path, collapse = " → ")
    cat(sprintf("  %s\n", path_str))
  }
}

cat("\n VISUALIZATION COMPLETE \n\n")