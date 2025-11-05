# 
# 04_reconstruct_sequence.R: Sequence Reconstruction from Eulerian Path
# 
# PURPOSE: Convert the Eulerian node path back into the original DNA sequence.
# WHAT: Concatenate (k-1)-mer nodes and the last character of each edge.
# WHY: Each node is a (k-1)-mer; moving along edges adds one character per step.
#      Walking the Eulerian path reconstructs the full sequence.
# 

library(stringr)

# STEP 1: Verify Eulerian path is loaded
# WHAT: Check that 03_eulerian_path.R was sourced successfully.
# WHY: We need the path, original k-mers, and k value.

if (!exists("euler_data")) {
  stop("Eulerian path not found. Run 03_eulerian_path.R first.")
}

path <- euler_data$path
kmers <- euler_data$kmers
k <- euler_data$k
original_sequence <- kmer_data$sequence

cat("Reconstructing Sequence...\n")
cat("Path length (nodes):", length(path), "\n")
cat("K-mer length:", k, "\n")

# STEP 2: Initialize reconstructed sequence
# WHAT: Start with the first node (a (k-1)-mer).
# WHY: The first node is the (k-1)-mer prefix; subsequent nodes add one char.

reconstructed <- path[1]  # First node is a (k-1)-mer

cat("Starting node (first ", k - 1, "-mer):", reconstructed, "\n")

# STEP 3: Append one character per subsequent node.
# WHAT: For each subsequent node in the path, append its last character.
# WHY: The overlap between consecutive nodes is k-2 characters; the new
#      character is at position k (last). This builds the full sequence.

if (length(path) > 1) {
  for (i in 2:length(path)) {
    # WHAT: Extract the current node.
    # WHY: It's a (k-1)-mer; we only need its last character.
    current_node <- path[i]
    
    # WHAT: Get the last character of the current node.
    # WHY: This represents the new base added when moving from node i-1 to node i.
    new_char <- str_sub(current_node, -1, -1)
    
    # WHAT: Append this character to the reconstructed sequence.
    # WHY: This incrementally builds the full sequence as we walk the path.
    reconstructed <- paste0(reconstructed, new_char)
  }
}

cat("\n========== SEQUENCE RECONSTRUCTION ==========\n")
cat("Original Sequence    :", original_sequence, "\n")
cat("Reconstructed Seq    :", reconstructed, "\n")
cat("Lengths Match        :", nchar(original_sequence) == nchar(reconstructed), "\n")

# STEP 4: Verify reconstruction accuracy
# WHAT: Compare reconstructed to original sequence character-by-character.
# WHY: Validates that the De Bruijn graph and Eulerian path correctly encode
#      the original sequence. In real assembly, this check may fail if data
#      is noisy or incomplete.

sequences_match <- (original_sequence == reconstructed)

if (sequences_match) {
  cat("\n✓ PERFECT RECONSTRUCTION: Sequences match exactly!\n")
} else {
  cat("\n✗ MISMATCH DETECTED:\n")
  cat("  Original:     ", original_sequence, "\n")
  cat("  Reconstructed:", reconstructed, "\n")
}

# STEP 5: Prepare output for metrics calculation
# WHAT: Store the reconstructed sequence and comparison results.
# WHY: The next script uses this to compute assembly quality metrics.

recon_data <- list(
  original_sequence = original_sequence,
  reconstructed_sequence = reconstructed,
  sequences_match = sequences_match,
  path = path,
  start_node = euler_data$start_node,
  end_node = euler_data$end_node,
  path_type = euler_data$path_type,
  k = k
)

cat("\n= RECONSTRUCTION COMPLETE =\n")
cat("Ready for metrics calculation.\n")
