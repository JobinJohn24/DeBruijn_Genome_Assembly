# 
# 01_kmers.R: K-mer Extraction from DNA Sequence
# 
# PURPOSE: Extract all k-mers of length k from a DNA sequence.
# WHAT: Generate overlapping substrings of fixed length k.
# WHY: K-mers are the building blocks of De Bruijn graphs; each k-mer becomes
#      an edge, and their overlaps define the graph structure.
# INPUT: A DNA sequence (string) and k (integer).
# OUTPUT: A list of k-mers and their frequency counts.

library(stringr)

# STEP 1: Load or define the toy sequence
# WHAT: Define the ground-truth DNA sequence we want to reconstruct.
# WHY: We'll extract k-mers from this to simulate a sequencing experiment.
#      In real applications, k-mers come from a sequencer; here we generate
#      them to ensure we have perfect coverage for proof-of-concept.

sequence <- "ACGTAC"  # Toy sequence
k <- 3                # K-mer length (adjust to 5, 7, etc. for longer sequences)

cat("Original Sequence:", sequence, "\n")
cat("K-mer Length (k):", k, "\n")

# STEP 2: Extract all overlapping k-mers
# 
# WHAT: Slide a window of size k across the sequence, one position at a time.
#       Collect each substring.
# WHY: This generates all possible k-mers. In graph terms, each k-mer will
#      become an edge from its (k-1)-mer prefix to its (k-1)-mer suffix.

kmers <- character(nchar(sequence) - k + 1)

for (i in 1:(nchar(sequence) - k + 1)) {
  # WHAT: Extract substring of length k starting at position i.
  # WHY: str_sub uses 1-based indexing; we get the i-th to (i+k-1)-th character.
  kmer <- str_sub(sequence, i, i + k - 1)
  kmers[i] <- kmer
}

cat("\nExtracted K-mers:\n")
print(kmers)

# STEP 3: Count k-mer frequencies
# WHAT: Create a table of k-mer counts (how many times each appears).
# WHY: In real data, low-frequency k-mers are often errors; we may filter
#      them out. Here all should appear once (perfect data), but this prepares
#      us for error-handling in later extensions.

kmer_counts <- table(kmers)
kmer_df <- data.frame(
  kmer = names(kmer_counts),
  count = as.numeric(kmer_counts),
  row.names = NULL
)

cat("\nK-mer Frequency Table:\n")
print(kmer_df)

# 
# STEP 4: Extract (k-1)-mers (prefixes and suffixes)

# WHAT: For each k-mer, extract the (k-1)-mer prefix and suffix.

# WHY: These become the nodes in the De Bruijn graph. The prefix and suffix
#      overlap by k-2 characters, encoding the sequence information.

prefixes <- str_sub(kmers, 1, k - 1)         # First k-1 characters
suffixes <- str_sub(kmers, 2, k)             # Last k-1 characters

cat("\nPrefixes (k-1 = ", k - 1, "-mers):\n", sep = "")
print(prefixes)

cat("\nSuffixes (k-1 = ", k - 1, "-mers):\n", sep = "")
print(suffixes)

# STEP 5: Prepare output for next script
# WHAT: Store extracted k-mers and their structure for the next script.
# WHY: These will be used to build the De Bruijn graph in 02_debruijn_graph.R.
#      We pass them as global variables accessible to the next script.

# Save results in global environment for next script
kmer_data <- list(
  sequence = sequence,
  k = k,
  kmers = kmers,
  kmer_counts = kmer_counts,
  kmer_df = kmer_df,
  prefixes = prefixes,
  suffixes = suffixes
)

cat("\n========== K-MER EXTRACTION COMPLETE ==========\n")
cat("Total k-mers extracted:", length(kmers), "\n")
cat("Unique k-mers:", nrow(kmer_df), "\n")
cat("Ready for De Bruijn graph construction.\n")
