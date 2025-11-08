# 👨‍🔬🧬 DeBruijn Genome Assembly 👨‍🔬🧬

A complete implementation of De Bruijn graph-based genome assembly in R, demonstrating how real sequencing tools (SPAdes, Velvet, ABySS) reconstruct DNA from overlapping fragments.

## 📋Overview

This project implements the complete De Bruijn graph assembly pipeline:

1. **Extract k-mers** from DNA sequence
2. **Build a directed graph** where edges represent k-mers
3. **Find an Eulerian path** through all edges (Hierholzer's algorithm)
4. **Reconstruct the sequence** from the path
5. **Calculate quality metrics** (identity, contiguity, complexity)

## 📖 How It Works

### The Algorithm: 4 Steps

#### Step 1: K-mer Extraction

Break DNA into overlapping k-mers (words of length k):
```
Input DNA:    ACGTAC
K-value:      3

Extract sliding windows:
  Position 1: ACG
  Position 2:  CGT
  Position 3:   GTA
  Position 4:    TAC

Output k-mers: [ACG, CGT, GTA, TAC]
```

#### Step 2: Build De Bruijn Graph

Create a directed graph where:
- **Nodes** = (k-1)-mers (prefixes and suffixes of k-mers)
- **Edges** = k-mers (connections between nodes)
```
For k=3, the k-mers connect (k-1)=2-mers:

K-mers → Connections
ACG    → AC → CG
CGT    → CG → GT
GTA    → GT → TA
TAC    → TA → AC

Graph visualization:
    AC
   /  \
  CG  TA
   \  /
    GT
```

#### Step 3: Find Eulerian Path

Use Hierholzer's algorithm to find a path visiting every edge exactly once:
```
Start: AC (arbitrary)
Path: AC → CG → GT → TA → AC
Edges: ACG → CGT → GTA → TAC
```

#### Step 4: Reconstruct Sequence

Concatenate the path, using overlaps:
```
Path:     AC → CG → GT → TA → AC
Extract:  AC + (G)  + (T) + (A) + (C)
Result:   ACGTAC ✓
```

### 🧮 Algorithm 

| Operation | Time | Space |
|-----------|------|-------|
| K-mer extraction | O(n × k) | O(n) |
| Graph construction | O(n) | O(V + E) |
| Eulerian path | O(V + E) | O(V + E) |
| **Total** | **O(n × k)** | **O(n)** |

Where: n = sequence length, k = k-mer size, V = nodes, E = edges

## 📚Results

### Toy Data Assembly

**Input:** ACGTAC (6 bp), K=3

| Metric | Result |
|--------|--------|
| Reconstructed sequence | ACGTAC |
| Match | ✓ TRUE |
| Identity | 100% |
| Mismatches | 0 |
| Nodes | 4 |
| Edges | 4 |
| Eulerian circuit | ✓ Found |

**Graph Properties:**
- All 4 nodes balanced (in-degree = out-degree)
- Single connected component
- Perfect cycle: AC → CG → GT → TA → AC

### Real mtDNA Analysis

**Input:** Human mitochondrial DNA (16,565 bp), tested with k = 5, 7, 10, 15, 20, 25

| K | # K-mers | # Nodes | # Edges | Balanced Nodes | Connected | Status |
|---|----------|---------|---------|----------------|-----------|--------|
| 5 | 16,565 | 260 | 1,022 | 248/260 (95%) | ✓ YES | Ambiguous |
| 7 | 16,563 | 3,499 | 8,070 | 2,113/3,499 (60%) | ✓ YES | Medium |
| 10 | 16,560 | 15,136 | 16,129 | 14,585/15,136 (96%) | ✓ YES | Good |
| **15** | **16,555** | **16,554** | **16,550** | **16,550/16,554 (99.98%)** | **✓ YES** | **✅ OPTIMAL** |
| 20 | 16,550 | 16,551 | 16,550 | 16,549/16,550 (99.99%) | ✓ YES | Good |
| 25 | 16,545 | 16,546 | 16,545 | 16,544/16,546 (99.99%) | ✓ YES | Good |

**Key Finding:** K=15 is optimal for mtDNA
- Achieves 99.98% balanced nodes
- 99.99% k-mer uniqueness (high specificity)
- Maintains graph connectivity
- All Eulerian paths found successfully

## 📊 Visualizations

### De Bruijn Graph Structures

#### Circular Layout (Toy Data - K=3)

Shows the Eulerian circuit as a perfect cycle:
```
The 4 nodes (AC, CG, GT, TA) form a directed cycle.
Each edge represents a k-mer:
- ACG: AC → CG
- CGT: CG → GT
- GTA: GT → TA
- TAC: TA → AC

Following the cycle reconstructs: ACGTAC
```
## Quick Start

### Installation
```bash
cd toy_genome
R
```

In R:
```r
# Install required packages
install.packages(c("igraph", "stringr", "stringdist"))

# Run the complete pipeline
source("RUN_ME_FIRST.R")
```
## Project Structure
```
toy_genome/
├── RUN_ME_FIRST.R              # Main entry point - run this!
├── 01_kmers.R                  # K-mer extraction
├── 02_debruijn_graph.R         # Graph construction
├── 03_eulerian_path.R          # Eulerian path finding (Hierholzer)
├── 04_sequence_reconstruct.R   # Sequence assembly
├── 05_metrics.R                # Quality metrics calculation
├── toy_data.R                  # Sample sequences
├── visualize_graph.R           # Graph visualization
├── test_k_values.R             # K-value analysis (k=2-5)
├── test_real_mtdna.R           # Real data testing (mtDNA, k=5-25)
├── RESULTS.md                  # Detailed analysis & findings
├── README.md                   # This file
├── data/
│   └── mtDNA_fragment.fasta    # Real human mtDNA (16.5k bp)
└── example_output/
    ├── debruijn_graph_labeled.png          # Toy data visualization
    └── debruijn_graph_visualization.png    # mtDNA visualization
```

### Requirements

- **R version:** 3.6 or later
- **Required packages:** igraph, stringr, stringdist
- **RAM:** 1 GB minimum (laptop-friendly)
- **OS:** Mac, Linux, Windows

### Dependencies

| Package | Purpose | Installation |
|---------|---------|--------------|
| igraph | Graph construction & algorithms | `install.packages("igraph")` |
| stringr | String manipulation | `install.packages("stringr")` |
| stringdist | Edit distance calculations | `install.packages("stringdist")` |

### Hierholzer's Algorithm

Finds Eulerian paths in O(V + E) time:
```
1. Start at source node (out-degree > in-degree)
   or any node if all balanced

2. Use depth-first search with stack:
   - Follow unused edges from current node
   - Mark each edge as used
   - When stuck, backtrack and record node

3. Reverse the recorded path

4. Result: Path visiting every edge exactly once
```


