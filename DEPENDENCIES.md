# Dependencies & Environment Setup

## Required Software

- **R**: 4.0 or later
- **Operating System**: Linux, macOS, or Windows

## R Packages

### Required Packages

| Package | Version | Purpose |
|---------|---------|---------|
| `igraph` | ≥ 1.3.0 | Graph data structure and algorithms |
| `stringr` | ≥ 1.4.0 | String manipulation (k-mer extraction) |
| `stringdist` | ≥ 0.9.0 | Edit distance calculations (metrics) |

### Installation

Install all required packages at once:

```r
install.packages(c("igraph", "stringr", "stringdist"))
```

Or install individually:

```r
install.packages("igraph")
install.packages("stringr")
install.packages("stringdist")
```

### Verify Installation

```r
# Load packages to verify
library(igraph)
library(stringr)
library(stringdist)

# Check versions
packageVersion("igraph")
packageVersion("stringr")
packageVersion("stringdist")
```

## Environment Setup

### Option 1: RStudio (Recommended for beginners)

1. Download and install [RStudio Desktop](https://www.rstudio.com/products/rstudio/download/) (Free version)
2. Open RStudio
3. Create a new project in the cloned repository directory
4. Run in RStudio Console:
   ```r
   install.packages(c("igraph", "stringr", "stringdist"))
   source("RUN_ME_FIRST.R")
   ```

### Option 2: Command Line R

```bash
# Navigate to project directory
cd debruijn-assembly-r

# Open R
R

# Inside R console:
install.packages(c("igraph", "stringr", "stringdist"))
source("RUN_ME_FIRST.R")
quit()
```

### Option 3: Rscript (Non-interactive)

```bash
# Navigate to project directory
cd debruijn-assembly-r

# Run all scripts non-interactively
Rscript RUN_ME_FIRST.R
```

## System Requirements

### Minimum
- **CPU**: 1 core (single-threaded)
- **RAM**: 512 MB
- **Storage**: 50 MB for R + packages

### Recommended
- **CPU**: 2+ cores
- **RAM**: 2-4 GB
- **Storage**: 500 MB for R + packages

## Platform-Specific Notes

### Windows

1. Install R from [CRAN](https://cloud.r-project.org/bin/windows/base/)
2. Install RStudio Desktop
3. Open RStudio and run:
   ```r
   install.packages(c("igraph", "stringr", "stringdist"))
   ```

**Note**: Some users report graphviz issues on Windows; igraph falls back gracefully, so visualization works without Graphviz.

### macOS

1. Install R via Homebrew:
   ```bash
   brew install r
   ```
   Or download from [CRAN](https://cloud.r-project.org/bin/macosx/)

2. Install RStudio Desktop from [rstudio.com](https://www.rstudio.com/products/rstudio/download/)

3. In R Console:
   ```r
   install.packages(c("igraph", "stringr", "stringdist"))
   ```

### Linux (Ubuntu/Debian)

```bash
# Install R
sudo apt-get update
sudo apt-get install r-base r-base-dev

# Install RStudio (optional, recommended)
wget https://download1.rstudio.org/desktop/bionic/amd64/rstudio-1.4.1106-amd64.deb
sudo dpkg -i rstudio-1.4.1106-amd64.deb

# In R Console:
install.packages(c("igraph", "stringr", "stringdist"))
```

## Reproducible Setup with renv (Advanced)

For strict reproducibility across environments, use `renv`:

```r
# Install renv
install.packages("renv")

# Initialize renv in project
renv::init()

# Install project dependencies
renv::install(c("igraph", "stringr", "stringdist"))

# Create lock file (commit to git)
renv::snapshot()

# Others can restore environment
renv::restore()
```

