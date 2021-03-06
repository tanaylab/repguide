
# Repguide <img src="man/figures/logo.png" align="right" alt="" width="120" />

The Repguide R package is a *beta* development version to facilitate the
design of guideRNAs for CRISPR/dCas9 targeting of repetitive DNA
sequences, such as those derived from transposable elements.

### Functionality

The basic workflow consists of exploring and selecting target sites,
computing the guideRNA universe, and finding the optimal combination of
guides that maximizes target-specific
coverage.

<img src="man/figures/schematic.png" width="90%" style="display: block; margin: auto;" />

#### Features include:

  - Multi-species support
  - Flexible genome annotation
  - Black- and whitelisting
  - Multiple sequence alignment
  - Promoter identification
  - Consensus positional targeting
  - Genome-wide mapping
  - On and off-target scoring
  - Combinatorial optimization (multiple guides)
  - High-quality reports

More details on the usage of Repguide is available in the package
[vignette](https://tanaylab.github.io/Repguide/index.html).

### Installation

``` r
# Install BiocManager (in case you haven't already)
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

# Install Repguide
remotes::install_github(repo = 'tanaylab/repguide')
```

**Note**: Repguide requires R version 3.5 or higher. The package was
tested on linux and macbook machines, but is currently not compatible on
Windows. A typical application will require at least 4 GB RAM, but
heavier use cases (e.g. large guideRNA universe) may require 8 GB RAM or
more. For improved speed performance, we recommend a multi-CPU/core
workstation.
