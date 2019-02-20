
Repguide <img src="man/figures/logo.png" align="right" alt="" width="120" />
============================================================================

The Repguide R package is a *beta* development version to facilitate the design of guideRNAs for CRISPR/dCas9 targeting of repetitive DNA sequences, such as those derived from transposable elements.

### Functionality

The basic workflow consists of exploring and selecting target sites, computing the guideRNA universe, and finding the optimal combination of guides that maximizes user-defined scoring:

<img src="man/figures/schematic.png" width="90%" style="display: block; margin: auto;" />

#### Features include:

-   Multi-species support
-   Flexible genome annotation
-   Black- and whitelisting
-   Multiple sequence alignment
-   Promoter identification
-   Consensus targeting
-   Genomic mapping
-   On and off-target scoring
-   Combinatorial optimization for multiple guides
-   Quality control reports

More details on the usage of Repguide is available in the package [vignette](https://tanaylab.bitbucket.io/Repguide/articles/Repguide.html).

### Installation

``` r
# Install BiocManager, devtools, and tgstat (in case you haven't already)
install.packages('BiocManager')
install.packages('devtools')
install.packages('tgstat', repos=c(getOption('repos'), 'https://tanaylab.bitbucket.io/repo'))

# Install Repguide
options(repos = c(getOption("repos"), BiocManager::repositories()))
devtools::install_bitbucket('tanaylab/repguide', ref='default')
```

**Note**: Repguide currently requires Unix environment. In particular it uses the Tanay group tgstat library that utilizes shared memory and distributed computing (as well as some specific optional CPU features).
