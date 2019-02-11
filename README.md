---
output:
  md_document:
    variant: markdown_github
---

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

# Repguide

The Repguide R package facilitates the design of guideRNAs for CRISPR/dCas9 targeting of repetitive DNA sequences, such as those derived from transposable elements.
The basic workflow consists of chosing the target sites, computating and scoring guideRNAs, and finding the optimal combination of multiple guides to maximize and minimize on- and off-targeting, respectively.
More details on the usage of the Repguide pipeline is available in the package vignette.

#### Installation

Make sure you have BiocManager and devtools installed. Then run:
```{r, eval=FALSE}
options(repos = c(getOption("repos"), BiocManager::repositories()))
install.packages('tgstat', repos=c(getOption('repos'), 'https://tanaylab.bitbucket.io/repo'))
devtools::install_bitbucket('tanaylab/repguide', ref='default')
```

NOTE: Repguide currently requires Unix environment.