
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ADViSELipidomics <img src="man/figures/NewLogoAL.png" align="right" height="139"/>

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-stable-succes.svg)](https://www.tidyverse.org/lifecycle/#stable)
<!-- badges: end -->

ADViSELipidomics is a novel Shiny app for the preprocessing, analysis,
and visualization of lipidomics data. It copes with the outputs from
LipidSearch and LIQUID for lipid identification and quantification, and
with data available from the Metabolomics Workbench. ADViSELipidomics
extracts information by parsing lipid species (using LIPID MAPS
classification) and, together with information available on the samples,
allows performing several exploratory and statistical analyses. In the
presence of internal lipid standards in the experiment, ADViSELipidomics
can normalize the data matrix, providing absolute values of
concentration per lipid and sample. Moreover, it allows the
identification of differentially abundant lipids in simple and complex
experimental designs, dealing with batch effect correction.

## Installation

ADViSELipidomics is a stand-alone application implemented using the R
language (R \> 4.0) and the Shiny libraries. It can be installed as any
other R package on several operating systems (Windows, macOS and Linux).
Before installing the package you have to perform few supplementary
steps based on your operating systems:

-   **Windows (tested on Windows 10 64bit)** 

Before installing the package you need also to install Rtools from the
following link: <https://cran.r-project.org/bin/windows/Rtools>.

-   **MacOS**  

``` r
brew install imagemagick@6
brew install cairo
```

-   **Ubuntu (tested on su 18.04).**  

<!-- -->

    sudo apt install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev
    sudo apt-get install libcairo2-dev
    sudo apt-get install libxt-dev
    sudo apt install libmagick++-dev
    sudo apt-get install libc6
    sudo apt-get install libnlopt-dev
