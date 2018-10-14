# BioGeoBEARS
BioGeography with Bayesian (and likelihood) Evolutionary Analysis with R Scripts

# UPDATES, September 2018

The GitHub repository contains all updates formerly posted on PhyloWiki. This includes functions for Biogeographical Stochastic Mapping (BSM) and trait-dependent dispersal models. Further updates will be posted to GitHub.

To install the GitHub version of BioGeoBEARS, first:

1. Install the new versions of [rexpokit](https://CRAN.R-project.org/package=rexpokit) and [cladoRcpp](https://CRAN.R-project.org/package=cladoRcpp), both available on CRAN (which gives you binaries for easy installation; the GitHub versions have to be compiled from source):

```
install.packages("rexpokit")
install.packages("cladoRcpp")
```

2. Install the new version of [BioGeoBEARS from GitHub](https://github.com/nmatzke/BioGeoBEARS), using [devtools](https://CRAN.R-project.org/package=devtools):

```
library(devtools)
devtools::install_github(repo="nmatzke/BioGeoBEARS")
```

**Benefit:** you will no longer have to use the tedious source() commands that were in the BioGeoBEARS example script.


# DESCRIPTION

BioGeoBEARS allows probabilistic inference of both historical biogeography (ancestral geographic ranges on a phylogeny) as well as comparison of different models of range evolution. It reproduces the model available in LAGRANGE (Ree and Smith 2008), as well as making available numerous additional models. For example, LAGRANGE as typically run has two free parameters, d (dispersal rate, i.e. the rate of range addition along a phylogenetic branch) and e (extinction rate, really the rate of local range loss along a phylogenetic branch). LAGRANGE also has a fixed cladogenic model which gives equal probability to a number of allowed range inheritance events, e.g.: (1) vicariance, (2) a new species starts in a subset of the ancestral range, (3) the ancestral range is copied to both species; in all cases, at least one species must have a starting range of size 1. LAGRANGE assigns equal probability to each of these events, and zero probability to other events. BioGeoBEARS adds an additional cladogenic event: founder-event speciation (the new species jumps to a range outside of the ancestral range), and also allows the relative weighting of the different sorts of events to be made into free parameters, allowing optimization and standard model choice procedures to pick the best model. The relative probability of different descendent range sizes is also parameterized and thus can also be specified or estimated. The flexibility available in BioGeoBEARS also enables the natural incorporation of (1) imperfect detection of geographic ranges in the tips, and (2) inclusion of fossil geographic range data, when the fossils are tips on the phylogeny. Bayesian analysis has been implemented through use of the "LaplacesDemon" package, however this package is now maintained off of CRAN, so its usage is not formally included in BioGeoBEARS at the current time. CITATION INFO: This package is the result of my Ph.D. research, please cite the package if you use it! Type: citation(package="BioGeoBEARS") to get the citation information.


# UNIT TESTS

BioGeoBEARS version 1.1 also includes 156 unit-tests in the "tests" directory, using the R package "testthat". These check the likelihood calculations, ancestral state probabilities, and ML optimizations for regular, time-stratified, and trait-dependent models.

The tests all run successfully on my Mac, but they take too long to build in Travis-CI (which has a time-limit of 50 minutes). See [Build #24 for a successful build](https://travis-ci.org/nmatzke/BioGeoBEARS/builds/439942601) that does only some of the tests.

(The present version of BioGeoBEARS will install from GitHub regardless of the unit-tests issue.)
