#' BioGeoBEARS: BioGeography with Bayesian (and Likelihood) Evolutionary Analysis in R Scripts
#'
#' \tabular{ll}{
#' Package: \tab BioGeoBEARS\cr
#' Type: \tab Package\cr
#' Version: \tab 0.2.1\cr
#' Date: \tab 2012-07-27\cr
#' License: \tab GPL (>= 3)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' Summary: This package performs model-based statistical inference for historical biogeography. This
#' includes inference of model parameters, ancestral states, and model comparison. This package 
#' performs ML (maximum-likelihood) based inference, but the same functions can easily be 
#' integrated into a Bayesian analysis via use of MCMC sampling functions from other packages.
#' 
#' Details: BioGeoBEARS allows probabilistic inference of both historical
#' biogeography (ancestral geographic ranges on a phylogeny) as well as
#' comparison of different models of range evolution.  It reproduces
#' the model available in LAGRANGE (Ree and Smith 2008), as well as
#' making available numerous additional models. For example, LAGRANGE
#' as typically run has two free parameters, d (dispersal rate, i.e.
#' the rate of range addition along a phylogenetic branch) and e
#' (extinction rate, really the rate of local range loss along a
#' phylogenetic branch). LAGRANGE also has a fixed cladogenic model
#' which gives equal probability to a number of allowed range
#' inheritance events, e.g.: (1) vicariance, (2) a new species starts
#' in a subset of the ancestral range, (3) the ancestral range is
#' copied to both species; in all cases, at least one species must have
#' a starting range of size 1.  LAGRANGE assigns equal probability to
#' each of these events, and zero probability to other events.
#' BioGeoBEARS adds an additional cladogenic event: founder-event
#' speciation (the new species jumps to a range outside of the
#' ancestral range), and also allows the relative weighting of the
#' different sorts of events to be made into free parameters, allowing
#' optimization and standard model choice procedures to pick the best
#' model. The relative probability of different descendent range sizes
#' is also parameterized and thus can also be specified or estimated.
#' The flexibility available in BioGeoBEARS also enables the natural
#' incorporation of (1) imperfect detection of geographic ranges in the
#' tips, and (2) inclusion of fossil geographic range data, when the
#' fossils are tips on the phylogeny.  Bayesian analysis has been implemented
#' through use of the "LaplacesDemon" package, however this package is 
#' now maintained off of CRAN, so its usage is not formally included in 
#' BioGeoBEARS at the current time.  
#' 
#' CITATION INFO: This package is 
#' the result of my Ph.D. research, please cite the package if you use 
#' it! Type: citation(package="BioGeoBEARS") to get the citation information.
#' 
#' See also the citation information for the sister packages,
#' \code{citation(package="rexpokit")} and \code{citation(package="cladoRcpp")}.
#'
#' @name BioGeoBEARS-package
#' @aliases BioGeoBEARS
#' @docType package
#' @title BioGeography with Bayesian (and likelihood) Evolutionary Analysis of RangeS
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @keywords package Rcpp rexpokit cladoRcpp
#' @seealso \code{\link[rexpokit]{rexpokit}} \code{\link[cladoRcpp]{cladoRcpp}}
#' @examples
#' test=1
#' 
#' # To get citation information for BioGeoBEARS, type:
#' citation(package="BioGeoBEARS")
#' 
#' # Please also cite the accessory packages I created to make BioGeoBEARS work:
#' citation(package="cladoRcpp")
#' citation(package="rexpokit")		# Roger Sidje is a coauthor of rexpokit 
#'                                  # and author of the FORTRAN EXPOKIT
#' 
 
