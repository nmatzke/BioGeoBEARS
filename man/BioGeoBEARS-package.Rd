\docType{package}
\name{BioGeoBEARS-package}
\alias{BioGeoBEARS}
\alias{BioGeoBEARS-package}
\title{BioGeography with Bayesian (and likelihood) Evolutionary Analysis of RangeS}
\description{
  BioGeoBEARS: BioGeography with Bayesian (and likelihood)
  Evolutionary Analysis in R Scripts
}
\details{
  \tabular{ll}{ Package: \tab BioGeoBEARS\cr Type: \tab
  Package\cr Version: \tab 0.1\cr Date: \tab 2012-05-30\cr
  License: \tab GPL (>= 3)\cr LazyLoad: \tab yes\cr }

  This package performs model-based statistical inference
  for historical biogeography. This includes inference of
  model parameters, ancestral states, and model comparison.
  Currently, the package performs ML (maximum-likelihood)
  based inference, but will soon be expanded to Bayesian
  analysis.

  \bold{Details:} BioGeoBEARS allows probabilistic
  inference of both historical biogeography (ancestral
  geographic ranges on a phylogeny) as well as comparison
  of different models of range evolution.  It reproduces
  the model available in LAGRANGE (Ree and Smith 2008), as
  well as making available numerous additional models. For
  example, LAGRANGE as typically run has two free
  parameters, d (dispersal rate, i.e. the rate of range
  addition along a phylogenetic branch) and e (extinction
  rate, really the rate of local range loss along a
  phylogenetic branch). LAGRANGE also has a fixed
  cladogenic model which gives equal probability to a
  number of allowed range inheritance events, e.g.: (1)
  vicariance, (2) a new species starts in a subset of the
  ancestral range, (3) the ancestral range is copied to
  both species; in all cases, at least one species must
  have a starting range of size 1.  LAGRANGE assigns equal
  probability to each of these events, and zero probability
  to other events. BioGeoBEARS adds an additional
  cladogenic event: founder-event speciation (the new
  species jumps to a range outside of the ancestral range),
  and also allows the relative weighting of the different
  sorts of events to be made into free parameters, allowing
  optimization and standard model choice procedures to pick
  the best model. The relative probability of different
  descendent range sizes is also parameterized and thus can
  also be specified or estimated. The flexibility available
  in BioGeoBEARS also enables the natural incorporation of
  (1) imperfect detection of geographic ranges in the tips,
  and (2) inclusion of fossil geographic range data, when
  the fossils are tips on the phylogeny.
}
\note{
  Go BEARS!
}
\examples{
test=1
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
  \url{https://code.google.com/p/lagrange/}

  Matzke N (2012). "Founder-event speciation in BioGeoBEARS
  package dramatically improves likelihoods and alters
  parameter inference in Dispersal-Extinction-Cladogenesis
  (DEC) analyses." _Frontiers of Biogeography_, *4*(suppl.
  1), pp. 210. ISSN 1948-6596, Poster abstract published in
  the Conference Program and Abstracts of the International
  Biogeography Society 6th Biannual Meeting, Miami,
  Florida. Poster Session P10: Historical and
  Paleo-Biogeography. Poster 129B. January 11, 2013, <URL:
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}>.

  Ree RH and Smith SA (2008). "Maximum likelihood inference
  of geographic range evolution by dispersal, local
  extinction, and cladogenesis." _Systematic Biology_,
  *57*(1), pp. 4-14. <URL:
  http://dx.doi.org/10.1080/10635150701883881>, <URL:
  http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&dopt=Citation&list_uids=18253896>.
}
\seealso{
  \code{\link{rexpokit}} \code{\link{cladoRcpp}}
}
\keyword{cladoRcpp}
\keyword{package}
\keyword{Rcpp}
\keyword{rexpokit}

