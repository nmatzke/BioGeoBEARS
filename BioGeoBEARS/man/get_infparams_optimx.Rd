\name{get_infparams_optimx}
\alias{get_infparams_optimx}
\title{Get the inferred parameters from an ML optimization}
\usage{
  get_infparams_optimx(results_object, inffn)
}
\arguments{
  \item{results_object}{The results returned by
  \code{\link{bears_2param_standard_fast}} or a similar
  function.}

  \item{inffn}{The filename holding the results_object,
  which specifies which model was run.}
}
\value{
  \code{infparams} The vector of inferred parameters.
}
\description{
  This function extracts the ML parameter values, and
  associated statistics and codes, from the
  \code{results_object} returned by
  \code{\link{bears_2param_standard_fast}} and similar
  functions.
}
\details{
  The function has subroutines for recognizing a variety of
  currently-implemented models, assuming they used
  \code{\link[optimx]{optimx}} internally to do the ML
  search.  New models would require addition of new
  subroutines.

  \code{\link{get_infparams_optimx}} and
  \code{\link{get_infparams_optimx_nosim}} differ only in
  the format of the filenames.
}
\note{
  Go BEARS!
}
\examples{
testval=1
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

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
}
\seealso{
  \code{\link{get_infparams_optimx_nosim}},
  \code{\link{bears_2param_standard_fast}},
  \code{\link{get_inf_LgL_etc_optimx}}
}

