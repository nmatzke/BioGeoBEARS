\name{get_inf_LgL_etc_optimx}
\alias{get_inf_LgL_etc_optimx}
\title{Get the inferred parameters from a results object (utility function)}
\usage{
  get_inf_LgL_etc_optimx(results_object)
}
\arguments{
  \item{results_object}{The results returned by
  \code{\link{bears_2param_standard_fast}} or a similar
  function.}
}
\value{
  \code{infparams} The vector of inferred parameters.
}
\description{
  This function extracts the ML parameter values from the
  \code{results_object} returned by
  \code{\link{bears_2param_standard_fast}} and similar
  functions.
}
\details{
  This is primarily a utility function for
  \code{\link{get_infparams_optimx}}.
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

  Matzke_2012_IBS
}
\seealso{
  \code{\link{bears_2param_standard_fast}},
  \code{\link{get_infparams_optimx}}
}

