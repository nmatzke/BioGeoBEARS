\name{check_if_state_is_allowed}
\alias{check_if_state_is_allowed}
\title{Check if a geographic range/state is allowed, given an areas-allowed matrix.}
\usage{
  check_if_state_is_allowed(state_0based_indexes,
    areas_allowed_mat)
}
\arguments{
  \item{state_0based_indexes}{The input state is a 0-based
  vector of area indices.}

  \item{areas_allowed_mat}{A matrix (number of areas x
  number of areas) with 1s indicating allowed connections
  between areas, and 0s indicating disallowed connections.}
}
\value{
  \code{TRUE} or \code{FALSE}
}
\description{
  If the user has specified a matrix stating which areas
  are allowed to be connected (and thus have a species with
  a range in both areas), this function checks if the input
  list of areas (as a 0-based vector of areas) in a single
  state/geographic range is consistent with the
  areas-allowed matrix.
}
\details{
  This function may be used by e.g.
  \code{\link[base]{apply}}.
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

  Matzke_2012_IBS
}
\seealso{
  \code{\link[base]{apply}}
}

