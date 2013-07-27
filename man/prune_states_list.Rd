\name{prune_states_list}
\alias{prune_states_list}
\title{Cut down the states list according to areas_allowed_mat}
\usage{
  prune_states_list(states_list_0based_index,
    areas_allowed_mat)
}
\arguments{
  \item{states_list_0based_index}{A states_list, 0-based,
  e.g. from
  \code{\link[cladoRcpp]{rcpp_areas_list_to_states_list}}}

  \item{areas_allowed_mat}{The matrix of area combinations
  allowed (represented by 1s)}
}
\value{
  \code{states_list_0based_index_new} A 0-based list of
  allowed states/ranges
}
\description{
  Go through a list of states.  Remove states that
  represent areas disallowed according to
  areas_allowed_mat.  It is assumed (crucial!) that the
  areas in the \code{states_list}, and in the
  \code{areas_allowed_mat}, have the same order.
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
  \code{\link[cladoRcpp]{rcpp_areas_list_to_states_list}}
}

