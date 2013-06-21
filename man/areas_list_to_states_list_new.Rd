\name{areas_list_to_states_list_new}
\alias{areas_list_to_states_list_new}
\title{Convert a list of areas to a list of geographic ranges (states); R version}
\usage{
  areas_list_to_states_list_new(areas = c("A", "B", "C"),
    maxareas = length(areas), include_null_range = TRUE,
    split_ABC = TRUE)
}
\arguments{
  \item{areas}{a list of areas (character or number; the
  function converts these to numbers, starting with 0)}

  \item{maxareas}{maximum number of areas in this analyses}

  \item{include_null_range}{\code{TRUE} or \code{FALSE},
  should the \code{NULL} range be included in the possible
  states? (e.g., LAGRANGE default is yes)}

  \item{split_ABC}{\code{TRUE} or \code{FALSE} If
  \code{TRUE} the output will consist of a list of lists
  (c("A","B","C"), c("A","B"), c("A","D"), etc.); if
  \code{FALSE}, the list of areas will be collapsed ("ABC",
  "AB", "AD", etc.).}
}
\value{
  \code{states_list} A list of the states.
}
\description{
  R version of areas_list_to_states_list_old, which makes
  use of \code{\link[cladoRcpp]{cladoRcpp}}'s
  \code{\link[cladoRcpp]{rcpp_areas_list_to_states_list}}.
}
\details{
  This is the original R version of the function which
  converts a list of possible areas to a list of all
  possible states (geographic ranges).  This gets slow for
  large numbers of areas.

  The function is mostly replaced by
  \code{\link[cladoRcpp]{rcpp_areas_list_to_states_list}}
  in optimized code, but is still used in some places for
  display purposes.
}
\note{
  Go BEARS!
}
\examples{
areas = c("A","B","C")
areas_list_to_states_list_new(areas=areas, maxareas=length(areas),
include_null_range=TRUE, split_ABC=TRUE)
areas_list_to_states_list_new(areas=areas, maxareas=length(areas),
include_null_range=TRUE, split_ABC=FALSE)
areas_list_to_states_list_new(areas=areas, maxareas=length(areas),
include_null_range=FALSE, split_ABC=TRUE)
areas_list_to_states_list_new(areas=areas, maxareas=length(areas),
include_null_range=FALSE, split_ABC=FALSE)
areas_list_to_states_list_new(areas=areas, maxareas=2,
include_null_range=TRUE, split_ABC=TRUE)
areas_list_to_states_list_new(areas=areas, maxareas=2,
include_null_range=TRUE, split_ABC=FALSE)
areas_list_to_states_list_new(areas=areas, maxareas=2,
include_null_range=FALSE, split_ABC=TRUE)
areas_list_to_states_list_new(areas=areas, maxareas=2,
include_null_range=FALSE, split_ABC=FALSE)
areas_list_to_states_list_new(areas=areas, maxareas=1,
include_null_range=TRUE, split_ABC=TRUE)
areas_list_to_states_list_new(areas=areas, maxareas=1,
include_null_range=TRUE, split_ABC=FALSE)
areas_list_to_states_list_new(areas=areas, maxareas=1,
include_null_range=FALSE, split_ABC=TRUE)
areas_list_to_states_list_new(areas=areas, maxareas=1,
include_null_range=FALSE, split_ABC=FALSE)
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
  \url{https://code.google.com/p/lagrange/}

  Matzke_2012_IBS

  ReeSmith2008
}
\seealso{
  \code{\link{numstates_from_numareas}},
  \code{\link{rcpp_areas_list_to_states_list}}
}

