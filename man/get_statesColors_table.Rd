\name{get_statesColors_table}
\alias{get_statesColors_table}
\title{Make a color table for each area and their combinations}
\usage{
  get_statesColors_table(areanames = c("K", "O", "M", "H"))
}
\arguments{
  \item{areanames}{A list of the area names.}
}
\value{
  \code{statesColors_table} A table giving the colors for
  each state.
}
\description{
  Given a list of areas, make a color table for the various
  combinations.
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

  Matzke_2012_IBS

  ReeSmith2008
}
\seealso{
  \code{\link{get_lagrange_nodenums}},
  \code{\link{LGpy_splits_fn_to_table}},
  \code{\link{LGcpp_splits_fn_to_table}}
}

