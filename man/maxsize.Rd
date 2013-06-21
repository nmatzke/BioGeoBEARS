\name{maxsize}
\alias{maxsize}
\title{Get the maximum rangesize for a given ancestral rangesize}
\usage{
  maxsize(areasizes_possible_01)
}
\arguments{
  \item{areasizes_possible_01}{A list of 0/1 values,
  indicating whether an range of that size (rangesize =
  1-based index = 1, 2, 3...) is possible (1) or not (0).}
}
\value{
  \code{max_number_of_areas} The maximum number of areas
}
\description{
  This function returns the maximum descendant rangesize
  for a given ancestral rangesize, given a list of 0/1
  values specifying the possibility of each descendant
  rangesizes.
}
\details{
  This is mostly a utility function used within
  \code{\link[base]{apply}} within other functions.
}
\note{
  Go BEARS!
}
\examples{
testval=1
areasizes_possible_01 = c(1,1,1,0,0)
maxsize(areasizes_possible_01)
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

