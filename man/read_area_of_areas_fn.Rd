\name{read_area_of_areas_fn}
\alias{read_area_of_areas_fn}
\title{Read in the area areas by time}
\usage{
  read_area_of_areas_fn(inputs = NULL,
    area_of_areas_fn = NULL)
}
\arguments{
  \item{inputs}{The inputs list}

  \item{area_of_areas_fn}{The area-of-areas filename.}
}
\value{
  \code{list_of_area_areas_mats} A list object
}
\description{
  area_areas file is just a list of distance matrices,
  separated by blank lines, from youngest to oldest.
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
  \code{\link[stats]{convolve}}
}

