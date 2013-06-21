\name{read_areas_allowed_fn}
\alias{read_areas_allowed_fn}
\title{Read in the area areas by time}
\usage{
  read_areas_allowed_fn(inputs = NULL,
    areas_allowed_fn = NULL)
}
\arguments{
  \item{inputs}{The inputs list}

  \item{areas_allowed_fn}{The areas-allowed filename.}
}
\value{
  \code{list_of_areas_allowed_mats} A list object
}
\description{
  areas_allowed file is just a list of distance matrices,
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

