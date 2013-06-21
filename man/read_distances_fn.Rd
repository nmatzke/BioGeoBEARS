\name{read_distances_fn}
\alias{read_distances_fn}
\title{Read in the distances by time}
\usage{
  read_distances_fn(inputs = NULL, distsfn = NULL)
}
\arguments{
  \item{inputs}{The inputs list}

  \item{distsfn}{The distances filename.}
}
\value{
  \code{list_of_distances_mats} A list object
}
\description{
  Distances file is just a list of distance matrices,
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

