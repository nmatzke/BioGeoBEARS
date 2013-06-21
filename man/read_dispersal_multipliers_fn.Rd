\name{read_dispersal_multipliers_fn}
\alias{read_dispersal_multipliers_fn}
\title{Read in the hard-coded dispersal multipliers from file}
\usage{
  read_dispersal_multipliers_fn(inputs = NULL,
    dispersal_multipliers_fn = NULL)
}
\arguments{
  \item{inputs}{The inputs list}

  \item{dispersal_multipliers_fn}{The dispersal multipliers
  filename.}
}
\value{
  \code{list_of_dispersal_multipliers_mats} A list object
}
\description{
  dispersal_multipliers file is just a list of distance
  matrices, separated by blank lines, from youngest to
  oldest
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

