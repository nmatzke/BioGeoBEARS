\name{read_times_fn}
\alias{read_times_fn}
\title{Read in the stratification time breakpoints}
\usage{
  read_times_fn(inputs = NULL, timesfn = NULL)
}
\arguments{
  \item{inputs}{The inputs list}

  \item{timesfn}{The times filename.}
}
\value{
  \code{timeperiods} A list object
}
\description{
  The timeperiods file is just a list of times, 1 per line,
  from youngest to oldest.
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

