\name{sourceall}
\alias{sourceall}
\title{Source all .R files in a directory, except "compile" and "package" files}
\usage{
  sourceall(path = path, pattern = "\\\\.R", ...)
}
\arguments{
  \item{path}{The path to source}

  \item{pattern}{Default is .R}

  \item{...}{Additional arguments to source}
}
\value{
  \code{path} The path that was sourced.
}
\description{
  Utility function.
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
}
\seealso{
  \code{\link[base]{source}}
}

