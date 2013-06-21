\name{meval}
\alias{meval}
\title{\code{eval()} function for use in \code{sapply}}
\usage{
  meval(equation_txt)
}
\arguments{
  \item{equation_txt}{The text of the equation to run
  \code{\link{eval}} on -- e.g., from a cell of a
  text-based transition matrix.}
}
\value{
  \code{outval} The numeric result of \code{\link{eval}}.
}
\description{
  \code{meval} is a wrapper for \code{\link{eval}}, to
  allow use in \code{sapply}.
}
\details{
  This is an attempt to speed up the use of
  \code{\link{eval}}; in general use of \code{\link{eval}}
  to convert a text version of a transition matrix to a
  numeric version with probabilities is a poor, slow
  choice; but it can be useful for examples and display
  purposes.

  See \code{\link[cladoRcpp]{cladoRcpp}} for fast C++
  implementations of transition matrix setup.
}
\note{
  Go BEARS!
}
\examples{
testval=1
d = 0.1
equation_txt = "1*d+1*d"
meval(equation_txt)
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

