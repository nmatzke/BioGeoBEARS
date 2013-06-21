\name{moref}
\alias{moref}
\title{print to screen the header of a file}
\usage{
  moref(fn, printnotcat = FALSE)
}
\arguments{
  \item{fn}{A filename.}

  \item{printnotcat}{If \code{TRUE}, use
  \code{\link[base]{print}} instead of
  \code{\link[base]{cat}}. Default \code{FALSE}.}
}
\value{
  Nothing returned.
}
\description{
  This does the rough equivalent of the \code{UNIX}
  function \code{more}, but within R.
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
  \code{\link[base]{scan}}
}

