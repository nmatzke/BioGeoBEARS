\name{strsplit_whitespace}
\alias{strsplit_whitespace}
\title{Split strings on whitespace}
\usage{
  strsplit_whitespace(tmpline)
}
\arguments{
  \item{tmpline}{A string containing text.}
}
\value{
  \code{list_of_strs}
}
\description{
  This function splits strings on whitespace (spaces and
  tabs), so you don't have to remember the
  \code{regexp}/\code{grep} format codes.
}
\note{
  Go BEARS!
}
\examples{
tmpline = "Hello world see	my	tabs."
strsplit_whitespace(tmpline)
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link[base]{strsplit}}
}

