\name{extract_numbers}
\alias{extract_numbers}
\title{Extract just the numbers from a string, including decimal points}
\usage{
  extract_numbers(tmpstr)
}
\arguments{
  \item{tmpstr}{An input string.}
}
\value{
  \code{x2} The list of numbers
}
\description{
  This function extracts numbers from a string.  Contiguous
  digits, including decimal points, are made into a single
  number. A list of numbers is returned.
}
\details{
  This saves you having to remember the
  \code{regexp}/\code{\link[base]{gregexpr}} code for this
  sort of thing, and makes it much easier to parse numbers
  out of the text output of various programs.
}
\note{
  Go BEARS!
}
\examples{
tmpstr = "190Ma - 65Ma"
extract_numbers(tmpstr)

tmpstr = "190.1Ma - 65.5Ma"
extract_numbers(tmpstr)
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link[base]{gregexpr}}
}

