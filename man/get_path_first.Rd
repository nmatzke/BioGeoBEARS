\name{get_path_first}
\alias{get_path_first}
\title{Get the text that comes before the last slash}
\usage{
  get_path_first(inpath, addslash = "FALSE")
}
\arguments{
  \item{inpath}{A string of class
  \code{\link[base]{character}}.}

  \item{addslash}{If \code{TRUE}, add a slash at the end of
  the path.}
}
\value{
  \code{outpath} A string with the full path, without the
  file.
}
\description{
  Extracts the path from a full path, removing the
  filename.
}
\note{
  Go BEARS!
}
\examples{
get_path_first("/Users/nickm/Library/Psychotria_geog.data")
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link{get_path_last}}
}

