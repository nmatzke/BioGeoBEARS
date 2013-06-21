\name{get_path_last}
\alias{get_path_last}
\title{Get the text that comes after the last slash}
\usage{
  get_path_last(path)
}
\arguments{
  \item{path}{A string of class
  \code{\link[base]{character}}.}
}
\value{
  \code{lastword} A string with the filename, without the
  path.
}
\description{
  Extracts the filename from a full path.
}
\note{
  Go BEARS!
}
\examples{
get_path_last("/Users/nickm/Psychotria_geog.data")
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
  \url{https://code.google.com/p/lagrange/}

  Matzke_2012_IBS

  ReeSmith2008

  FosterIdiots
}
\seealso{
  \code{\link{get_path_first}}
}

