\name{get_fn_prefix}
\alias{get_fn_prefix}
\title{Get everything BEFORE the last suffix (.nex or whatever)}
\usage{
  get_fn_prefix(fn)
}
\arguments{
  \item{fn}{The input filename.}
}
\value{
  \code{prefix} The output string.
}
\description{
  Extracts the string from before the last suffix.  I.e.,
  "filename.nex" becomes "filename".
}
\note{
  Go BEARS!
}
\examples{
get_fn_prefix("/Users/nickm/Library/R/Psychotria_geog.data")
get_fn_prefix("Psychotria_geog.data")
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link{get_path_last}}, \code{\link{get_path_first}}
}

