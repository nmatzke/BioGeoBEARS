\name{cornerlabels}
\alias{cornerlabels}
\title{Make labels for plotting ranges on corners}
\usage{
  cornerlabels(text, coords, bg = "green3", col = "black",
    adj = c(0.5, 0.5), ...)
}
\arguments{
  \item{text}{The text to put at the corners.}

  \item{coords}{The coordinates at which to plot the
  labels}

  \item{bg}{The background color}

  \item{col}{The text color}

  \item{adj}{Position adjustment; default
  \code{adj=c(0.5,0.5)}}

  \item{...}{Additional arguments to standard functions}
}
\value{
  nothing
}
\description{
  This function makes labels for plotting ranges on
  corners.
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

  ReeSmith2008
}
\seealso{
  \code{\link{cornerpies}}, \code{\link{corner_coords}},
  \code{\link{get_lagrange_nodenums}},
  \code{\link{LGpy_splits_fn_to_table}},
  \code{\link{LGcpp_splits_fn_to_table}}
}

