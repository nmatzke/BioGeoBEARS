\name{cornerpies}
\alias{cornerpies}
\title{Make pie charts for plotting ranges on corners}
\usage{
  cornerpies(pievals, coords, piecol, adj = c(0.5, 0.5),
    ...)
}
\arguments{
  \item{pievals}{The matrix (numnodes x numstates) of
  probabilities to plot.}

  \item{coords}{The coordinates at which to plot the
  labels.}

  \item{piecol}{The color for each possible state.}

  \item{adj}{Position adjustment; default
  \code{adj=c(0.5,0.5)}}

  \item{...}{Additional arguments to standard functions}
}
\value{
  nothing
}
\description{
  This function makes pie charts for plotting ranges on
  corners.  It makes use of \code{ape:::floating.pie.asp}
  to plot the pie charts on the corners.
}
\details{
  To get the corner coordinates, use
  \code{\link{corner_coords}}.  Please note the special
  input required in that function to get it to access a
  corner-coordinates function in the extensions data
  (\code{extdata}) directory.
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
  \code{\link{cornerlabels}}, \code{\link{corner_coords}},
  \code{\link{get_lagrange_nodenums}},
  \code{\link{LGpy_splits_fn_to_table}},
  \code{\link{LGcpp_splits_fn_to_table}}
}

