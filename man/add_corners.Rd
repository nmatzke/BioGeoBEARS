\name{add_corners}
\alias{add_corners}
\title{Iterate up through a plotted tree, getting the coordinates of the corners}
\usage{
  add_corners(startnode, tr, nodecoords, corners_list)
}
\arguments{
  \item{startnode}{The node to start at (this is a
  recursive function)}

  \item{tr}{A tree object in \code{\link[ape]{phylo}}
  format.}

  \item{nodecoords}{The accumulating list of node
  coordinates}

  \item{corners_list}{The accumulating list of corners}
}
\value{
  \code{corners_list}
}
\description{
  What it says.
}
\note{
  Go BEARS!
}
\examples{
blah=1
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link[ape]{phylo}}, \code{\link{get_nodenums}}
}

