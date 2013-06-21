\name{add_to_downpass_labels}
\alias{add_to_downpass_labels}
\title{Iterate up and down a tree in C++ LAGRANGE downpass order}
\usage{
  add_to_downpass_labels(tr, downpass_node_matrix,
    currnode)
}
\arguments{
  \item{tr}{A \code{\link[ape]{phylo}} tree object.}

  \item{downpass_node_matrix}{A matrix (\code{tr$Nnode}
  rows, 2 columns). Column 1 has R's native internal
  numbering scheme, and column 2 has the node numbers in a
  LAGRANGE downpass.}

  \item{currnode}{The current node being viewed}
}
\value{
  \code{downpass_node_matrix} A matrix containing node
  numbers.
}
\description{
  This is the utility function for
  \code{\link{get_lagrange_nodenums}}, which traces a tree
  down and up in C++ LAGRANGE's downpass order.
}
\details{
  This returns a matrix containing (column 1) R's native
  internal numbering scheme, and (column 2) the node
  numbers in a LAGRANGE downpass.  Note that this is
  different from \code{LAGRANGE}'s downpass ordering (see
  \code{\link{get_lagrange_nodenums}}).
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
  \code{\link{get_lagrange_nodenums}}
}

