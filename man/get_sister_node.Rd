\name{get_sister_node}
\alias{get_sister_node}
\title{Get the node sister to two nodes}
\usage{
  get_sister_node(tr, nodepair)
}
\arguments{
  \item{tr}{A \code{\link[ape]{phylo}} tree object.}

  \item{nodepair}{A vector (length 2) with the node numbers
  of two nodes/tips.}
}
\value{
  \code{moms_sister} The aunt node.
}
\description{
  Input two sister nodes, returns their "aunt".  Assumes a
  binary tree.
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
  \code{\link{add_to_downpass_labels}}
}

