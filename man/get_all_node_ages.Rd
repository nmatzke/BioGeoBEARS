\name{get_all_node_ages}
\alias{get_all_node_ages}
\title{Get the ages of all the nodes in the tree (above the root)}
\usage{
  get_all_node_ages(obj)
}
\arguments{
  \item{obj}{An ape phylo object}
}
\value{
  \code{TF_tips} The age (from the root) of each node.
}
\description{
  A utility function. Use of \code{\link[ape]{dist.nodes}}
  may be slow.
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
  \code{\link{prt}}, \code{\link{chainsaw2}}
}

