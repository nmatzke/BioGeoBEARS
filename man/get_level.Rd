\name{get_level}
\alias{get_level}
\title{Get a node's level in the tree}
\usage{
  get_level(nodenum, t, tmplevel = 0)
}
\arguments{
  \item{nodenum}{The node number to get the parent of}

  \item{t}{An ape phylo object}

  \item{tmplevel}{A starting level (the function is
  recursive)}
}
\value{
  \code{tmplevel} The level of the node.
}
\description{
  Finds how many nodes deep a node is.
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

