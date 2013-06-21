\name{get_daughters}
\alias{get_daughters}
\title{Get all the direct daughters nodes of a node}
\usage{
  get_daughters(nodenum, t)
}
\arguments{
  \item{nodenum}{The node number to get the daughters of}

  \item{t}{An ape phylo object}
}
\value{
  \code{daughter_nodenums} List of the daughter node
  numbers
}
\description{
  Get all the direct daughters nodes of a node
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
  \code{\link{findall}}, \code{\link{chainsaw2}}
}

