\name{get_indices_of_tip_nodes}
\alias{get_indices_of_tip_nodes}
\title{Get TRUE/FALSE for nodes being tips}
\usage{
  get_indices_of_tip_nodes(obj)
}
\arguments{
  \item{obj}{An ape phylo object}
}
\value{
  \code{tip_indices} The node numbers of the tips.
}
\description{
  A utility function that returns indices (node numbers) of
  the tips. This mostly saves typing.
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
  \code{\link{prt}}, \code{\link{chainsaw2}},
  \code{\link[ape]{phylo}},
  \code{\link{get_indices_of_branches_under_tips}}
}

