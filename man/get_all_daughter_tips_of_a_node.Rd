\name{get_all_daughter_tips_of_a_node}
\alias{get_all_daughter_tips_of_a_node}
\title{Get all the daughter tips of a node}
\usage{
  get_all_daughter_tips_of_a_node(nodenum, t)
}
\arguments{
  \item{nodenum}{The node to find}

  \item{t}{A \code{\link[ape]{phylo}} tree object.}
}
\value{
  \code{temp_tips} The list of daughter tipnodes
}
\description{
  Like it says. Utility function.
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
  \code{\link{add_to_downpass_labels}},
  \code{\link[ape]{extract.clade}}
}

