\name{traverse_up}
\alias{traverse_up}
\title{Traverse the tree from node up to the tips}
\usage{
  traverse_up(tr4, startnode, traverse_records)
}
\arguments{
  \item{tr4}{A tree object in
  \code{\link[phylobase]{phylo4}} format.}

  \item{startnode}{The node number to start the uppass at.}

  \item{traverse_records}{A list of the nodes visited.}
}
\value{
  \code{traverse_records}
}
\description{
  This is a utility function for
  \code{\link{nodenums_bottom_up}}.
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
  \code{\link[phylobase]{phylo4}},
}

