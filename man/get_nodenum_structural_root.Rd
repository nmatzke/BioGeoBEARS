\name{get_nodenum_structural_root}
\alias{get_nodenum_structural_root}
\title{Gets the root node}
\usage{
  get_nodenum_structural_root(t, print_nodenum = FALSE)
}
\arguments{
  \item{t}{A tree object in \code{\link[ape]{phylo}}
  format.}

  \item{print_nodenum}{Print the node numbers as you go
  through the list? Default FALSE.}
}
\value{
  \code{root_nodenums_list}
}
\description{
  This function gets the root node by finding the node not
  in the descendants list (edge[,2]). This may be more
  reliable than e.g. assuming length(tr$tip.label)+1.
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

