\name{get_indices_of_branches_under_tips}
\alias{get_indices_of_branches_under_tips}
\title{Get the indices of the branches (row number in edge matrix) below each tip}
\usage{
  get_indices_of_branches_under_tips(obj)
}
\arguments{
  \item{obj}{An \code{\link[ape]{ape}}
  \code{\link[ape]{phylo}} object}
}
\value{
  \code{branchnums_under_tips} The indices of the branches
  (row number in edge matrix) below each tip.
}
\description{
  A utility function. Gets the indices of the branches (row
  number in edge matrix) below each tip.
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
  \code{\link{get_indices_of_tip_nodes}},
  \code{\link{get_indices_where_list1_occurs_in_list2_noNA}}
}

