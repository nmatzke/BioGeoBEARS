\name{chainsaw2}
\alias{chainsaw2}
\title{Saw a tree off at a particular time before present}
\usage{
  chainsaw2(tr, timepoint = 10, return_pieces = TRUE)
}
\arguments{
  \item{tr}{An ape phylo object.}

  \item{timepoint}{The time at which the tree should be
  "chopped".}

  \item{return_pieces}{Default TRUE, which means pieces
  should be returned}
}
\value{
  \code{chainsaw_result} (a list object with the pieces) or
  \code{tree_to_chainsaw}, just the leftover tree
}
\description{
  This function chops a tree like a hedge-trimmer, cutting
  straight across at a particular timepoint. The pieces are
  returned, as is the leftover tree, with branches
  shortened appropriately.  Pieces that are mini-trees are
  returned as ape objects, whereas single branches are just
  lengths.
}
\details{
  This function is used during stratification, but could
  have other uses as well.
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
  \code{\link{section_the_tree}}
}

