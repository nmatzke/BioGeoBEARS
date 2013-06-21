\name{get_pruningwise_nodenums}
\alias{get_pruningwise_nodenums}
\title{Get internal node numbers in pruningwise order}
\usage{
  get_pruningwise_nodenums(tr)
}
\arguments{
  \item{tr}{A \code{\link[ape]{phylo}} tree object}
}
\value{
  \code{node_numbers_matrix} A matrix of node numbers
}
\description{
  There are many ways of numbering nodes in a tree.  This
  returns a matrix containing (column 1) R's native
  internal numbering scheme, and (column 2) the node
  numbers in a pruningwise downpass.  Note that this is
  different from \code{LAGRANGE}'s downpass ordering (see
  \code{\link{get_lagrange_nodenums}}).
}
\note{
  Go BEARS!
}
\examples{
extdata_dir = system.file("extdata", package="BioGeoBEARS")
trfn = slashslash(paste(extdata_dir,
"/examples/Psychotria_M0/LGcpp/Psychotria_5.2.newick", sep=""))
tr = read.tree(trfn)
node_numbers_matrix = get_pruningwise_nodenums(tr)
node_numbers_matrix
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link{get_lagrange_nodenums}}, \code{\link{prt}}
}

