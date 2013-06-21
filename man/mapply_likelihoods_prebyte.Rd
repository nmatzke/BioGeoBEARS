\name{mapply_likelihoods_prebyte}
\alias{mapply_likelihoods_prebyte}
\title{Use mapply on matrix exponentiations -- pre-byte-compiling}
\usage{
  mapply_likelihoods_prebyte(Qmat, phy2, transpose_needed)
}
\arguments{
  \item{Qmat}{an input Q transition matrix.}

  \item{phy2}{A phylogenetic tree.}

  \item{transpose_needed}{If TRUE (default), matrix will be
  transposed (apparently EXPOKIT needs the input matrix to
  be transposed compared to normal).}
}
\value{
  \code{independent_likelihoods_on_each_branch} The output
  matrix of the likelihoods for each state on each branch.
}
\description{
  During the likelihood calculations from the tips to the
  root of a tree, the transition matrix Qmat needs to be
  exponentiated for each branch length in the tree.  This
  is the slowest step of the likelihood calculation,
  especially for large matrices.  This function performs
  this with mapply.
}
\details{
  Byte-compiling is supposed to speed up functions; this is
  an attempt to do this on the
  \code{\link[rexpokit]{rexpokit}} function
  \code{\link[rexpokit]{expokit_dgpadm_Qmat}}.  It is also
  possible to byte-compile everything during package
  installation (via \code{ByteCompile: true} in the
  DESCRIPTION file), which is implemented in
  \code{BioGeoBEARS}, so this may be redundant.

  \code{\link{mapply_likelihoods_prebyte}} gets
  byte-compiled into \code{\link{mapply_likelihoods}}.

  See \url{http://dirk.eddelbuettel.com/blog/2011/04/12/}
  for discussion of the \code{\link[compiler]{compile}}
  package.
}
\note{
  Go BEARS!
}
\examples{
testval=1
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link[base]{mapply}},
  \code{\link[rexpokit]{expokit_dgpadm_Qmat}},
  \code{\link{expokit_dgpadm_Qmat2}},
  \code{\link[compiler]{compile}},
  \code{\link[compiler]{cmpfun}}
}

