\name{make_relprob_nummatrix_sp1}
\alias{make_relprob_nummatrix_sp1}
\title{Convert a observed-speciation transition matrix to an unobserved-speciation transition matrix (numeric version)}
\usage{
  make_relprob_nummatrix_sp1(probmat, spPmat,
    split = "\\\\|")
}
\arguments{
  \item{probmat}{A matrix of text, describing each of the
  allowed range-inheritance events.  Assumes that column
  names are in the "A|B" format.}

  \item{spPmat}{A matrix of numbers, where each cell
  contains the conditional probability of that
  ancestor-->(Left descendant,Right descendant) range
  inheritance scenario.}

  \item{split}{The value to split Left/Right pairs on
  (e.g., "A|B" --> "A", "B")}
}
\value{
  \code{newmat} A new square matrix.
}
\description{
  Convert a cladogenesis/speciation transition matrix
  (specifying the probability of each Left/Right descendant
  range pair, conditional on each ancestral state) of
  dimensions \code{numstates} by \code{numstates^2} to a
  square transition matrix of dimensions \code{numstates}
  by \code{numstates}, representing the probability of a
  transition when only one daughter survives in the tree.
}
\details{
  This matrix could be used to quantify the probability of
  range-change along a branch due to unobserved speciation
  events; all that would be required would be an estimate
  of the number of unobserved speciation events on the
  branch, and treating this as a Poisson process. (Note:
  this assumes that the probability of either branch
  surviving is identical, which might not be the case. See
  the GeoSSE (\cite{Goldberg et al. (2011)}) and ClaSSE
  (\cite{{"}Goldberg et al. (2012)}) for the beginnings of
  work on this, with 2 and 3 geographic areas,
  respectively.
}
\note{
  Go BEARS!
}
\examples{
testval=1
spmat = make_relprob_matrix_bi(states_list=list("_", c("A"),
c("B"), c("C"), c("A","B"), c("B","C"), c("A","C"), c("A","B","C")),
split_ABC=FALSE, splitval="", code_for_overlapping_subsets=NA, printwarn=1)
spmat

spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\\\\\\\+",
mergesym="*", ys=1, j=0, v=1, maxent_constraint_01=0.0001,
maxent_constraint_01v=0.0001, max_numareas=3)
spPmat

newmat = make_relprob_nummatrix_sp1(probmat=spmat, spPmat=spPmat, split="\\\\\\\\|")
newmat
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
  \url{https://code.google.com/p/lagrange/}
  \url{http://tigger.uic.edu/~eeg/code/code.html}

  Goldberg_etal_2011_GeoSSE

  Goldberg_Igic_2012_ClaSSE

  Matzke_2012_IBS

  ReeSmith2008
}
\seealso{
  \code{\link{make_relprob_matrix_bi}},
  \code{\link{make_relprob_txtmatrix_sp1}},
  \code{\link{paste_rows_without_zeros}}
}

