\name{symbolic_cell_to_relprob_cell_sp}
\alias{symbolic_cell_to_relprob_cell_sp}
\title{Convert symbolic cell (a text equation) to relprob cell (a numeric value) -- speciation matrix version}
\usage{
  symbolic_cell_to_relprob_cell_sp(charcell,
    cellsplit = "\\\\+", mergesym = "*", ys = 1, j = 0,
    v = 1,
    relprob_subsets_matrix = relative_probabilities_of_subsets(6, 1e-04),
    relprob_vicar_matrix = relative_probabilities_of_vicariants(6, 1e-04),
    ...)
}
\arguments{
  \item{charcell}{The text formula.}

  \item{cellsplit}{The symbol to split the formulas on.
  Default "\\\\+" (plus symbol, with escape code).}

  \item{mergesym}{The symbol to merge the formulas with.
  Default "+".}

  \item{ys}{Relative weight of fully sympatric speciation
  (range-copying) and sympatric "subset" speciation.
  Default \code{s=1} mimics LAGRANGE model.}

  \item{v}{Relative weight of vicariant speciation. Default
  \code{v=1} mimics LAGRANGE model.}

  \item{j}{Relative weight of "founder event
  speciation"/jump speciation. Default \code{j=0} mimics
  LAGRANGE model.}

  \item{relprob_subsets_matrix}{A numeric matrix describing
  the relative probability of each smaller daughter range,
  conditional on the ancestral rangesize.}

  \item{relprob_vicar_matrix}{A numeric matrix describing
  the relative probability of each smaller daughter range,
  conditional on the ancestral rangesize.}

  \item{...}{Additional arguments to pass to
  \code{\link{relative_probabilities_of_subsets}} and
  \code{\link{relative_probabilities_of_vicariants}}, and
  thence to \code{\link[base]{strsplit}}.}
}
\value{
  \code{cellval} The output cell value.
}
\description{
  This does the equivalent of
  \code{\link{symbolic_to_P_matrix}}, but for a
  speciation/cladogenesis matrix.
}
\details{
  These are 1-event probability matrices, not instantaneous
  rate matrices.

  This function can be used in \code{\link[base]{sapply}}.
  It still will not be very fast compared to the
  calculations in \code{\link[cladoRcpp]{cladoRcpp}}, but
  can be useful for demonstrative purposes.
}
\note{
  Go BEARS!
}
\examples{
testval=1

charcell = "y1"
symbolic_cell_to_relprob_cell_sp(charcell, cellsplit="\\\\\\\\+", mergesym="*", ys=1,
j=0, v=1, relprob_subsets_matrix=relative_probabilities_of_subsets(max_numareas=3,
maxent_constraint_01=0.0001),
relprob_vicar_matrix=relative_probabilities_of_vicariants(max_numareas=3,
maxent_constraint_01v=0.0001))

charcell = "y1"
symbolic_cell_to_relprob_cell_sp(charcell, cellsplit="\\\\\\\\+", mergesym="*", ys=1,
j=1, v=1, relprob_subsets_matrix=relative_probabilities_of_subsets(max_numareas=3,
maxent_constraint_01=0.0001),
relprob_vicar_matrix=relative_probabilities_of_vicariants(max_numareas=3,
maxent_constraint_01v=0.0001))

charcell = "j"
symbolic_cell_to_relprob_cell_sp(charcell, cellsplit="\\\\\\\\+", mergesym="*", ys=1,
j=0, v=1, relprob_subsets_matrix=relative_probabilities_of_subsets(max_numareas=3,
maxent_constraint_01=0.0001),
relprob_vicar_matrix=relative_probabilities_of_vicariants(max_numareas=3,
maxent_constraint_01v=0.0001))

charcell = "j"
symbolic_cell_to_relprob_cell_sp(charcell, cellsplit="\\\\\\\\+", mergesym="*", ys=1,
j=1, v=1, relprob_subsets_matrix=relative_probabilities_of_subsets(max_numareas=3,
maxent_constraint_01=0.0001),
relprob_vicar_matrix=relative_probabilities_of_vicariants(max_numareas=3,
maxent_constraint_01v=0.0001))

charcell = "v1_2"
symbolic_cell_to_relprob_cell_sp(charcell, cellsplit="\\\\\\\\+", mergesym="*", ys=1,
j=0, v=1, relprob_subsets_matrix=relative_probabilities_of_subsets(max_numareas=3,
maxent_constraint_01=0.0001),
relprob_vicar_matrix=relative_probabilities_of_vicariants(max_numareas=3,
maxent_constraint_01v=0.0001))

charcell = "v1_2"
symbolic_cell_to_relprob_cell_sp(charcell, cellsplit="\\\\\\\\+", mergesym="*", ys=1,
 j=1, v=1, relprob_subsets_matrix=relative_probabilities_of_subsets(max_numareas=3,
maxent_constraint_01=0.0001),
relprob_vicar_matrix=relative_probabilities_of_vicariants(max_numareas=3,
maxent_constraint_01v=0.0001))

charcell = "s1_2"
symbolic_cell_to_relprob_cell_sp(charcell, cellsplit="\\\\\\\\+", mergesym="*", ys=1,
j=0, v=1, relprob_subsets_matrix=relative_probabilities_of_subsets(max_numareas=3,
maxent_constraint_01=0.0001),
relprob_vicar_matrix=relative_probabilities_of_vicariants(
max_numareas=3,
maxent_constraint_01v=0.0001))

charcell = "s1_2"
symbolic_cell_to_relprob_cell_sp(charcell, cellsplit="\\\\\\\\+", mergesym="*",
ys=1, j=1, v=1, relprob_subsets_matrix=relative_probabilities_of_subsets(
max_numareas=3, maxent_constraint_01=0.0001),
relprob_vicar_matrix=relative_probabilities_of_vicariants(
max_numareas=3,
maxent_constraint_01v=0.0001))
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link{symbolic_to_relprob_matrix_sp}},
  \code{\link{make_relprob_matrix_de}}
}

