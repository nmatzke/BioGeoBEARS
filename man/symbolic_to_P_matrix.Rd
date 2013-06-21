\name{symbolic_to_P_matrix}
\alias{symbolic_to_P_matrix}
\title{Convert symbolic matrix to relprob matrix}
\usage{
  symbolic_to_P_matrix(dedf, cellsplit = "\\\\+",
    mergesym = "+", diags_sum_to_1 = FALSE, d = 0.1,
    e = 0.01, ...)
}
\arguments{
  \item{dedf}{The transition matrix or dispersal-extinction
  data.frame (dedf), contains the actual text of the
  formulas by which the transition probability matrix would
  be calculated.}

  \item{cellsplit}{The symbol to split the formulas on.
  Default "\\\\+" (plus symbol, with escape code).}

  \item{mergesym}{The symbol to merge the formulas with.
  Default "+".}

  \item{diags_sum_to_1}{Calculate the diagonals such that,
  when added to the sum of the off-diagonals in a row, the
  entire row sums to 1. This creates a transition
  probability matrix where each row sums to 1, i.e. each
  cell represents the conditional probability of the column
  state, given the ancestral row state.  The diagonal
  values represent the probability of staying the same.}

  \item{d}{The dispersal/range expansion rate. Default
  \code{d=0.1}.}

  \item{e}{The extinction/range contraction rate. Default
  \code{e=0.01}.}

  \item{...}{Additional arguments to pass to
  \code{\link{symbolic_cell_to_relprob_cell}} via
  \code{\link[base]{sapply}}, and thence to
  cell\code{\link[base]{strsplit}}.}
}
\value{
  \code{dedf_vals} The output
  \code{\link[base]{data.frame}}, contains the numeric
  results of the formulas calculating the transition
  probability matrix.
}
\description{
  This function takes a transition probability matrix (in
  text form) and converts to numeric form, given values for
  \emph{d}, \emph{e}, or other parameters in the text
  formulas.
}
\details{
  This is not particularly fast, but good for illustrative
  purposes.
}
\note{
  Go BEARS!
}
\examples{
testval=1

states_list = list("_", c("A"), c("B"), c("C"), c("A","B"),
c("B","C"), c("A","C"), c("A","B","C"))

states_list = areas_list_to_states_list_new(areas=c("A","B","C"),
include_null_range=TRUE, split_ABC=TRUE)
states_list

dedf = make_relprob_matrix_de(states_list=states_list,
 split_ABC=FALSE, split="", remove_simultaneous_events=TRUE,
add_multiple_Ds=TRUE,
dispersal_multiplier_matrix=make_dispersal_multiplier_matrix(states_list=states_list))
dedf

# Defaults
Pmat = symbolic_to_P_matrix(dedf, cellsplit="\\\\\\\\+", mergesym="+",
diags_sum_to_1=FALSE, d=0.1, e=0.01)
Pmat

# Calculate diagonal
Pmat = symbolic_to_P_matrix(dedf, cellsplit="\\\\\\\\+", mergesym="+",
diags_sum_to_1=TRUE, d=0.1, e=0.01)
Pmat

# You don't have to split, if the formulas are directly parsable
Pmat = symbolic_to_P_matrix(dedf, cellsplit="yadda", mergesym="",
diags_sum_to_1=FALSE, d=0.1, e=0.01)
Pmat
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link{areas_list_to_states_list_new}},
  \code{\link[cladoRcpp]{areas_list_to_states_list_old}},
  \code{\link{make_relprob_matrix_de}}
}

