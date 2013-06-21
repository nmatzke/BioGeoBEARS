\name{symbolic_to_Q_matrix_exper}
\alias{symbolic_to_Q_matrix_exper}
\title{Experimental version of \code{symbolic_to_Q_matrix_exper}, including base frequencies}
\usage{
  symbolic_to_Q_matrix_exper(dedf, cellsplit = "\\\\+",
    mergesym = "*", d = 0.1, e = 0.01,
    basefreqs = rep(1, nrow(dedf))/nrow(dedf), ...)
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

  \item{d}{The dispersal/range expansion rate. Default
  \code{d=0.1}.}

  \item{e}{The extinction/range contraction rate. Default
  \code{e=0.01}.}

  \item{basefreqs}{Base frequencies, i.e. the equilibrium
  probabilities of the different states; the meaning of
  such an idea is debatable in the context of a
  LAGRANGE-like model where the null range (extinct
  everywhere) is included in the matrix and is a
  nonreversible absorbing state.  Default is
  \code{rep(1,nrow(dedf))/nrow(dedf)}.}

  \item{...}{Additional arguments to pass to
  \code{\link{symbolic_cell_to_relprob_cell}} via
  \code{\link[base]{sapply}}, and thence to
  cell\code{\link[base]{strsplit}}.}
}
\value{
  \code{dedf_vals} The output
  \code{\link[base]{data.frame}}, contains the Q matrix
}
\description{
  Still experimental.
}
\details{
  This function takes a transition probability matrix (in
  text form) and converts it to an instantaneous rate
  matrix (Q matrix), given values for \emph{d}, \emph{e},
  or other parameters in the text formulas.

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

dedf = make_relprob_matrix_de(states_list=states_list, split_ABC=FALSE,
split="", remove_simultaneous_events=TRUE, add_multiple_Ds=TRUE,
dispersal_multiplier_matrix=make_dispersal_multiplier_matrix(states_list=states_list))
dedf

# Right
Qmat = symbolic_to_Q_matrix_exper(dedf, cellsplit="\\\\\\\\+", mergesym="+", d=0.1, e=0.01)
Qmat

# Wrong
Qmat = symbolic_to_Q_matrix_exper(dedf, cellsplit="\\\\\\\\+", mergesym="*", d=0.1, e=0.01)
Qmat

# You don't have to split, if the formulas are directly parsable
Qmat = symbolic_to_Q_matrix_exper(dedf, cellsplit="yadda", mergesym="", d=0.1, e=0.01)
Qmat

# Compare to symbolic_to_Q_matrix
Qmat = symbolic_to_Q_matrix(dedf, cellsplit="yadda", mergesym="", d=0.1, e=0.01)
Qmat
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS

  FosterIdiots
}
\seealso{
  \code{\link[cladoRcpp]{areas_list_to_states_list_old}},
  \code{\link{make_relprob_matrix_de}}
  \code{\link{areas_list_to_states_list_new}}
}

