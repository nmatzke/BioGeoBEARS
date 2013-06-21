\name{given_a_starting_state_simulate_split}
\alias{given_a_starting_state_simulate_split}
\title{Given the state just below a node, simulate the states after speciation}
\usage{
  given_a_starting_state_simulate_split(index_Qmat_0based_of_starting_state = 1,
    COO_probs_columnar, numstates)
}
\arguments{
  \item{index_Qmat_0based_of_starting_state}{An integer
  index value, between 0 and \code{(numstates-1)}, which
  specifies what state is the starting point for the
  branch.}

  \item{COO_probs_columnar}{A speciation/cladogenesis
  transition matrix, in COO-like form, as produced by
  \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOweights_faster}}.}

  \item{numstates}{The number of states/geographic ranges.}
}
\value{
  \code{split_desc} 0-based indices of the descendant
  states in the two daughters.
}
\description{
  This function simulates a biogeographical history during
  a speciation/cladogenesis range inheritance event, given
  a cladogenesis probability transition matrix and a
  starting state.
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
  \url{https://code.google.com/p/lagrange/}

  Matzke_2012_IBS

  Matzke_Maguire_2011_SVP
}
\seealso{
  \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOweights_faster}},
  \code{\link{rcpp_calc_rowsums_for_COOweights_columnar}}
}

