\name{simulate_biogeog_history}
\alias{simulate_biogeog_history}
\title{Simulate a biogeographical history, given a transition matrix and cladogenesis model}
\usage{
  simulate_biogeog_history(phy, Qmat, COO_probs_columnar,
    index_Qmat_0based_of_starting_state)
}
\arguments{
  \item{phy}{An R \code{phylo} object.}

  \item{Qmat}{A (square, dense) Q transition matrix.  Using
  a sparse matrix would require writing another function.}

  \item{COO_probs_columnar}{A speciation/cladogenesis
  transition matrix, in COO-like form, as produced by
  \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOweights_faster}}.}

  \item{index_Qmat_0based_of_starting_state}{An integer
  index value, between 0 and \code{(numstates-1)}, which
  specifies what state will be the starting point for the
  simulation.}
}
\value{
  \code{simulated_states_by_node} A numeric matrix, giving
  the 0-based index of the state at each node and tip in
  the simulated history.  Getting a more detailed history
  would require a version of stochastic mapping
  (\cite{Huelsenbeck et al. (2003)}, \cite{Bollback
  (2005)}, \cite{Bollback (2006)}), but customized for the
  nonreversible and cladogenic aspects of biogeographical
  range evolution models.
}
\description{
  This function simulates a biogeographical history, given
  a Q transition matrix, a cladogenesis model giving the
  relative probability of different range inheritance
  scenarios, a phylogeny, and a 0-based index value
  deciding the starting state (which could be randomly
  generated according to a prior distribution of states).
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

  Huelsenbeck_etal_2003_stochastic_mapping

  Bollback_2005

  Bollback_2006_SIMMAP

  Matzke_2012_IBS
}
\seealso{
  \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOweights_faster}}
}

