\name{get_infprobs_of_simstates}
\alias{get_infprobs_of_simstates}
\title{Get the probabilities of the true (simulated) states}
\usage{
  get_infprobs_of_simstates(relprobs_matrix, simhist_row)
}
\arguments{
  \item{relprobs_matrix}{A relative probabilities matrix
  returned by \code{\link{bears_2param_standard_fast}} or a
  similar function. The user should specify WHICH matrix in
  the results_object -- i.e., scaled conditional
  likelihoods on downpass or uppass, or actual marginal
  probabilities of ancestral states.  (The latter is the
  main thing of interest.)  This specification is done via
  e.g. \code{relprobs_matrix =
  results_object$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS}..}

  \item{simhist_row}{A row from a table, which must have a
  column named \code{simulated_states_by_node_txt}.}
}
\value{
  \code{infprobs_of_simstates} The probability of each
  state at each node (all 1s and 0s).
}
\description{
  Basically this function assigns probability 1 to the
  simulated state/geographic range, and probability 0 for
  the other states/geographic ranges.  These data -- the
  simulated truth -- can then be compared to the inferred
  probabilities for the states, from e.g.
  \code{\link{get_ML_probs}}.
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
  \code{\link{simulate_biogeog_history}},
  \code{\link{infprobs_to_probs_of_each_area}}
}

