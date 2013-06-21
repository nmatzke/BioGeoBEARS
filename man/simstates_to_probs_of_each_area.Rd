\name{simstates_to_probs_of_each_area}
\alias{simstates_to_probs_of_each_area}
\title{Convert simulated states to probabilities of each area}
\usage{
  simstates_to_probs_of_each_area(simulated_states_by_node,
    states_list, relprobs_matrix)
}
\arguments{
  \item{simulated_states_by_node}{The simulated states by
  node (0-based indices).}

  \item{states_list}{A list of the possible
  states/geographic ranges, in 0-based index form.}

  \item{relprobs_matrix}{A relative probabilities matrix
  returned by \code{\link{bears_2param_standard_fast}} or a
  similar function. The user should specify WHICH matrix in
  the results_object -- i.e., scaled conditional
  likelihoods on downpass or uppass, or actual marginal
  probabilities of ancestral states.  (The latter is the
  main thing of interest.)  This specification is done via
  e.g. \code{relprobs_matrix =
  results_object$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS}..}
}
\value{
  \code{area_probs} The probability of presence in each
  area.
}
\description{
  Basically this function assigns probability 1 to occupied
  areas according to the simulated state for a node, and
  probability 0 for the other areas.  These data -- the
  simulated truth -- can then be compared to the inferred
  probabilities of presence in each area, from
  \code{\link{infprobs_to_probs_of_each_area}}.
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

