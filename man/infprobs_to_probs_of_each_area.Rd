\name{infprobs_to_probs_of_each_area}
\alias{infprobs_to_probs_of_each_area}
\title{Convert probabilities of each state, to the probabilities of presence in each area}
\usage{
  infprobs_to_probs_of_each_area(relprobs_matrix,
    states_list)
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

  \item{states_list}{A list of the possible
  states/geographic ranges, in 0-based index form.}
}
\value{
  \code{area_probs} The probability of presence in each
  area.
}
\description{
  Biogeographic inference in LAGRANGE and DIVA has focused
  heavily on inference of the exact ancestral
  state/geographic range.  However, when the state space is
  large, there is often considerable uncertainty in the
  exact ancestral range.  Even the ancestral state that
  confers the maximum likelihood on the data, and thus is
  the most probable ancestor, may have less than 50%
  probability, or even less (25%, 5%...), depending on the
  size of the state space.  This function converts the
  probability of specific states/geographic ranges into the
  probability of presence/absence in each area.  This can
  typically be inferred with much higher confidence.
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
  \code{\link{bears_2param_standard_fast}},
  \code{\link{get_ML_states}}, \code{\link{get_ML_probs}},
  \code{\link{infprobs_to_probs_of_each_area_from_relprobs}}
}

