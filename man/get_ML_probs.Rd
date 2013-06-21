\name{get_ML_probs}
\alias{get_ML_probs}
\title{Get the probability of the ML state for each node, from a BioGeoBEARS model results list}
\usage{
  get_ML_probs(relprobs_matrix, unlist_TF = TRUE)
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

  \item{unlist_TF}{Unlist the output? Default TRUE.}
}
\value{
  \code{inf_probsvec} The inferred vector of probabilities
  of ML states.
}
\description{
  This function extracts the probability of the ML states
  from the results list produced by
  \code{\link{bears_2param_standard_fast}} or a similar ML
  search function.
}
\details{
  This is useful for displaying e.g. pie charts of the
  probability of the ML ancestral state at each node.

  Note, though, that it is somewhat peculiar and arbitrary
  to focus on the ancestral states just at nodes,
  particularly in the context of fossils with time ranges
  and geographic ranges.
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
  \url{http://blog.phytools.org/2013/03/marginal-ancestral-state-reconstruction.html}
  \url{https://code.google.com/p/lagrange/}

  Matzke_2012_IBS
}
\seealso{
  \code{\link{get_ML_probs}},
  \code{\link{bears_2param_standard_fast}},
  \code{\link{get_ML_state_indices}}
}

