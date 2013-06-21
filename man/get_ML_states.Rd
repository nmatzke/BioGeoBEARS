\name{get_ML_states}
\alias{get_ML_states}
\title{Get ML states from a BioGeoBEARS model results list}
\usage{
  get_ML_states(relprobs_matrix, unlist_TF = TRUE)
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
  results_object$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS}.}

  \item{unlist_TF}{Unlist the output? Default TRUE.}
}
\value{
  \code{inf_statesvec} The inferred vector of states.
}
\description{
  This function extracts the ML states from the results
  list produced by \code{\link{bears_2param_standard_fast}}
  or a similar ML search function.
}
\details{
  Currently, the scaled conditional probabilities are used
  to determine the optimum states.  However, this is not
  strictly correct, as these use only tips-down information
  (\cite{Felsenstein (2004)}; see also this post by Revell:
  \url{http://blog.phytools.org/2013/03/marginal-ancestral-state-reconstruction.html}).
  This is what \code{LAGRANGE} seems to do when reporting
  ancestral states, also (personal observation, perhaps
  imperfect, especially if the scaled conditional
  likelihoods and the marginal ancestral state
  probabilities turn out to be very close). What is desired
  is the marginal ancestral state reconstructions.  Most
  authors discuss ML ancestral state reconstruction as
  being a matter of re-rooting the tree at each node,
  yielding the marginal estimate for that node, conditional
  on the rest of the tree.  However, this procedure assumes
  a time-reversible model on both branches and cladogenesis
  events, and we have neither in biogeography.  Probably,
  the solution is just an up-pass from the root,
  calculating the probabilities on the forward model and
  multiplying by likelihoods from the downpass. However,
  this has not yet been implemented.
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

  Felsenstein2004

  Matzke_2012_IBS
}
\seealso{
  \code{\link{get_ML_probs}},
  \code{\link{bears_2param_standard_fast}},
  \code{\link{get_ML_state_indices}}
}

