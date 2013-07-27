\name{get_MLsplitprobs_from_results}
\alias{get_MLsplitprobs_from_results}
\title{Extract the ML probs for the base of each branch above a split}
\usage{
  get_MLsplitprobs_from_results(results_object)
}
\arguments{
  \item{results_object}{The results from a BioGeoBEARS ML
  search.}
}
\value{
  results_object with
  results_object$ML_marginal_prob_each_split_at_branch_bottom_BELOW_node
  added
}
\description{
  This function takes a BioGeoBEARS results_object from a
  ML search, extracts the downpass and uppass likelihoods
  of the data for each possible state at the base of each
  left and right branch, and produces the ML ancestral
  split estimates for the bottom of each branch.
}
\note{
  Go BEARS!
}
\examples{
test=1
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
  \url{https://code.google.com/p/lagrange/}

  Matzke_2012_IBS

  ReeSmith2008
}
\seealso{
  \code{\link{get_lagrange_nodenums}},
  \code{\link{LGpy_splits_fn_to_table}},
  \code{\link{LGcpp_splits_fn_to_table}}
}

