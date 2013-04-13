\name{get_ML_probs}
\alias{get_ML_probs}
\title{Get the probability of the ML state for each node, from a BioGeoBEARS model results list}
\usage{
  get_ML_probs(results_object)
}
\arguments{
  \item{results_object}{The results returned by
  \code{\link{bears_2param_standard_fast}} or a similar
  function.}
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

  Matzke N (2012). "Founder-event speciation in BioGeoBEARS
  package dramatically improves likelihoods and alters
  parameter inference in Dispersal-Extinction-Cladogenesis
  (DEC) analyses." _Frontiers of Biogeography_, *4*(suppl.
  1), pp. 210. ISSN 1948-6596, Poster abstract published in
  the Conference Program and Abstracts of the International
  Biogeography Society 6th Biannual Meeting, Miami,
  Florida. Poster Session P10: Historical and
  Paleo-Biogeography. Poster 129B. January 11, 2013, <URL:
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}>.
}
\seealso{
  \code{\link{get_ML_probs}},
  \code{\link{bears_2param_standard_fast}},
  \code{\link{get_ML_state_indices}}
}

