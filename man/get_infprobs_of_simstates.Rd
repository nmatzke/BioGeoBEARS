\name{get_infprobs_of_simstates}
\alias{get_infprobs_of_simstates}
\title{Get the probabilities of the true (simulated) states}
\usage{
  get_infprobs_of_simstates(results_object, simhist_row)
}
\arguments{
  \item{results_object}{The results returned by
  \code{\link{bears_2param_standard_fast}} or a similar
  function.}

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
  \code{\link{simulate_biogeog_history}},
  \code{\link{infprobs_to_probs_of_each_area}}
}

