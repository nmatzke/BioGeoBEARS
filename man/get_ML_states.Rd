\name{get_ML_states}
\alias{get_ML_states}
\title{Get ML states from a BioGeoBEARS model results list}
\usage{
  get_ML_states(results_object)
}
\arguments{
  \item{results_object}{The results returned by
  \code{\link{bears_2param_standard_fast}} or a similar
  function.}
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

  Felsenstein J (2004). _Inferring Phylogenies_. Sinauer
  Associates, Sunderland, Massachusetts. ISBN 0878931775
  (pbk.)/ISSN, <URL:
  http://www.loc.gov/catdir/toc/ecip043/2003008942.html>.

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

