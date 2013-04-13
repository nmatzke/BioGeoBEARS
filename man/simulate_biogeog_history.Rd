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

  Huelsenbeck JP, Nielsen R and Bollback JP (2003).
  "Stochastic Mapping of Morphological Characters."
  _Systematic Biology_, *52*(2), pp. 131-158. <URL:
  http://dx.doi.org/10.1080/10635150390192780>, <URL:
  http://sysbio.oxfordjournals.org/content/52/2/131.abstract>.

  Bollback J (2005). "Posterior Mapping and Posterior
  Predictive Distributions." In series Statistics for
  Biology and Health, chapter 16, pp. 439-462. Springer New
  York. ISBN 978-0-387-22333-9/ISSN, <URL:
  http://dx.doi.org/10.1007/0-387-27733-1_16>, <URL:
  http://dx.doi.org/10.1007/0-387-27733-1_16>.

  Bollback J (2006). "SIMMAP: Stochastic character mapping
  of discrete traits on phylogenies." _BMC Bioinformatics_,
  *7*(1), pp. 88. <URL:
  http://www.biomedcentral.com/1471-2105/7/88>.

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
  \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOweights_faster}}
}

