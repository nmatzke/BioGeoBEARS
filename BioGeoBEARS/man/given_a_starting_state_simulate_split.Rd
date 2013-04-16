\name{given_a_starting_state_simulate_split}
\alias{given_a_starting_state_simulate_split}
\title{Given the state just below a node, simulate the states after speciation}
\usage{
  given_a_starting_state_simulate_split(index_Qmat_0based_of_starting_state = 1,
    COO_probs_columnar, numstates)
}
\arguments{
  \item{index_Qmat_0based_of_starting_state}{An integer
  index value, between 0 and \code{(numstates-1)}, which
  specifies what state is the starting point for the
  branch.}

  \item{COO_probs_columnar}{A speciation/cladogenesis
  transition matrix, in COO-like form, as produced by
  \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOweights_faster}}.}

  \item{numstates}{The number of states/geographic ranges.}
}
\value{
  \code{split_desc} 0-based indices of the descendant
  states in the two daughters.
}
\description{
  This function simulates a biogeographical history during
  a speciation/cladogenesis range inheritance event, given
  a cladogenesis probability transition matrix and a
  starting state.
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

  Matzke NJ and Maguire KC (2011). "Effect of fossil range
  and extinction data in DEC analyses: a case study with
  Canidae." _Society of Vertebrate Paleontology Annual
  Meeting, November 2011, Las Vegas_. <URL:
  http://www.vertpaleo.org/AM/Template.cfm?Section=Past_Meetings&Template=/CM/ContentDisplay.cfm&ContentID=2200>.
}
\seealso{
  \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOweights_faster}},
  \code{\link{rcpp_calc_rowsums_for_COOweights_columnar}}
}

