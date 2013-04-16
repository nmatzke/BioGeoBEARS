\name{given_a_starting_state_simulate_branch_end}
\alias{given_a_starting_state_simulate_branch_end}
\title{Given the state at the start of a branch, simulate the state at the end of the branch}
\usage{
  given_a_starting_state_simulate_branch_end(index_Qmat_0based_of_starting_state = 1,
    Qmat, branchlength = 1, all_tips_living = TRUE)
}
\arguments{
  \item{index_Qmat_0based_of_starting_state}{An integer
  index value, between 0 and \code{(numstates-1)}, which
  specifies what state is the starting point for the
  branch.}

  \item{Qmat}{A (square, dense) Q transition matrix.  Using
  a sparse matrix would require writing another function.}

  \item{branchlength}{The length of the branch, or branch
  segment if you are dealing with a stratified phylogeny.}

  \item{all_tips_living}{Currently this is the only
  assumption.  If, hypothetically, you had a phylogeny with
  extinct tips (representing the ends of the ranges of
  fossil taxa), you might want to treat them differently,
  IF you think that the time-invariant geographic range
  addition/subtraction process is the same one that made
  lineages go extinct (it could be something else, e.g.
  mass extinction).  False attribution of extinctions to
  the range loss process will dramatically elevate the rate
  of range loss, and also range expansion to compensate,
  and the resulting high rates can substantially degrade
  inference (\cite{Matzke_Maguire_2011_SVP}).}
}
\value{
  \code{state_desc} 0-based index of the descendant state
  (just before cladogenesis, if below a node).
}
\description{
  This function simulates a biogeographical history, given
  a Q transition matrix, a starting state, and a branch
  length.  All this involves is exponentiating the Q
  transition matrix, producing a P transition probability
  matrix, and then producing a random draw from this P
  matrix, conditional on the ancestor.
}
\details{
  This could be sped up in various ways, if needed.
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
  \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOweights_faster}}
}

