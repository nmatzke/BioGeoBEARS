\name{make_relprob_matrix_bi}
\alias{make_relprob_matrix_bi}
\title{Make a relative probability matrix for a single speciation (bifurcation) event}
\usage{
  make_relprob_matrix_bi(states_list = default_states_list(),
    split_ABC = FALSE, splitval = "",
    code_for_overlapping_subsets = NA, printwarn = 1)
}
\arguments{
  \item{states_list}{A list of states, where each state
  consists of a list of areas. A default example list is
  provided.}

  \item{split_ABC}{\code{TRUE} or \code{FALSE} If
  \code{TRUE} then each state/range in the input geographic
  ranges (\code{states_list}) will be split on the argument
  contained in \code{split}.}

  \item{splitval}{The character to split on.}

  \item{code_for_overlapping_subsets}{Hypothetically, there
  is no reason that a vicariance event could happen, e.g.
  ABC-->AB, BC.  This is disallowed in \code{LAGRANGE}
  BioGeoBEARS defaults, and, if one is going to employ the
  construct of discrete areas in the first place, overlaps
  should probably be avoided.  But this parameter will
  allow experimentation.  Here,
  \code{code_for_overlapping_subsets=NA} equals the
  default, and any other value means that overlapping
  vicariance events are included, with a number describing
  the number of areas in the overlap.  Users could then
  manually convert this to a probability according to some
  function.}

  \item{printwarn}{If \code{printwarn>0}
  (\code{printwarn=1} by default), then print to screen a
  message describing the size of the cladogenesis matrix.}
}
\value{
  \code{probmat} A matrix of strings, where each cell
  contains the parameters describing the conditional
  probability of that ancestor-->(Left descendant,Right
  descendant) range inheritance scenario.
}
\description{
  Given the identity of the states/geographic ranges on the
  left branch (\code{Lstates}), right branch
  (\code{Rstates}), and ancestral areas
  (\code{ancareas_txt_tmp}), construct the (text version)
  of the row of transition probabilities.  This means that
  each nonzero cell gets a \emph{v} for a vicariance event,
  a \emph{y} for a sympatric speciation/range-copying
  event, a \emph{j} for a founder-event/jump speciation
  event, and an \emph{s} for a sympatric-subset event.
}
\details{
  This function is utilized by \code{\link[base]{apply}} in
  other functions (e.g. ) in an attempt to speed up
  calculation over rows.  However, processing of text
  formulas via \code{\link[base]{apply}} will never be fast
  enough for large matrices; see
  \code{\link[cladoRcpp]{cladoRcpp}} for optimized
  functions.

  This text-based matrix later gets evaluated by other
  functions to calculate the numerical probabilities.
  I.e., if j=0 and the other forms of speciation have
  weights equal to each other, this is the \code{LAGRANGE}
  cladogenesis model.

  NOTE: This function is veeeeeeery slow, even for only 3
  areas (i.e. \code{2^3=8} geographic ranges).  It is
  mostly useful for illustration. See
  \code{\link[cladoRcpp]{cladoRcpp}} for drastic
  improvements in calculating cladogenesis models.
}
\note{
  Go BEARS!
}
\examples{
testval=1
probmat = make_relprob_matrix_bi(states_list=list("_",
c("A"), c("B"), c("C"), c("A","B"), c("B","C"), c("A","C"),
c("A","B","C")), split_ABC=FALSE, splitval="",
code_for_overlapping_subsets=NA, printwarn=1)
probmat
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
  \code{\link{size_species_matrix}},
  \code{\link{make_spmat_row}}
}

