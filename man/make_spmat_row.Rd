\name{make_spmat_row}
\alias{make_spmat_row}
\title{Construct a (text) cell of the cladogenesis/speciation matrix}
\usage{
  make_spmat_row(Lstates, Rstates, ancareas_txt_tmp,
    splitval = "", code_for_overlapping_subsets = NA)
}
\arguments{
  \item{Lstates}{A string listing the possible left states,
  which will be split by \code{splitval}.}

  \item{Rstates}{A string listing the possible right
  states, which will be split by \code{splitval}.}

  \item{ancareas_txt_tmp}{A string listing the possible
  ancestral states, which will be split by
  \code{splitval}.}

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
}
\value{
  \code{returncell} The text specifying the type of
  transition.
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

  Matzke_2012_IBS

  ReeSmith2008
}
\seealso{
  \code{\link{size_species_matrix}},
  \code{\link{make_relprob_matrix_bi}}
}

