\name{check_BioGeoBEARS_run}
\alias{check_BioGeoBEARS_run}
\title{Check the inputs for various problems}
\usage{
  check_BioGeoBEARS_run(inputs, allow_huge_ranges = FALSE)
}
\arguments{
  \item{inputs}{The inputs list}

  \item{allow_huge_ranges}{Default FALSE, which will stop
  the run if there are more than 500 states. If TRUE, this
  will just print a warning, and continue, at which point
  you will wait for weeks or forever for the analysis to
  finish.  See \code{\link[cladoRcpp]{cladoRcpp}}'s
  \code{\link[cladoRcpp]{numstates_from_numareas}} function
  to calculate the size of the state space ahead of time,
  and links therein to see how the number of states scales
  with areas (2^number of areas, in an unconstrained
  analysis), how the size of the transition matrix you will
  be exponentiating scales (size = numstates * numstates),
  and the size of the
  ancestor/left-descendant/right-descendant cladogenesis
  matrix scales (numstates * numstates * numstates).  At
  500 states, this is 500^3 = 125,000,000 combinations of
  ancestor/left/right to check at every cladogenesis event,
  although \code{\link[cladoRcpp]{cladoRcpp}}'s tricks
  speed this up substantially.}
}
\value{
  \code{TRUE} if no errors found; otherwise a stop() is
  called.
}
\description{
  Numerous subtle mistakes in the input files for a
  BioGeoBEARS run can cause the run to crash.  As I come
  across these, I am putting in error checks for them.
}
\details{
  Some include:

  - Trees with negative branchlengths (as produced
  sometimes by e.g. BEAST MCC consensus trees (MCC =
  majority clade consensus).  These trees are always fully
  resolved, but the median node heights can sometimes be
  behind the node position in the tree.  Users should fix
  this manually, pathological results or crashes will
  result otherwise.

  - Trees with polytomies.  \code{BioGeoBEARS} (and
  \code{LAGRANGE}, and \code{DIVA}) assume a model where
  lineages bifurcate, and never multifurcate.  Users can
  convert multifurcating trees to bifurcating trees with
  \code{APE}'s \code{\link[ape]{multi2di}} (they will have
  to decide what branchlength to use for the new branches;
  it should be small, but bigger than the minimum
  branchlength used to identify fossils hooks (as hooks are
  considered to be anagenetic members of a lineage, and
  thus are connected to the tree without a cladogenesis
  event invoked). Users can then run their analysis several
  times on differently-resolved trees.

  NOTE: After the above correction, users may wish to
  correct the tip branchlengths (or make some other
  adjustment) so that all the tips are at age 0 my before
  present, as in an ultrametric tree. (However, note that
  trees with fossil tips are not ultrametric according to
  APE's \code{\link[ape]{is.ultrametric}}, even though they
  are time-scaled.  To make living (nonfossil) tips line up
  to zero, see \code{\link{average_tr_tips}} or the
  (different!) . They should be used with care.
  Alternately, a small amount of error in tip heights will
  make very little difference in the likelihood
  calculations (e.g. if some tips are 0.1 my too high, but
  the tree spans 200 my), which would be an argument for
  not requiring perfection after the (crucial) corrections
  of negative branchlengths, zero-branchlengths, and
  polytomies have been made.

  - Check for an absurdly large number of states.  I've set
  the limit at 500 (it starts getting slow around 200),
  users can override with \code{allow_huge_ranges=TRUE}.

  - Geography tipranges files should have same number of
  area labels as columns.

  - Geography tipranges files should have same number of
  taxa as the tree, and with the (exact!!) same names.
  This can be the source of many headaches, as different
  programs (\code{Mesquite}, etc.) treat spaces, periods,
  etc. in different ways, and re-write tipnames
  with/without quotes, underscores, etc.; and in my
  experience, my biologist colleagues find it very
  difficult to guarantee that the tipnames in their tree
  and their data tables will match exactly.  The SAFEST
  approach is to NEVER use these characters in tipnames or
  table names: space, comma, semicolon, dashes,
  parentheses, brackets, apostrophes or quote marks, or
  periods. Use ONLY letters, numbers, and underscores (_).
  When plotting trees, \code{APE} automatically reads
  underscores as spaces, which is nice for display.

  - There must be the same or more timeperiods than the
  other stratified items (distances matrices, etc.)
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

  Matzke_2012_IBS
}
\seealso{
  \code{\link{average_tr_tips}},
}

