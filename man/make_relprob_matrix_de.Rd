\name{make_relprob_matrix_de}
\alias{make_relprob_matrix_de}
\title{Make a relative dispersal probability matrix (in text form)}
\usage{
  make_relprob_matrix_de(states_list = default_states_list(),
    split_ABC = FALSE, split = "",
    remove_simultaneous_events = TRUE,
    add_multiple_Ds = TRUE,
    dispersal_multiplier_matrix = make_dispersal_multiplier_matrix(states_list))
}
\arguments{
  \item{states_list}{A list of states, where each state
  consists of a list of areas. A default example list is
  provided.}

  \item{split_ABC}{\code{TRUE} or \code{FALSE} If
  \code{TRUE} then each state/range in the input geographic
  ranges (\code{states_list}) will be split on the argument
  contained in \code{split}.}

  \item{split}{The character to split on.}

  \item{remove_simultaneous_events}{If \code{TRUE}
  (default, as in \code{LAGRANGE} and almost all
  phylogenetic Markov models), then it is assumed that all
  changes in geographic range along branches must happen
  one event at a time.  If \code{FALSE}, simultaneous
  events are not excluded; this is not recommended.
  However, notably, a commonly-used biogeographic model
  (treating biogeography as a multistate discrete character
  in an ML framework, where every species/lineage inhabits
  one and only one area at any point in time) effectively
  is invoking a simultaneous event: e.g., A->B is a
  simultaneous range gain and range loss, from the
  perspective of the dispersal-extinction framework.}

  \item{add_multiple_Ds}{If \code{TRUE} (default, as in
  \code{LAGRANGE}), the probabilities of dispersal from
  each possible source area are added together.}

  \item{dispersal_multiplier_matrix}{A user-provided
  dispersal multiplier matrix; the default is a matrix of
  1s from
  \code{\link{make_dispersal_multiplier_matrix}}(states_list=states_list).}
}
\value{
  \code{dedf} The output \code{\link[base]{data.frame}},
  termed \code{dedf} (dedf=dispersal-extinction
  data.frame), contains the actual text of the formulas by
  which the transition probability matrix would be
  calculated.
}
\description{
  This function takes a list of states/geographic ranges,
  and makes a relative probability matrix describing the
  probability of transition between each state.  These
  probabilities are described in terms of d, "dispersal"
  (actually range expansion) and "extinction" (actually
  local extirpation, or range contraction), as done in the
  program \code{LAGRANGE} (\cite{Ree et al. (2008)},
  \cite{Smith et al. (2010)}).
}
\details{
  The output \code{\link[base]{data.frame}}, termed
  \code{dedf} (dedf=dispersal-extinction data.frame),
  contains the actual text of the formulas by which the
  transition probability matrix would be calculated.  E.g.,
  the example calculates the matrix corresponding to
  Equation 1 on p. 6 of Ree & Smith (2008).

  Note that the geographic range-change process described
  here is a continuous-time process, where the probability
  of change is a function of branch length, and all
  transitions occur because of dispersal and extinction.
  LAGRANGE also implements a cladogenesis model (thus DEC
  -- dispersal-extinction-cladogenesis) which describes an
  "instantaneous" process of geographic range change at
  speciation/lineage-splitting events. \code{BioGeoBEARS}
  allows users to turn on, turn off, or otherwise customize
  both the continuous-time model and the cladogenesis
  model.
}
\note{
  Go BEARS!
}
\examples{
testval=1

states_list = list("_", c("A"), c("B"), c("C"), c
("A","B"), c("B","C"), c("A","C"), c("A","B","C"))

states_list = areas_list_to_states_list_new(
areas=c("A","B","C"), include_null_range=TRUE, split_ABC=TRUE)
states_list

dedf = make_relprob_matrix_de(states_list=states_list,
split_ABC=FALSE, split="", remove_simultaneous_events=TRUE,
add_multiple_Ds=TRUE,
dispersal_multiplier_matrix=make_dispersal_multiplier_matrix(states_list=states_list))

dedf
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
  \url{https://code.google.com/p/lagrange/}

  ReeSmith2008

  SmithRee2010_CPPversion

  Matzke_2012_IBS

  FosterIdiots
}
\seealso{
  \code{\link{make_dispersal_multiplier_matrix}}
}

