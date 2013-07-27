\name{size_species_matrix}
\alias{size_species_matrix}
\title{Calculate the dimensions of the cladogenesis/speciation matrix}
\usage{
  size_species_matrix(states_list = default_states_list(),
    printwarn = 1)
}
\arguments{
  \item{states_list}{A list of states, where each state
  consists of a list of areas. A default example list is
  provided.}

  \item{printwarn}{If \code{printwarn>0}
  (\code{printwarn=1} by default), then print to screen a
  message describing the size of the cladogenesis matrix.}
}
\value{
  \code{spmat_dimensions} The dimensions of the
  cladogenesis matrix.
}
\description{
  This function calculates the dimensions of the
  cladogenesis/speciation matrix describing the transition
  probabilities between ancestral geographic ranges and
  descendant geographic range pairs on Left (L) and Right
  (R) branches.
}
\details{
  Under a cladogenesis model of geographic range change,
  the model will give the conditional probability of each
  possible combination of geographic ranges on the Left (L)
  and Right (R) descendant branches, conditional on a
  particular ancestral state.  A matrix representing these
  transitions will have \code{numstates} ancestral states,
  and \code{numstates*numstates} possible descendant pairs.
  Many of these will have 0 conditional probability under
  the model, but, for visualization or experimental
  purposes it can be useful to display them all.

  However, because \code{numstates = 2^numareas} under
  default conditions, and the number of cells the processor
  has to consider (without optimization tricks) is
  \code{numstates^3}, this transition matrix can very
  quickly become cumbersome to explicitly calculate or
  display.  \code{size_species_matrix} allows the user to
  check this ahead of time.

  See \code{\link[cladoRcpp]{numstates_from_numareas}} for
  the details of calculating \code{numstates}.

  At various points in \code{BioGeoBEARS} code, the text
  and numeric versions of the cladogenesis matrix are named
  \code{spmat} and \code{spPmat}, respectively.
}
\note{
  Go BEARS!
}
\examples{
testval=1
spmat_dimensions = size_species_matrix(
states_list=list("_", c("A"), c("B"), c("C"), c("A","B"),
c("B","C"), c("A","C"), c("A","B","C")), printwarn=1)
spmat_dimensions
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
  \code{\link{make_relprob_matrix_de}},
  \code{\link{make_spmat_row}}
}

