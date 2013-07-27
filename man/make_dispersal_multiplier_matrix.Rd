\name{make_dispersal_multiplier_matrix}
\alias{make_dispersal_multiplier_matrix}
\title{Make a default matrix of relative dispersal probabilities between areas}
\usage{
  make_dispersal_multiplier_matrix(areas = NULL,
    states_list = default_states_list(),
    dispersal_multipliers_matrix = NULL,
    distances_mat = NULL, x_exponent = 0)
}
\arguments{
  \item{areas}{A list of areas; if \code{NULL}, the states
  list will be used.}

  \item{states_list}{A list of states, where each state
  consists of a list of areas. A default example list is
  provided.}

  \item{dispersal_multipliers_matrix}{Default NULL.}

  \item{distances_mat}{Default NULL.}

  \item{x_exponent}{Default 0.}
}
\value{
  \code{dispersal_multiplier_matrix} A square matrix, with
  1s for all cells.
}
\description{
  Given either a list of areas, or a list of states, this
  function provides a square dispersal matrix giving the
  relative probability of dispersal between areas.  The
  function fills in these dispersals probabilities with the
  value 1.  The user can then modify this as desired.
  dispersal_multipliers_matrix Default NULL distances_mat
  Default NULL x_exponent Default 0
}
\details{
  If only a states list is given, the list of areas is
  calculated by getting \code{\link[base]{unique}} values
  from the concatenated states list.
}
\note{
  Go BEARS!
}
\examples{
testval=1
make_dispersal_multiplier_matrix(areas=NULL,
states_list=list("_", c("A"), c("B"), c("C"),
c("A","B"), c("B","C"), c("A","C"), c("A","B","C")))
make_dispersal_multiplier_matrix(areas=c("A","B","C","D"))
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS

  FosterIdiots
}
\seealso{
  \code{\link{make_relprob_matrix_de}}
}

