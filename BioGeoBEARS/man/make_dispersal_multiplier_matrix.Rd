\name{make_dispersal_multiplier_matrix}
\alias{make_dispersal_multiplier_matrix}
\title{Make a default matrix of relative dispersal probabilities between areas}
\usage{
  make_dispersal_multiplier_matrix(areas = NULL,
    states_list = list("_", c("A"), c("B"), c("C"), c("A", "B"), c("B", "C"), c("A", "C"), c("A", "B", "C")))
}
\arguments{
  \item{areas}{A list of areas; if \code{NULL}, the states
  list will be used.}

  \item{states_list}{A list of states, where each state
  consists of a list of areas. A default example list is
  provided.}
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
make_dispersal_multiplier_matrix(areas=NULL, states_list=list("_", c("A"), c("B"), c("C"), c("A","B"), c("B","C"), c("A","C"), c("A","B","C")))
make_dispersal_multiplier_matrix(areas=c("A","B","C","D"))
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

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

  Foster PG (2001). "The Idiot's Guide to the Zen of
  Likelihood in a Nutshell in Seven Days for Dummies,
  Unleashed." Online PDF, widely copied, <URL:
  \url{http://www.bioinf.org/molsys/data/idiots.pdf}>.
}
\seealso{
  \code{\link{make_relprob_matrix_de}}
}

