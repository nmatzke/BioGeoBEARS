\name{calcP_n}
\alias{calcP_n}
\title{Calculate Z (part of equation 6.4 of Harte 2011)}
\usage{
  calcP_n(n, lambda1, Z)
}
\arguments{
  \item{n}{Value of the state (e.g., which of a number of
  faces on a die, or number of different size classes of
  geographic range).}

  \item{lambda1}{Lambda parameter (\cite{Harte2011}).}

  \item{Z}{numeric values from \code{\link{calcZ_part}}.}
}
\value{
  \code{Prob_n}, numeric value of the probability of state
  \code{n}.
}
\description{
  This function is a used by \code{\link{get_probvals}},
  which calculates the Maximum Entropy (\cite{Harte
  (2011)}) discrete probability distribution of a number of
  ordered states (e.g., faces of a 6-sided die) given the
  mean of many rolls.  Here, this is merely used so that a
  single parameter can control the probability distribution
  of small versus large descendant areas during
  cladogenesis.
}
\details{
  See also: Maximum Entropy probability distribution for
  discrete variable with given mean (and discrete uniform
  flat prior)
  \url{http://en.wikipedia.org/wiki/Maximum_entropy_probability_distribution}
}
\note{
  Go BEARS!
}
\examples{
testval=1
n = 6
lambda1 = 0.5
Z = 1
calcP_n(n, lambda1, Z)
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
  \url{http://en.wikipedia.org/wiki/Maximum_entropy_probability_distribution}

  Harte J (2011). _Maximum entropy and ecology : a theory
  of abundance, distribution, and energetics_, series
  Oxford series in ecology and evolution. Oxford University
  Press, Oxford ; New York. ISBN 9780199593415 0199593418
  9780199593422 (pbk.) 0199593426 (pbk.)/ISSN, 2011284493
  John Harte.. ill. ; 24 cm. Includes bibliographical
  references (p. [244]-252) and index. Oxford series in
  ecology and evolution., <URL:
  http://socrates.berkeley.edu/~hartelab/MaxEnt.html>.

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
  \code{\link{calcZ_part}}, \code{\link[FD]{maxent}},
  \code{\link{symbolic_to_relprob_matrix_sp}}
}

