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

  Harte2011

  Matzke_2012_IBS
}
\seealso{
  \code{\link{calcZ_part}}, \code{\link[FD]{maxent}},
  \code{\link{symbolic_to_relprob_matrix_sp}}
}

