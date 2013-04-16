\name{meval}
\alias{meval}
\title{\code{eval()} function for use in \code{sapply}}
\usage{
  meval(equation_txt)
}
\arguments{
  \item{equation_txt}{The text of the equation to run
  \code{\link{eval}} on -- e.g., from a cell of a
  text-based transition matrix.}
}
\value{
  \code{outval} The numeric result of \code{\link{eval}}.
}
\description{
  \code{meval} is a wrapper for \code{\link{eval}}, to
  allow use in \code{sapply}.
}
\details{
  This is an attempt to speed up the use of
  \code{\link{eval}}; in general use of \code{\link{eval}}
  to convert a text version of a transition matrix to a
  numeric version with probabilities is a poor, slow
  choice; but it can be useful for examples and display
  purposes.

  See \code{\link[cladoRcpp]{cladoRcpp}} for fast C++
  implementations of transition matrix setup.
}
\note{
  Go BEARS!
}
\examples{
testval=1
d = 0.1
equation_txt = "1*d+1*d"
meval(equation_txt)
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
}
\seealso{
  \code{\link[stats]{convolve}}
}

