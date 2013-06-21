\name{make_all_tips_come_to_0}
\alias{make_all_tips_come_to_0}
\title{Take a tree, ensure all tips end at 0}
\usage{
  make_all_tips_come_to_0(trfn, outfn)
}
\arguments{
  \item{trfn}{The Newick tree filename}

  \item{outfn}{The filename for the resulting file}
}
\value{
  \code{phy} The corrected phylogeny
}
\description{
  Makes tree precisely ultrametric (unlike e.g. the default
  <i>Psychotria</i> tree)
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
  \url{https://code.google.com/p/lagrange/}

  Matzke_2012_IBS
}
\seealso{
  \code{\link[ape]{read.tree}}
}

