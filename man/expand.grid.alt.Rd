\name{expand.grid.alt}
\alias{expand.grid.alt}
\title{A faster version of expand.grid}
\usage{
  expand.grid.alt(seq1, seq2)
}
\arguments{
  \item{seq1}{A sequence of elements}

  \item{seq2}{A sequence of elements}
}
\value{
  \code{matrix_of_combinations} A matrix of all the
  possible combinations.
}
\description{
  This should be faster than
  \code{\link[base]{expand.grid}}, which "[c]reate[s] a
  data frame from all combinations of the supplied vectors
  or factors" (R documentation).
}
\details{
  The source of this function was this discussion thread:
  \url{http://stackoverflow.com/questions/10405637/use-outer-instead-of-expand-grid}
}
\note{
  Go BEARS!
}
\examples{
testval=1
seq1 = c("A","B","C")
seq2 = seq1
expand.grid(seq1,seq2)
expand.grid.alt(seq1,seq2)
expand.grid.jc(seq1,seq2)
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link[stats]{convolve}}, \code{\link{expand.grid}},
  \code{\link{expand.grid.jc}}
}

