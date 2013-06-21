\name{get_rownum_ref_model}
\alias{get_rownum_ref_model}
\title{Get rownum of named model}
\usage{
  get_rownum_ref_model(AICvals, ref_model = "best")
}
\arguments{
  \item{AICvals}{A vector of AIC values.}

  \item{ref_model}{What is the row of the reference model?
  "best", "worst", or a row number.}
}
\value{
  \code{ref_model_num} The
}
\description{
  Find the row number of the best model according to AIC,
  the worst model according to AIC, or just takes the row
  number if that is what was input.
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

  ReeSmith2008

  FosterIdiots
}
\seealso{
  \code{\link[stats]{convolve}}
}

