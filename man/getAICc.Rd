\name{getAICc}
\alias{getAICc}
\title{Calculate AICc}
\usage{
  getAICc(LnL, numparams, samplesize)
}
\arguments{
  \item{LnL}{The log-likelihood (typically negative, but
  may not be for continuous data).}

  \item{numparams}{The number of parameters for each
  model.}

  \item{samplesize}{The number of data on which the model
  conferred likelihood.}
}
\value{
  \code{AICcval} A vector of AICc results.
}
\description{
  Calculate AICc (Akaike Information Criterion).  Lower
  values of AICc indicate some combination of better fit to
  the data and more parsimony in the model (fewer free
  parameters).
}
\details{
  See \cite{Burnham et al. (2002)} and
  \url{http://www.brianomeara.info/tutorials/AICc} for
  discussion of AICc and its uses.
}
\note{
  Go BEARS!
}
\examples{
LnL = -34.5
numparams = 2
samplesize = 20
getAICc(LnL, numparams, samplesize)

LnL = -20.9
numparams = 3
samplesize = 20
getAICc(LnL, numparams, samplesize)

LnL = -34.5
numparams = 2
samplesize = 5
getAICc(LnL, numparams, samplesize)

LnL = -20.9
numparams = 3
samplesize = 5
getAICc(LnL, numparams, samplesize)
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
  \url{http://www.brianomeara.info/tutorials/AICc}

  Burnham_Anderson_2002

  Matzke_2012_IBS
}
\seealso{
  \code{\link{calc_AICc_column}},
  \code{\link{calc_AICc_column}}
}

