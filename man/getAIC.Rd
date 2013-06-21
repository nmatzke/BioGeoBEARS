\name{getAIC}
\alias{getAIC}
\title{Calculate AIC}
\usage{
  getAIC(LnL, numparams)
}
\arguments{
  \item{LnL}{The log-likelihood (typically negative, but
  may not be for continuous data).}

  \item{numparams}{The number of parameters for each
  model.}
}
\value{
  \code{AICval} A vector of AIC results.
}
\description{
  Calculate AIC (Akaike Information Criterion).  Lower
  values of AIC indicate some combination of better fit to
  the data and more parsimony in the model (fewer free
  parameters).
}
\details{
  See \cite{Burnham et al. (2002)} and
  \url{http://www.brianomeara.info/tutorials/aic} for
  discussion of AIC and its uses.
}
\note{
  Go BEARS!
}
\examples{
LnL = -34.5
numparams = 2
getAIC(LnL, numparams)

LnL = -20.9
numparams = 3
getAIC(LnL, numparams)

# It turns out to work on lists, also
LnL = c(-34.5, -20.9)
numparams = c(2, 3)
getAIC(LnL, numparams)
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
  \url{http://www.brianomeara.info/tutorials/aic}

  Burnham_Anderson_2002

  Matzke_2012_IBS
}
\seealso{
  \code{\link{calc_AIC_column}},
  \code{\link{calc_AIC_column}}
}

