\name{calc_AIC_column}
\alias{calc_AIC_column}
\title{Calculate AICs to make a column in a table}
\usage{
  calc_AIC_column(LnL_vals, nparam_vals)
}
\arguments{
  \item{LnL_vals}{A vector of log-likelihoods (typically
  negative, but may not be for continuous data).}

  \item{nparam_vals}{A vector of the number of parameters
  for each model.}
}
\value{
  \code{AIC_col} A \code{\link[base]{data.frame}} column of
  AIC results.
}
\description{
  A list of AICs (Akaike Information Criterion) is
  calculated from two input lists.  Lower values of AIC
  indicate some combination of better fit to the data and
  more parsimony in the model (fewer free parameters).
}
\details{
  The two input lists are:

  \bold{1.} A list of data likelihoods under a variety of
  models.\cr \bold{2.} A list of the number of free
  parameters under each model.\cr

  See \cite{Burnham et al. (2002)} and
  \url{http://www.brianomeara.info/tutorials/aic} for
  discussion of AIC and its uses.
}
\note{
  Go BEARS!
}
\examples{
test=1

LnL_vals = c(-20.9, -20.9, -20.9, -20.9)
nparam_vals = c(3, 4, 5, 6)
calc_AIC_column(LnL_vals, nparam_vals)
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
  \code{\link{calc_AIC_vals}}, \code{\link{calc_AICc_vals}}
}

