\name{calc_AICc_vals}
\alias{calc_AICc_vals}
\title{Calculate AICc values for a list of models}
\usage{
  calc_AICc_vals(LnL_vals, nparam_vals, samplesize)
}
\arguments{
  \item{LnL_vals}{A vector of log-likelihoods (typically
  negative, but may not be for continuous data).}

  \item{nparam_vals}{A vector of the number of parameters
  for each model.}

  \item{samplesize}{A single samplesize, or a vector of the
  samplesizes each model.  However, samplesize should
  always be the same for all comparisons, since maximum
  likelihood and AIC/AICc model-selection methods are
  always comparing different models on the \emph{same}
  data, not different data on the same mode.}
}
\value{
  \code{AICc_vals} A vector of AICc results.
}
\description{
  A list of AICc values (second order Akaike Information
  Criterion) is calculated from two input lists.  Lower
  values of AICc indicate some combination of better fit to
  the data and more parsimony in the model (fewer free
  parameters). AICc contains a correction for sample size.
}
\details{
  The two input lists are:

  \bold{1.} A list of data likelihoods under a variety of
  models.\cr \bold{2.} A list of the number of free
  parameters under each model.\cr

  \code{samplesize} can be a scalar or vector; but see
  below.

  See \cite{Burnham et al. (2002)} and
  \url{http://www.brianomeara.info/tutorials/aic} for
  discussion of AIC, AICc and their uses.
}
\note{
  Go BEARS!
}
\examples{
LnL_vals = c(-34.5, -20.9)
nparam_vals = c(2, 3)
calc_AICc_vals(LnL_vals, nparam_vals, samplesize=20)

LnL_vals = c(-20.9, -20.9, -20.9, -20.9)
nparam_vals = c(3, 4, 5, 6)
calc_AICc_vals(LnL_vals, nparam_vals, samplesize=20)
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
  \code{\link{calc_AIC_vals}},
  \code{\link{calc_AICc_column}}
}

