\name{AICstats_2models}
\alias{AICstats_2models}
\title{Calculate all the AIC and LRT stats between two models}
\usage{
  AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
}
\arguments{
  \item{LnL_1}{Log-likelihood of more complex model.}

  \item{LnL_2}{Log-likelihood of simpler complex model.}

  \item{numparams1}{Number of free parameters of the more
  complex model.}

  \item{numparams2}{Number of free parameters of the less
  complex model.}
}
\value{
  \code{LRT_AIC_results} A table of LRT and AIC results.
}
\description{
  The Likelihood Ratio Test (LRT) is a standard method for
  testing whether or not the data likelihood conferred by a
  more complex is significantly better than the data
  likelihood conferred by the simpler model.  See
  \code{\link{lrttest}} and
  \code{\link{lrttest_on_summary_table}} for more
  discussion.
}
\details{
  See \cite{Burnham et al. (2002)} and
  \url{http://www.brianomeara.info/tutorials/aic} for
  discussion of AIC and its uses.

  This function assumes that \code{LnL_1} and
  \code{numparams1} refer to the more complex model, and
  that \code{LnL_2} and \code{numparams2} refer to the
  simpler model nested within the more complex one.
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

  Burnham_Anderson_2002

  Matzke_2012_IBS
}
\seealso{
  \code{\link{lrttest}},
  \code{\link{lrttest_on_summary_table}}
}

