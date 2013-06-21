\name{get_AICweight_ratio_model1_over_model2}
\alias{get_AICweight_ratio_model1_over_model2}
\title{Calculate ratio of Akaike Weights}
\usage{
  get_AICweight_ratio_model1_over_model2(AICval_1,
    AICval_2)
}
\arguments{
  \item{AICval_1}{The AIC of the model of interest.}

  \item{AICval_2}{The AIC of another model of interest, for
  a pairwise comparison.}
}
\value{
  \code{AICweight_ratio_model1} Ratio of Akaike Weights.
}
\description{
  See \cite{Burnham et al. (2002)} and
  \url{http://www.brianomeara.info/tutorials/aic} for
  discussion of AIC and its uses.
}
\note{
  Go BEARS!
}
\examples{
test=1

AICval_1 = 20
AICval_2 = 30
get_AICweight_ratio_model1_over_model2(AICval_1, AICval_2)
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
  \code{\link{AkaikeWeights_on_summary_table}}
}

