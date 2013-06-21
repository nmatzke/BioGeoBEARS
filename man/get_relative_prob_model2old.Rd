\name{get_relative_prob_model2old}
\alias{get_relative_prob_model2old}
\title{Calculate relative probability of model 1 (Akaike Weight)}
\usage{
  get_relative_prob_model2old(AICval_1, AICval_2)
}
\arguments{
  \item{AICval_1}{The AIC of the model of interest.}

  \item{AICval_2}{The AIC of another model of interest, for
  a pairwise comparison.}
}
\value{
  \code{relative_prob_model1} Akaike Weight of model 1.
}
\description{
  See \cite{Burnham et al. (2002)} and
  \url{http://www.brianomeara.info/tutorials/aic} for
  discussion of AIC and its uses.
}
\details{
  This is an older version of
  \code{\link{get_relative_prob_model1old}}, kept for
  back-compatibility.
}
\note{
  Go BEARS!
}
\examples{
test=1

AICval_1 = 20
AICval_2 = 30
get_relative_prob_model1old(AICval_1, AICval_2)
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
  \code{\link{AkaikeWeights_on_summary_table}},
  \code{\link{get_relative_prob_model1old}}
}

