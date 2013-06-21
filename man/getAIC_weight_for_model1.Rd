\name{getAIC_weight_for_model1}
\alias{getAIC_weight_for_model1}
\title{Calculate Akaike Weight}
\usage{
  getAIC_weight_for_model1(AICval_1, AICvals)
}
\arguments{
  \item{AICval_1}{The AIC of the model of interest.}

  \item{AICvals}{The AICs of all the models being
  compared.}
}
\value{
  \code{AICweight} AICweight for the models.
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
AICvals = c(20,30,40)
getAIC_weight_for_model1(AICval_1, AICvals)
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

