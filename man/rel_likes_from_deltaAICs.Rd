\name{rel_likes_from_deltaAICs}
\alias{rel_likes_from_deltaAICs}
\title{Calculate the relative likelihoods of the models, from the deltaAIC}
\usage{
  rel_likes_from_deltaAICs(deltaAICs)
}
\arguments{
  \item{deltaAICs}{A vector of deltaAIC values.}
}
\value{
  \code{rel_likes_AIC} A vector of relative likelihoods.
}
\description{
  Given deltaAIC (Akaike Information Criterion), the
  absolute difference between the best model (lowest AIC)
  and other models, calculate the relative likelihoods of
  the models.
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
test=1

AICvals = c(40, 50, 60)
deltaAICs = get_deltaAIC(AICvals)
deltaAICs

rel_likes_AIC = rel_likes_from_deltaAICs(deltaAICs)
rel_likes_AIC
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
  \code{\link{get_Akaike_weights_from_rel_likes}},
  \code{\link{rel_likes_from_deltaAICs}},
  \code{\link{getAIC}}
}

