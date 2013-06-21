\name{get_Akaike_weights_from_rel_likes}
\alias{get_Akaike_weights_from_rel_likes}
\title{Calculate the Akaike Weights, from the relative likelihoods of the models}
\usage{
  get_Akaike_weights_from_rel_likes(rel_likes_AIC)
}
\arguments{
  \item{rel_likes_AIC}{A vector of relative likelihoods.}
}
\value{
  \code{Akaike_weights} A vector of Akaike Weights.
}
\description{
  Given the relative likelihoods of the models, calculate
  the Akaike weight of the models. Akaike weights sum to 1.
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

Akaike_weights = rel_likes_from_deltaAICs(deltaAICs)
Akaike_weights
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

