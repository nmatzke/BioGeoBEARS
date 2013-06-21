\name{rel_likes_from_deltaAICs_pairwise}
\alias{rel_likes_from_deltaAICs_pairwise}
\title{Calculate the relative likelihoods of the models, from the deltaAICs, pairwise}
\usage{
  rel_likes_from_deltaAICs_pairwise(deltaAICs_pairwise)
}
\arguments{
  \item{deltaAICs_pairwise}{A vector of AIC values.}
}
\value{
  \code{rel_likes_AIC_pairwise} A
  \code{\link[base]{data.frame}} of relative likelihoods
  for each row (column 1) and the reference model (column
  2).
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
deltaAICs = get_deltaAIC_pairwise_w_ref_model(AICvals, ref_model="best")
deltaAICs

rel_likes_AIC_pairwise = rel_likes_from_deltaAICs_pairwise(deltaAICs)
rel_likes_AIC_pairwise
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

