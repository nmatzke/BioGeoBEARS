\name{get_Akaike_weights_from_rel_likes_pairwise}
\alias{get_Akaike_weights_from_rel_likes_pairwise}
\title{Calculate the Akaike Weights, from the relative likelihoods of the models}
\usage{
  get_Akaike_weights_from_rel_likes_pairwise(rel_likes_AIC_pairwise)
}
\arguments{
  \item{rel_likes_AIC_pairwise}{A 2-column
  \code{data.frame} of relative likelihoods of each pair of
  models.}
}
\value{
  \code{Akaike_weights_pairwise} A
  \code{\link[base]{data.frame}} of Akaike Weights for each
  row (column 1) and the reference model (column 2). Note
  that only 2 models are being compared in each row, not
  all of them, as in
  \code{\link{get_Akaike_weights_from_rel_likes}}.
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
deltaAICs = get_deltaAIC_pairwise_w_ref_model(AICvals, ref_model="best")
deltaAICs

rel_likes_AIC_pairwise = rel_likes_from_deltaAICs_pairwise(deltaAICs)
rel_likes_AIC_pairwise

Akaike_weights_pairwise = get_Akaike_weights_from_rel_likes_pairwise(rel_likes_AIC_pairwise)
Akaike_weights_pairwise
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

