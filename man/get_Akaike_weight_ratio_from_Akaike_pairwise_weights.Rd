\name{get_Akaike_weight_ratio_from_Akaike_pairwise_weights}
\alias{get_Akaike_weight_ratio_from_Akaike_pairwise_weights}
\title{Get the ratio between the pairwise Akaike Weights}
\usage{
  get_Akaike_weight_ratio_from_Akaike_pairwise_weights(Akaike_weights_pairwise)
}
\arguments{
  \item{Akaike_weights_pairwise}{A 2-column
  \code{data.frame} of Akaike Weights for each pair of
  models.}
}
\value{
  \code{Akaike_weight_ratios_pairwise} A
  \code{\link[base]{data.frame}} of Akaike Weight Ratios
  for each row (column 1) and the reference model (column
  2). Note that only 2 models are being compared in each
  row, not all of them.
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

Akaike_weights_pairwise = get_Akaike_weights_from_rel_likes_pairwise(
rel_likes_AIC_pairwise)
Akaike_weights_pairwise

Akaike_weight_ratios_pairwise = get_Akaike_weight_ratio_from_Akaike_pairwise_weights(
Akaike_weights_pairwise)
Akaike_weight_ratios_pairwise
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
  \code{\link{get_Akaike_weights_from_rel_likes_pairwise}},
  \code{\link{get_Akaike_weights_from_rel_likes}},
  \code{\link{rel_likes_from_deltaAICs}},
  \code{\link{getAIC}}
}

