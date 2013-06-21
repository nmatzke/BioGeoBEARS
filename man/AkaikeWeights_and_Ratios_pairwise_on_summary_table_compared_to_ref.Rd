\name{AkaikeWeights_and_Ratios_pairwise_on_summary_table_compared_to_ref}
\alias{AkaikeWeights_and_Ratios_pairwise_on_summary_table_compared_to_ref}
\title{Get the ratio between the pairwise Akaike Weights}
\usage{
  AkaikeWeights_and_Ratios_pairwise_on_summary_table_compared_to_ref(restable,
    colname_to_use = "AIC", ref_model = "best",
    add_to_table = TRUE)
}
\arguments{
  \item{restable}{A \code{\link[base]{data.frame}} with at
  least columns named "LnL" and "nparams".}

  \item{colname_to_use}{The name of the column containing
  AIC values.}

  \item{ref_model}{What is the row of the reference model?
  "best", "worst", or a row number.}

  \item{add_to_table}{If TRUE, add to the main table and
  return the main table. If FALSE, return just the Akaike
  Weights results.}
}
\value{
  \code{restable}, the modified table, or
  \code{AICstats_pairwise}, the pairwise Akaike statistics.
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

tmptable = adf(c(40, 50, 60))
names(tmptable) = "AIC"
AkaikeWeights_and_Ratios_pairwise_on_summary_table_compared_to_ref(
restable=tmptable, colname_to_use="AIC", ref_model="best", add_to_table=TRUE)
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

