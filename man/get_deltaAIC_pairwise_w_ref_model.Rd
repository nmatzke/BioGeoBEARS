\name{get_deltaAIC_pairwise_w_ref_model}
\alias{get_deltaAIC_pairwise_w_ref_model}
\title{Calculate deltaAIC}
\usage{
  get_deltaAIC_pairwise_w_ref_model(AICvals,
    ref_model = "best")
}
\arguments{
  \item{AICvals}{A vector of AIC values.}

  \item{ref_model}{What is the row of the reference model?
  "best", "worst", or a row number.}
}
\value{
  \code{deltaAICs_pairwise} A 2-column
  \code{\link[base]{data.frame}} of pairwise deltaAICs for
  each row (column 1) and the reference model (column 2).
}
\description{
  Calculate deltaAIC (Akaike Information Criterion), the
  absolute difference between the best model (lowest AIC)
  and other models.  This function does it pairwise only,
  with a reference model.
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
get_deltaAIC(AICvals)
get_deltaAIC_pairwise_w_ref_model(AICvals, ref_model="best")
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
  \code{\link{get_deltaAIC}},
  \code{\link{rel_likes_from_deltaAICs}},
  \code{\link{getAIC}}
}

