\name{get_deltaAIC}
\alias{get_deltaAIC}
\title{Calculate deltaAIC}
\usage{
  get_deltaAIC(AICvals)
}
\arguments{
  \item{AICvals}{A vector of AIC values.}
}
\value{
  \code{deltaAICs} A vector of deltaAICs.
}
\description{
  Calculate deltaAIC (Akaike Information Criterion), the
  absolute difference between the best model (lowest AIC)
  and other models.
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
  \code{\link{rel_likes_from_deltaAICs}},
  \code{\link{getAIC}}
}

