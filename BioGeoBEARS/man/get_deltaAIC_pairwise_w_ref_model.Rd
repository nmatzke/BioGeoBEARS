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

  Burnham KP and Anderson DR (2002). _Model selection and
  multimodel inference : a practical information-theoretic
  approach_, 2nd edition. Springer, New York. ISBN
  0387953647 (alk. paper)/ISSN, 2001057677 Kenneth P.
  Burnham, David R. Anderson. ill. ; 24 cm. Rev. ed. of:
  Model selection and inference. c1998. Includes
  bibliographical references (p. [455]-484) and index.,
  <URL:
  http://www.loc.gov/catdir/enhancements/fy0816/2001057677-d.htmlhttp://www.loc.gov/catdir/enhancements/fy0816/2001057677-t.html>.

  Matzke N (2012). "Founder-event speciation in BioGeoBEARS
  package dramatically improves likelihoods and alters
  parameter inference in Dispersal-Extinction-Cladogenesis
  (DEC) analyses." _Frontiers of Biogeography_, *4*(suppl.
  1), pp. 210. ISSN 1948-6596, Poster abstract published in
  the Conference Program and Abstracts of the International
  Biogeography Society 6th Biannual Meeting, Miami,
  Florida. Poster Session P10: Historical and
  Paleo-Biogeography. Poster 129B. January 11, 2013, <URL:
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}>.
}
\seealso{
  \code{\link{get_deltaAIC}},
  \code{\link{rel_likes_from_deltaAICs}},
  \code{\link{getAIC}}
}

