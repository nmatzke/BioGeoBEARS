\name{get_relative_prob_model2old}
\alias{get_relative_prob_model2old}
\title{Calculate relative probability of model 1 (Akaike Weight)}
\usage{
  get_relative_prob_model2old(AICval_1, AICval_2)
}
\arguments{
  \item{AICval_1}{The AIC of the model of interest.}

  \item{AICval_2}{The AIC of another model of interest, for
  a pairwise comparison.}
}
\value{
  \code{relative_prob_model1} Akaike Weight of model 1.
}
\description{
  See \cite{Burnham et al. (2002)} and
  \url{http://www.brianomeara.info/tutorials/aic} for
  discussion of AIC and its uses.
}
\details{
  This is an older version of
  \code{\link{get_relative_prob_model1old}}, kept for
  back-compatibility.
}
\note{
  Go BEARS!
}
\examples{
AICval_1 = 20
AICval_2 = 30
get_relative_prob_model1old(AICval_1, AICval_2)
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
  \code{\link{AkaikeWeights_on_summary_table}},
  \code{\link{get_relative_prob_model1old}}
}

