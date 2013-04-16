\name{getAIC}
\alias{getAIC}
\title{Calculate AIC}
\usage{
  getAIC(LnL, numparams)
}
\arguments{
  \item{LnL}{The log-likelihood (typically negative, but
  may not be for continuous data).}

  \item{numparams}{The number of parameters for each
  model.}
}
\value{
  \code{AICval} A vector of AIC results.
}
\description{
  Calculate AIC (Akaike Information Criterion).  Lower
  values of AIC indicate some combination of better fit to
  the data and more parsimony in the model (fewer free
  parameters).
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
LnL = -34.5
numparams = 2
getAIC(LnL, numparams)

LnL = -20.9
numparams = 3
getAIC(LnL, numparams)

# It turns out to work on lists, also
LnL = c(-34.5, -20.9)
numparams = c(2, 3)
getAIC(LnL, numparams)
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
  \code{\link{calc_AIC_column}},
  \code{\link{calc_AIC_column}}
}

