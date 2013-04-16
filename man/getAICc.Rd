\name{getAICc}
\alias{getAICc}
\title{Calculate AICc}
\usage{
  getAICc(LnL, numparams, samplesize)
}
\arguments{
  \item{LnL}{The log-likelihood (typically negative, but
  may not be for continuous data).}

  \item{numparams}{The number of parameters for each
  model.}

  \item{samplesize}{The number of data on which the model
  conferred likelihood.}
}
\value{
  \code{AICcval} A vector of AICc results.
}
\description{
  Calculate AICc (Akaike Information Criterion).  Lower
  values of AICc indicate some combination of better fit to
  the data and more parsimony in the model (fewer free
  parameters).
}
\details{
  See \cite{Burnham et al. (2002)} and
  \url{http://www.brianomeara.info/tutorials/AICc} for
  discussion of AICc and its uses.
}
\note{
  Go BEARS!
}
\examples{
LnL = -34.5
numparams = 2
samplesize = 20
getAICc(LnL, numparams, samplesize)

LnL = -20.9
numparams = 3
samplesize = 20
getAICc(LnL, numparams, samplesize)

LnL = -34.5
numparams = 2
samplesize = 5
getAICc(LnL, numparams, samplesize)

LnL = -20.9
numparams = 3
samplesize = 5
getAICc(LnL, numparams, samplesize)
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
  \url{http://www.brianomeara.info/tutorials/AICc}

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
  \code{\link{calc_AICc_column}},
  \code{\link{calc_AICc_column}}
}

