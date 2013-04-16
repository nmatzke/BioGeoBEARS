\name{calc_AIC_column}
\alias{calc_AIC_column}
\title{Calculate AICs to make a column in a table}
\usage{
  calc_AIC_column(LnL_vals, nparam_vals)
}
\arguments{
  \item{LnL_vals}{A vector of log-likelihoods (typically
  negative, but may not be for continuous data).}

  \item{nparam_vals}{A vector of the number of parameters
  for each model.}
}
\value{
  \code{AIC_col} A \code{\link[base]{data.frame}} column of
  AIC results.
}
\description{
  A list of AICs (Akaike Information Criterion) is
  calculated from two input lists.  Lower values of AIC
  indicate some combination of better fit to the data and
  more parsimony in the model (fewer free parameters).
}
\details{
  The two input lists are:

  \bold{1.} A list of data likelihoods under a variety of
  models.\cr \bold{2.} A list of the number of free
  parameters under each model.\cr

  See \cite{Burnham et al. (2002)} and
  \url{http://www.brianomeara.info/tutorials/aic} for
  discussion of AIC and its uses.
}
\note{
  Go BEARS!
}
\examples{
LnL_vals = c(-20.9, -20.9, -20.9, -20.9)
nparam_vals = c(3, 4, 5, 6)
calc_AIC_column(LnL_vals, nparam_vals)
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
  \code{\link{calc_AIC_vals}}, \code{\link{calc_AICc_vals}}
}

