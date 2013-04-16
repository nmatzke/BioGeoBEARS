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
tmptable = adf(c(40, 50, 60))
names(tmptable) = "AIC"
AkaikeWeights_and_Ratios_pairwise_on_summary_table_compared_to_ref(restable=tmptable, colname_to_use="AIC", ref_model="best", add_to_table=TRUE)
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
  \code{\link{get_Akaike_weights_from_rel_likes_pairwise}},
  \code{\link{get_Akaike_weights_from_rel_likes}},
  \code{\link{rel_likes_from_deltaAICs}},
  \code{\link{getAIC}}
}

