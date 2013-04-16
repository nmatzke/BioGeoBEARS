\name{lrttest_on_summary_table}
\alias{lrttest_on_summary_table}
\title{Calculate Likelihood Ratio Test (LRT) results, and add to table}
\usage{
  lrttest_on_summary_table(restable, row_to_use_as_null,
    rows_to_exclude, returnwhat = "pval",
    add_to_table = TRUE)
}
\arguments{
  \item{restable}{A \code{\link[base]{data.frame}} with at
  least columns named "LnL" and "nparams".}

  \item{row_to_use_as_null}{This is the row specifying the
  model to which the others will be compared in pairwise
  fashion.}

  \item{rows_to_exclude}{Some rows may have models that the
  simpler model cannot nest within.  These should be
  excluded.}

  \item{returnwhat}{If "pval", just return the p-value.  If
  "all", return all of the intermediate outputs.}

  \item{add_to_table}{If TRUE, add to the main table and
  return the main table. If FALSE, return just the Akaike
  Weights results.}
}
\value{
  \code{pval} or \code{LRTrow}, both
  \code{\link[base]{data.frame}}.  Depends on
  \code{returnwhat}.
}
\description{
  The Likelihood Ratio Test (LRT) is a standard method for
  testing whether or not the data likelihood conferred by a
  more complex is significantly better than the data
  likelihood conferred by the simpler model, given a
  certain number of extra free parameters for the complex
  model.  The null hypothesis is that there is no
  difference; rejection means that there is a statistically
  significant improvement in the more complex model.
}
\details{
  The LRT only works for situations in which the simpler
  model is nested within the more complex model (i.e., by
  taking some parameters of the more complex model and
  forcing them to be fixed to a specific value).  In
  addition, the LRT may be unreliable in data-poor
  situations, and inherits whatever difficulties there may
  be in ML searches.  See \cite{Burnham et al. (2002)} for
  discussion.

  This function assumes that the log-likelihoods are in the
  column "LnL", and the number of parameters is specified
  in "nparams"
}
\note{
  Go BEARS!
}
\examples{
test=1
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

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
  \code{\link{lrttest}}
}

