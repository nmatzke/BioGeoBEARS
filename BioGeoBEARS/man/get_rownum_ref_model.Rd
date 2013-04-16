\name{get_rownum_ref_model}
\alias{get_rownum_ref_model}
\title{Get rownum of named model}
\usage{
  get_rownum_ref_model(AICvals, ref_model = "best")
}
\arguments{
  \item{AICvals}{A vector of AIC values.}

  \item{ref_model}{What is the row of the reference model?
  "best", "worst", or a row number.}
}
\value{
  \code{ref_model_num} The
}
\description{
  Find the row number of the best model according to AIC,
  the worst model according to AIC, or just takes the row
  number if that is what was input.
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
  \url{https://code.google.com/p/lagrange/}

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

  Ree RH and Smith SA (2008). "Maximum likelihood inference
  of geographic range evolution by dispersal, local
  extinction, and cladogenesis." _Systematic Biology_,
  *57*(1), pp. 4-14. <URL:
  http://dx.doi.org/10.1080/10635150701883881>, <URL:
  http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&dopt=Citation&list_uids=18253896>.

  Foster PG (2001). "The Idiot's Guide to the Zen of
  Likelihood in a Nutshell in Seven Days for Dummies,
  Unleashed." Online PDF, widely copied, <URL:
  \url{http://www.bioinf.org/molsys/data/idiots.pdf}>.
}
\seealso{
  \code{\link[stats]{convolve}}
}

