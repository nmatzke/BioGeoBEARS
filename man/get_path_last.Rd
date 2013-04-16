\name{get_path_last}
\alias{get_path_last}
\title{Get the text that comes after the last slash}
\usage{
  get_path_last(path)
}
\arguments{
  \item{path}{A string of class
  \code{\link[base]{character}}.}
}
\value{
  \code{lastword} A string with the filename, without the
  path.
}
\description{
  Extracts the filename from a full path.
}
\note{
  Go BEARS!
}
\examples{
get_path_last("/Users/nickm/Library/R/2.15/library/BioGeoBEARS/extdata/Psychotria_geog.data")
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
  \code{\link{get_path_first}}
}

