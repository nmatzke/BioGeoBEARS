\name{get_path_first}
\alias{get_path_first}
\title{Get the text that comes before the last slash}
\usage{
  get_path_first(inpath, addslash = "FALSE")
}
\arguments{
  \item{inpath}{A string of class
  \code{\link[base]{character}}.}

  \item{addslash}{If \code{TRUE}, add a slash at the end of
  the path.}
}
\value{
  \code{outpath} A string with the full path, without the
  file.
}
\description{
  Extracts the path from a full path, removing the
  filename.
}
\note{
  Go BEARS!
}
\examples{
get_path_first("/Users/nickm/Library/R/2.15/library/BioGeoBEARS/extdata/Psychotria_geog.data")
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

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
  \code{\link{get_path_last}}
}

