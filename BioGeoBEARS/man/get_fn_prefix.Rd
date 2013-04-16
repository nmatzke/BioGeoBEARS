\name{get_fn_prefix}
\alias{get_fn_prefix}
\title{Get everything BEFORE the last suffix (.nex or whatever)}
\usage{
  get_fn_prefix(fn)
}
\arguments{
  \item{fn}{The input filename.}
}
\value{
  \code{prefix} The output string.
}
\description{
  Extracts the string from before the last suffix.  I.e.,
  "filename.nex" becomes "filename".
}
\note{
  Go BEARS!
}
\examples{
get_fn_prefix("/Users/nickm/Library/R/2.15/library/BioGeoBEARS/extdata/Psychotria_geog.data")
get_fn_prefix("Psychotria_geog.data")
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
  \code{\link{get_path_last}}, \code{\link{get_path_first}}
}

