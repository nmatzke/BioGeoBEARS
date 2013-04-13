\name{nullsym_to_NA}
\alias{nullsym_to_NA}
\title{Convert a specified null range code to NA}
\usage{
  nullsym_to_NA(mat, nullsym = "-")
}
\arguments{
  \item{mat}{A matrix.}

  \item{nullsym}{A character specifying the null symbol.}
}
\value{
  \code{mat} The revised matrix
}
\description{
  Takes a matrix \code{mat}, converts any instances of the
  \code{nullsym} symbol to \code{NA}.
}
\note{
  Go BEARS!
}
\examples{
testval=1
mat = matrix(c("-",1,1,1,"-",1,1,1,"-"), nrow=3)
mat
mat2 = nullsym_to_NA(mat, nullsym="-")
mat2
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
  \code{\link{remove_null_rowcols_from_mat}}
}

