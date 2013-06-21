\name{readfiles_BioGeoBEARS_run}
\alias{readfiles_BioGeoBEARS_run}
\title{Read in the stratification input files, if any}
\usage{
  readfiles_BioGeoBEARS_run(inputs)
}
\arguments{
  \item{inputs}{The inputs list}
}
\value{
  \code{inputs} The modified inputs list
}
\description{
  areas_allowed file is just a list of distance matrices,
  separated by blank lines, from youngest to oldest.
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

  Matzke_2012_IBS
}
\seealso{
  \code{\link[stats]{convolve}}
}

