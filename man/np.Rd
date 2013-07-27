\name{np}
\alias{np}
\title{normalizePath shortcut}
\usage{
  np(path = path, ...)
}
\arguments{
  \item{path}{The path to run
  \code{\link[base]{normalizePath}} on.}

  \item{...}{Additional arguments to
  \code{\link[base]{normalizePath}}.}
}
\value{
  \code{path} The path that was normalized.
}
\description{
  Utility function that runs
  \code{\link[base]{normalizePath}}. Useful for running on
  Mac vs. Windows.
}
\note{
  Go BEARS!
}
\examples{
# Get a path
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
extdata_dir

path = paste(extdata_dir, "//", "Psychotria_5.2.newick", sep="")
path

path = np(path)
path
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link[base]{normalizePath}}
}

