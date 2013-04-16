\name{slashslash}
\alias{slashslash}
\title{Remove double slash (slash a slash)}
\usage{
  slashslash(tmpstr)
}
\arguments{
  \item{tmpstr}{a path that you want to remove double
  slashes from}
}
\value{
  outstr a string of the fixed path
}
\description{
  Shortcut for: \code{gsub(pattern="//", replacement="/",
  x=tmpstr)}
}
\details{
  This function is useful for removing double slashes that
  can appear in full pathnames due to inconsistencies in
  trailing slashes in working directories etc.
}
\examples{
tmpstr = "/Library/Frameworks/R.framework/Versions/2.15/Resources/library/modiscloud/extdata/2002_raw_MYD//MYD03.A2002185.0645.005.2009192031332.hdf"

outstr = slashslash(tmpstr)
outstr
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\seealso{
  \code{\link[base]{getwd}}, \code{\link[base]{setwd}},
  \code{\link[base]{gsub}}
}

