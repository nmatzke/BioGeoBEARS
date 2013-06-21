\name{addslash}
\alias{addslash}
\title{Add a slash to a directory name if needed}
\usage{
  addslash(tmpstr)
}
\arguments{
  \item{tmpstr}{a path that you want to possibly add a
  slash to}
}
\value{
  outstr a string of the fixed path
}
\description{
  This function adds a slash to the end of the string, if
  one is not present. Handy for standardizing paths.
}
\examples{
tmpstr = "/Dropbox/_njm/__packages"
tmpstr
outstr = addslash(tmpstr)
outstr

# Annoying, getwd() often doesn't return the ending slash, which
# can make life hard for paste() later on
tmpstr = getwd()
tmpstr
outstr = addslash(tmpstr)
outstr
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\seealso{
  \code{\link[base]{getwd}}, \code{\link[base]{setwd}},
  \code{\link[base]{gsub}}
}

