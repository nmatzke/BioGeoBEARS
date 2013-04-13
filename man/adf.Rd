\name{adf}
\alias{adf}
\title{Convert to data.frame, without factors}
\usage{
  adf(x)
}
\arguments{
  \item{x}{matrix or other object transformable to
  data.frame}
}
\value{
  data.frame
}
\description{
  Shortcut for: \code{as.data.frame(x, row.names=NULL,
  stringsAsFactors=FALSE)}
}
\details{
  This function, and \code{\link{adf2}}, are useful for
  dealing with errors due to automatic conversion of some
  columns to factors.  Another solution may be to prepend
  \code{options(stringsAsFactors = FALSE)} at the start of
  one's script, to turn off all default stringsAsFactors
  silliness.
}
\examples{
x = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
adf(x)
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\seealso{
  \code{\link{adf2}}
}

