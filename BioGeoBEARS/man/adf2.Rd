\name{adf2}
\alias{adf2}
\title{Convert to data.frame, without factors}
\usage{
  adf2(x)
}
\arguments{
  \item{x}{matrix or other object transformable to
  data.frame}
}
\value{
  data.frame
}
\description{
  Shortcut for: \code{tmp_rownames = 1:nrow(x);
  as.data.frame(x, row.names=tmp_rownames,
  stringsAsFactors=FALSE)}
}
\details{
  This function, and \code{\link{adf2}}, are useful for
  dealing with errors due to automatic conversion of some
  columns to factors.  Another solution may be to prepend
  \code{options(stringsAsFactors = FALSE)} at the start of
  one's script, to turn off all default stringsAsFactors
  silliness.

  In adf2, rownames are forced to be numbers; this can
  prevent errors due to e.g. repeated rownames after an
  \code{rbind} operation.
}
\examples{
x = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
adf2(x)
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\seealso{
  \code{\link{adf}}
}

