\name{unlist_df3}
\alias{unlist_df3}
\title{Unlist the columns in a data.frame, with more checks and adf}
\usage{
  unlist_df3(df)
}
\arguments{
  \item{df}{matrix or other object transformable to
  data.frame}
}
\value{
  \code{outdf} data.frame
}
\description{
  Sometimes, matrices or data.frames will malfunction due
  to their having lists as columns and other weirdness.
  This runs \code{\link[base]{unlist}} and additional
  checks, and forces conversion to a
  \code{\link[base]{data.frame}} at the end.
}
\examples{
x = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
unlist_df3(x)
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\seealso{
  \code{\link{unlist_df}}
}

