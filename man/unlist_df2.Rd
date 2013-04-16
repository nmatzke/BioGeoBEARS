\name{unlist_df2}
\alias{unlist_df2}
\title{Unlist the columns in a data.frame, with more checks}
\usage{
  unlist_df2(df)
}
\arguments{
  \item{df}{matrix or other object transformable to
  data.frame}
}
\value{
  \code{outdf} A \code{\link[base]{matrix}}.
}
\description{
  Sometimes, matrices or data.frames will malfunction due
  to their having lists as columns and other weirdness.
  This runs \code{\link{unlist}} and additional checks.
}
\examples{
x = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
unlist_df2(x)
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\seealso{
  \code{\link{unlist_df}}
}

