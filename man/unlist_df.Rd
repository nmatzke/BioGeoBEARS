\name{unlist_df}
\alias{unlist_df}
\title{Unlist the columns in a data.frame}
\usage{
  unlist_df(df)
}
\arguments{
  \item{df}{matrix or other object transformable to
  data.frame}
}
\value{
  data.frame
}
\description{
  Sometimes, matrices or data.frames will malfunction due
  to their having lists as columns and other weirdness.
  This is a shortcut for \code{data.frame(lapply(df,
  function(x) unlist(x)))}.
}
\examples{
x = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
unlist_df2(x)
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\seealso{
  \code{\link{unlist_df2}}
}

