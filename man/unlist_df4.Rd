\name{unlist_df4}
\alias{unlist_df4}
\title{Unlist the columns in a data.frame, with more checks, adf, and dfnums_to_numeric}
\usage{
  unlist_df4(df, ...)
}
\arguments{
  \item{df}{matrix or other object transformable to
  data.frame}

  \item{...}{Additional options passed to
  \code{\link{dfnums_to_numeric}}.}
}
\value{
  \code{outdf} data.frame
}
\description{
  Sometimes, matrices or data.frames will malfunction due
  to their having lists as columns and other weirdness.
  This runs \code{\link[base]{unlist}} and additional
  checks, and forces conversion to a
  \code{\link[base]{data.frame}} at the end.  It also adds
  \code{\link{dfnums_to_numeric}} which should remove the
  problem of numbers columns being of class
  \code{\link[base]{character}}.
}
\details{
  See especially \code{\link[base]{data.matrix}} for a
  possibly simpler alternative.
}
\examples{
x = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
cls.df(x)
unlist_df4(x)

x = matrix(c(1,2,3,4,5,"A"), nrow=3, ncol=2)
cls.df(x)
unlist_df4(x)

x = adf(matrix(c(1,2,3,4,5,"A"), nrow=3, ncol=2))
names(x) = c("A","B")
cls.df(x)
unlist_df4(x)
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\seealso{
  \code{\link{unlist_df}}, \code{\link{dfnums_to_numeric}},
  \code{\link{cls.df}}, \code{\link[base]{data.matrix}}
}

