\name{cls.df}
\alias{cls.df}
\title{Get the class for each column in a list}
\usage{
  cls.df(dtf, printout = FALSE)
}
\arguments{
  \item{dtf}{Input \code{\link[base]{data.frame}}.}

  \item{printout}{Print the results to screen, if desired.}
}
\value{
  \code{dtf_classes} A \code{\link[base]{data.frame}}
  showing the column, column name, and column class.
}
\description{
  This function returns the \code{\link[base]{class}} of
  each column in a \code{\link[base]{data.frame}}.
}
\details{
  R does lots of weird and unpredictable things when you
  build up tables/matrices/data.frames by e.g.
  \code{\link[base]{cbind}} and \code{\link[base]{rbind}}
  on vectors of results.  The major problems are (1)
  columns get made into class \code{\link[base]{list}}; (2)
  \code{\link[base]{numeric}} columns are converted to
  class \code{\link[base]{factor}}; (3)
  \code{\link[base]{numeric}} columns are converted to
  class \code{\link[base]{character}}; (4) you have a
  \code{\link[base]{matrix}} when you think you have a
  \code{\link[base]{data.frame}}.

  All of this could be taken care of by detailed
  understanding and tracking of when R recasts values in
  vectors, matrices, and data frames...but this is a huge
  pain, it is easier to just have a function that jams
  everything back to a \code{\link[base]{data.frame}} with
  no lists, no factors, and with columns being numeric
  where possible.  See \code{\link{dfnums_to_numeric}} and
  \code{\link{unlist_df4}} for these options.
}
\note{
  Go BEARS!
}
\examples{
x = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
cls.df(x)
dfnums_to_numeric(adf(x))
unlist_df4(x)

x = matrix(c(1,2,3,4,5,"A"), nrow=3, ncol=2)
cls.df(x)
dfnums_to_numeric(adf(x))
unlist_df4(x)

x = adf(matrix(c(1,2,3,4,5,"A"), nrow=3, ncol=2))
names(x) = c("A","B")
cls.df(x)
dfnums_to_numeric(adf(x))
unlist_df4(x)
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link{dfnums_to_numeric}}, \code{\link{unlist_df4}}
}

