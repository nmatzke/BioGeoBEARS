\name{dfnums_to_numeric}
\alias{dfnums_to_numeric}
\title{Get the class for each column in a list}
\usage{
  dfnums_to_numeric(dtf, max_NAs = 0.5, printout = FALSE,
    roundval = NULL)
}
\arguments{
  \item{dtf}{Input \code{\link[base]{data.frame}}.}

  \item{max_NAs}{Non-numeric cells will get converted to
  NA, up to the fraction of cells specified by
  \code{max_NAs}.  Above this fraction, the column is
  converted to class \code{character}.}

  \item{printout}{Print the results to screen, if desired.}

  \item{roundval}{If not NULL, \code{\link[base]{round}}
  will be run using this for the number of digits.}
}
\value{
  \code{dtf} The output \code{\link[base]{data.frame}}.
}
\description{
  This function converts each column to class
  \code{\link[base]{numeric}} where possible, and class
  \code{\link[base]{character}} otherwise.
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
  where possible.  See \code{\link{unlist_df4}} for more,
  and \code{\link{cls.df}} to see the class of each column.

  \bold{WARNING: IF A COLUMN IS A MIX OF NUMBERS AND
  NON-NUMBERS, THE NON-NUMBERS WILL BE CONVERTED TO NA IF
  THE COLUMN IS MAJORITY NUMBERS (on default; see
  \code{max_NAs}).}
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
  \code{\link{cls.df}}, \code{\link{unlist_df4}}
}

