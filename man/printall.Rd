\name{printall}
\alias{printall}
\title{Print an entire table to screen}
\usage{
  printall(dtf, chunksize_toprint = 40, printflag = TRUE)
}
\arguments{
  \item{dtf}{The \code{\link[base]{data.frame}} to
  \code{\link[base]{print}}.}

  \item{chunksize_toprint}{Number of lines to print.
  Default 50.}

  \item{printflag}{For optional printing. Passed to
  \code{\link{prflag}}.}
}
\value{
  NULL
}
\description{
  Utility function.  This prints a table to screen in
  chunks of \code{chunksize_toprint} (default=40).  This
  avoids the annoying situation of not being able to see
  the bottom of a table. Note that if you print something
  huge, you will be waiting for awhile (try ESC or CTRL-C
  to cancel such an operation).
}
\details{
  Another option is to reset options to something like:
  \code{options(max.print=99999)}, but this is hard to
  remember.  Your current setting is
  \code{getOption("max.print")}.
}
\note{
  Go BEARS!
}
\examples{
test=1
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link[base]{print}}, \code{\link{prflag}}
}

