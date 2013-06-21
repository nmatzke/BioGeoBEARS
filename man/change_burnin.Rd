\name{change_burnin}
\alias{change_burnin}
\title{Change the burnin fraction when re-parameterizing LaplacesDemon MCMC output}
\usage{
  change_burnin(Fit, MyData, burnfract = 0.1,
    newsamps = NULL)
}
\arguments{
  \item{Fit}{The \code{\link[LaplacesDemon]{LaplacesDemon}}
  output object.}

  \item{MyData}{The
  \code{\link[LaplacesDemon]{LaplacesDemon}} input data.}

  \item{burnfract}{The new burnin fraction}

  \item{newsamps}{If you want to subset the samples.}
}
\value{
  \code{Fit} The modified MCMC output.
}
\description{
  \code{\link[LaplacesDemon]{LaplacesDemon}} likes to run
  its MCMC sampling on a simple number line.  Thus, the
  likelihood function etc. should transform the numbers
  into the range desired, e.g. 0-1.
}
\details{
  This function changes the burnin fraction
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

  LaplacesDemon_Tutorial

  Matzke_2012_IBS
}
\seealso{
  \code{\link[LaplacesDemon]{LaplacesDemon}}
}

