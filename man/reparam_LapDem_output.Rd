\name{reparam_LapDem_output}
\alias{reparam_LapDem_output}
\title{Re-parameterize LaplacesDemon MCMC output}
\usage{
  reparam_LapDem_output(Fit, MyData,
    transfun = transform_with_logistic)
}
\arguments{
  \item{Fit}{The \code{\link[LaplacesDemon]{LaplacesDemon}}
  output object.}

  \item{MyData}{The
  \code{\link[LaplacesDemon]{LaplacesDemon}} input data.}

  \item{transfun}{The function to use for the
  transformation.}
}
\value{
  \code{Fit} The transformed MCMC output.
}
\description{
  \code{\link[LaplacesDemon]{LaplacesDemon}} likes to run
  its MCMC sampling on a simple number line.  Thus, the
  likelihood function etc. should transform the numbers
  into the range desired, e.g. 0-1.
}
\details{
  This function transforms the
  \code{\link[LaplacesDemon]{LaplacesDemon}} output
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

