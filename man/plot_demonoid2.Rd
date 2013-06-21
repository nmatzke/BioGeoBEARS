\name{plot_demonoid2}
\alias{plot_demonoid2}
\title{Plot LaplacesDemon MCMC output}
\usage{
  plot_demonoid2(x, BurnIn = 0, Data = NULL, PDF = FALSE,
    Parms = NULL, ...)
}
\arguments{
  \item{x}{The MCMC output}

  \item{BurnIn}{The burnin amount}

  \item{Data}{The
  \code{\link[LaplacesDemon]{LaplacesDemon}} input data.}

  \item{PDF}{If TRUE, plot to PDF.}

  \item{Parms}{The parameters.}

  \item{...}{Additional arguments to standard functions}
}
\value{
  Nothing.
}
\description{
  Modified version of
  \code{\link[LaplacesDemon]{plot.demonoid}} to plot to
  PDF.
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

