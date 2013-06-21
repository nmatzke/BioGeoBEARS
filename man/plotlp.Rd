\name{plotlp}
\alias{plotlp}
\title{Plot LaplacesDemon plots after excluding some burnin}
\usage{
  plotlp(Fit, MyData, PDF = TRUE, Parms = NULL,
    burnfract = 0, ...)
}
\arguments{
  \item{Fit}{The \code{\link[LaplacesDemon]{LaplacesDemon}}
  output object.}

  \item{MyData}{The
  \code{\link[LaplacesDemon]{LaplacesDemon}} input data.}

  \item{PDF}{Plot to a PDF? (Recommended!)}

  \item{Parms}{The parameters}

  \item{burnfract}{The burnin fraction}

  \item{...}{Additional arguments to standard functions}
}
\value{
  The plot
}
\description{
  Takes the objects output by a
  \code{\link[LaplacesDemon]{LaplacesDemon}} MCMC search,
  and plots LnP, parameter values, etc.
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
  \code{\link[LaplacesDemon]{LaplacesDemon}}
}

