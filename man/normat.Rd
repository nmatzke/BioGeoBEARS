\name{normat}
\alias{normat}
\title{Utility functions to help deal with matrices}
\usage{
  normat(relative_matrix)
}
\arguments{
  \item{relative_matrix}{A square matrix giving the
  relative probabilities/weights of transitions.}
}
\value{
  \code{m} A Q matrix, i.e. normalized transition matrix
  (Qmat)
}
\description{
  Normalize a transition matrix
}
\details{
  \code{normat} normalizes a square transition matrix, such
  that each row sums to 0, and the diagonal equals the
  negative of the sum of the rest of the cells in the row.
  This matrix can then be exponentiated by values of
  \code{t} (time or another measure of branch length) to
  produce transition probabilities for any given value of
  \code{t}.

  See \cite{Foster (2001)} for a succinct summary of
  transition matrices and their exponentiation.
}
\note{
  Go BEARS!
}
\examples{
testval=1
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
  \url{https://code.google.com/p/lagrange/}

  Foster PG (2001). "The Idiot's Guide to the Zen of
  Likelihood in a Nutshell in Seven Days for Dummies,
  Unleashed." Online PDF, widely copied, <URL:
  \url{http://www.bioinf.org/molsys/data/idiots.pdf}>.
}

