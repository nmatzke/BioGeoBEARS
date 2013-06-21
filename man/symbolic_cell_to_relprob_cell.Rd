\name{symbolic_cell_to_relprob_cell}
\alias{symbolic_cell_to_relprob_cell}
\title{Convert symbolic cell (a text equation) to relprob matrix (a numeric value).}
\usage{
  symbolic_cell_to_relprob_cell(charcell, cellsplit = "",
    mergesym = "*", d = 0.1, e = 0.01, ...)
}
\arguments{
  \item{charcell}{The text formula.}

  \item{cellsplit}{The symbol to split the formulas on.
  Default "\\\\+" (plus symbol, with escape code).}

  \item{mergesym}{The symbol to merge the formulas with.
  Default "+".}

  \item{d}{The dispersal/range expansion rate. Default
  \code{d=0.1}.}

  \item{e}{The extinction/range contraction rate. Default
  \code{e=0.01}.}

  \item{...}{Additional arguments to pass to
  \code{\link[base]{strsplit}}.}
}
\value{
  \code{cellval} The output cell value.
}
\description{
  This is a utility function for
  \code{\link{symbolic_to_P_matrix}} and
  \code{\link{symbolic_to_Q_matrix}}.
}
\details{
  This function can be used in \code{\link[base]{sapply}}.
  It still will not be very fast compared to the
  calculations in \code{\link[cladoRcpp]{cladoRcpp}}, but
  can be useful for demonstrative purposes.
}
\note{
  Go BEARS!
}
\examples{
testval=1

charcell = "1*d+1*d"

# Right
cellval = symbolic_cell_to_relprob_cell(charcell, cellsplit="yadda",
mergesym="", d=0.1, e=0.01)
cellval

# Wrong
cellval = symbolic_cell_to_relprob_cell(charcell, cellsplit="\\\\+",
mergesym="*", d=0.1, e=0.01)
cellval

# Right
cellval = symbolic_cell_to_relprob_cell(charcell, cellsplit="\\\\+",
mergesym="+", d=0.1, e=0.01)
cellval
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link{symbolic_to_P_matrix}}
}

