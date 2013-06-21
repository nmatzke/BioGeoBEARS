\name{get_simstates}
\alias{get_simstates}
\title{Load the simulation information from an underscore delimited text string.}
\usage{
  get_simstates(simhist_row)
}
\arguments{
  \item{simhist_row}{A row from a table, which must have a
  column named \code{simulated_states_by_node_txt}.}
}
\value{
  \code{simulated_states_by_node} A numeric vector of
  0-based state indices.
}
\description{
  If the simulated states are stored in a big text file, it
  can be useful to store them as a single string in a
  single cell per row, so that the number of columns
  doesn't have to change with each different-sized tree.
  This function extracts the simulated states from this
  format.
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

  Matzke_2012_IBS
}
\seealso{
  \code{\link[utils]{read.table}}
}

