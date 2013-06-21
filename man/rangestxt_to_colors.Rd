\name{rangestxt_to_colors}
\alias{rangestxt_to_colors}
\title{Convert a list of ranges text (KOM, MH, KOMIH, etc.)}
\usage{
  rangestxt_to_colors(possible_ranges_list_txt,
    colors_list_for_states, MLstates)
}
\arguments{
  \item{possible_ranges_list_txt}{A list of the allowed
  ranges/states}

  \item{colors_list_for_states}{The corresponding colors}

  \item{MLstates}{The ML states for the internal nodes}
}
\value{
  \code{MLcolors} The colors for the ML states
}
\description{
  Like it says.
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
  \code{\link[stats]{optim}}
}

