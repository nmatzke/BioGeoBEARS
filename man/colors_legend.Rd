\name{colors_legend}
\alias{colors_legend}
\title{Plot a colors legend for geographic ranges}
\usage{
  colors_legend(possible_ranges_list_txt,
    colors_list_for_states, legend_ncol = NULL,
    legend_cex = 1)
}
\arguments{
  \item{possible_ranges_list_txt}{A list of the allowed
  ranges/states}

  \item{colors_list_for_states}{The corresponding colors}

  \item{legend_ncol}{The number of columns in the legend.
  If \code{NULL} (default), the function calculates
  \code{floor(sqrt(length(possible_ranges_list_txt) / 2))}.
  Note that when you have hundreds of states, there is
  probably no good way to have a coherent legend, and it is
  easier to just rely upon printing the character codes for
  the ML states in the plots, with the colors, and users
  can then see and trace the common colors/states by eye.}

  \item{legend_cex}{The cex (character expansion size) for
  the legend.  Defaults to 1, which means the
  \code{\link[graphics]{legend}} function determines the
  size.  The value 2.5 works well for 15 or 16
  states/ranges.}
}
\value{
  Nothing
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
  \code{\link[graphics]{legend}},
  \code{\link[base]{floor}}, \code{\link[base]{ceiling}},
  \code{\link[base]{floor}}
}

