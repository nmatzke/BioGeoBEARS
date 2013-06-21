\name{mix_colors_for_states}
\alias{mix_colors_for_states}
\title{Mix colors logically to produce colors for multi-area ranges}
\usage{
  mix_colors_for_states(colors_matrix,
    states_list_0based_index, exclude_null = TRUE)
}
\arguments{
  \item{colors_matrix}{A column with a color for each
  single area}

  \item{states_list_0based_index}{States list giving areas,
  0-based}

  \item{exclude_null}{If TRUE, null ranges are excluded
  (however coded). Default TRUE.}
}
\value{
  \code{colors_list_for_states} The colors for the ML
  states
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

