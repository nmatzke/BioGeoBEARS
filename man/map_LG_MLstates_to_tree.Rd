\name{map_LG_MLstates_to_tree}
\alias{map_LG_MLstates_to_tree}
\title{Map states to the nodes on a phylogeny}
\usage{
  map_LG_MLstates_to_tree(MLstates_LGcpp, tr, tipranges,
    removechar = NULL, type = "C++",
    statesColors_table = "default", bgcol = "green3",
    areanames = "default", newplot = TRUE, ...)
}
\arguments{
  \item{MLstates_LGcpp}{A data.frame containing the node
  numbers, states, and states probabilities.}

  \item{tr}{An ape phylo object}

  \item{tipranges}{Tipranges object}

  \item{removechar}{The character to remove, if needed.}

  \item{type}{The type of LAGRANGE input (default C++)}

  \item{statesColors_table}{If not default, a table with a
  color for each area combination.}

  \item{bgcol}{The background color}

  \item{areanames}{The area names, if different from those
  in the tipranges object}

  \item{newplot}{Default TRUE; should there be a new plot,
  or should the splits be added to another plot?}

  \item{...}{Additional arguments to standard functions}
}
\value{
  \code{MLstates_LGcpp} The states table, ordered
  appropriately.
}
\description{
  What it says.
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
  \url{https://code.google.com/p/lagrange/}

  Matzke_2012_IBS

  ReeSmith2008
}
\seealso{
  \code{\link{get_lagrange_nodenums}},
  \code{\link{LGpy_splits_fn_to_table}},
  \code{\link{LGcpp_splits_fn_to_table}}
}

