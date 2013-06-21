\name{LGcpp_MLstate_per_node}
\alias{LGcpp_MLstate_per_node}
\title{Get the ML states per node, from a states table}
\usage{
  LGcpp_MLstate_per_node(states)
}
\arguments{
  \item{states}{A data.frame containing the node numbers,
  states, and state probabilities.}
}
\value{
  \code{MLstates} A data.frame containing the node numbers,
  ML states, and state probabilities.
}
\description{
  Given a table of states probabilities from either
  \code{\link{LGcpp_states_fn_to_table}} or
  \code{\link{LGcpp_states_fn_to_table}}, get the ML state
  for each node.
}
\details{
  See \code{\link{get_lagrange_nodenums}} for connecting
  these node numbers to APE node numbers.
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
  \code{\link{LGcpp_states_fn_to_table}},
  \code{\link{LGcpp_states_fn_to_table}}
}

