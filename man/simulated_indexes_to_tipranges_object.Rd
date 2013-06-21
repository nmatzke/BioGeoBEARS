\name{simulated_indexes_to_tipranges_object}
\alias{simulated_indexes_to_tipranges_object}
\title{Convert simulated Qmat 0-based indexes to a tipranges object}
\usage{
  simulated_indexes_to_tipranges_object(simulated_states_by_node,
    areas_list, states_list, trfn)
}
\arguments{
  \item{simulated_states_by_node}{The simulated
  states/geographic ranges, in 0-based index form, ordered
  as the tips & nodes are ordered in a
  \code{pruningwise}-ordered \code{phylo} object in
  \code{APE}.}

  \item{areas_list}{A list of the desired area
  names/abbreviations/letters.}

  \item{states_list}{A list of the possible
  states/geographic ranges, in 0-based index form.}

  \item{trfn}{The filename of the source Newick tree.}
}
\value{
  \code{tipranges_object} An object of class
  \code{tipranges}.
}
\description{
  This function takes simulated state indices (ranging from
  0 to numstates-1, i.e. number of possible geographic
  ranges-1) and converts them to a tipranges object.  This
  can then be converted into a C++-\code{LAGRANGE}-style
  PHYLIP geographic ranges file.
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

  Matzke_2012_IBS

  ReeSmith2008

  SmithRee2010_CPPversion
}
\seealso{
  \code{\link{define_tipranges_object}},
  \code{\link{getareas_from_tipranges_object}},
  \code{\link{simulated_indexes_to_tipranges_file}}
}

