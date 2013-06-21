\name{simulated_indexes_to_tipranges_file}
\alias{simulated_indexes_to_tipranges_file}
\title{Convert simulated Qmat 0-based indexes to a tipranges file}
\usage{
  simulated_indexes_to_tipranges_file(simulated_states_by_node,
    areas_list, states_list, trfn,
    out_geogfn = "lagrange_area_data_file.data")
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

  \item{out_geogfn}{The output filename.}
}
\value{
  \code{out_geogfn} The output filename.
}
\description{
  This function takes simulated state indices (ranging from
  0 to numstates-1, i.e. number of possible geographic
  ranges-1) and converts them to a
  C++-\code{LAGRANGE}-style PHYLIP geographic ranges file.
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
  \code{\link{simulated_indexes_to_tipranges_object}}
}

