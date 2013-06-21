\name{postorder_nodes_phylo4_return_table}
\alias{postorder_nodes_phylo4_return_table}
\title{Get a table of node numbers, including DIVA node numbers}
\usage{
  postorder_nodes_phylo4_return_table(tr4)
}
\arguments{
  \item{tr4}{A tree object in \code{\link[ape]{phylo}} or
  \code{\link[phylobase]{phylo4}} format.}
}
\value{
  \code{postorder_table} A data.frame showing the various
  corresponding node numbers.
}
\description{
  Various programs (annoyingly) label internal nodes in
  different ways.  This function shows the corresponding
  node numbers for several different systems. This table
  can then be used to translate, when the user wishes to
  plot the output from various programs on the nodes of a
  tree.  In particular, the last column contains the DIVA
  node-numbering scheme (\cite{Ronquist (1996)},
  \cite{Ronquist (1997)}).
}
\details{
  There are many ways of numbering nodes in a tree.  This
  returns a matrix containing (column 1) R's native
  internal numbering scheme, and (column 2) the node
  numbers in the downpass numbering used by C++
  \code{LAGRANGE}, in particular in their .bgkey output
  file.  Note that this is different from
  \code{\link[ape]{ape}}'s \code{pruningwise} downpass
  ordering (see \code{\link{get_pruningwise_nodenums}}).

  The python version of LAGRANGE labels internal nodes
  differently (sigh), but they are in the same order at
  least, so can just be renumbered from 1 to
  \code{tr$Nnode} to get them to match the C++
  \code{LAGRANGE} node numbering.

  DIVA has yet a different node numbering scheme; see
  \code{\link{postorder_nodes_phylo4_return_table}}
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

  Ronquist1996_DIVA

  Ronquist_1997_DIVA

  Matzke_2012_IBS

  ReeSmith2008
}
\seealso{
  \code{\link{get_pruningwise_nodenums}},
  \code{\link{get_lagrange_nodenums}}, \code{\link{prt}}
}

