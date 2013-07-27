\name{get_lagrange_nodenums}
\alias{get_lagrange_nodenums}
\title{Get internal node numbers in LAGRANGE's downpass order}
\usage{
  get_lagrange_nodenums(tr)
}
\arguments{
  \item{tr}{A \code{\link[ape]{phylo}} tree object}
}
\value{
  \code{downpass_node_matrix} A matrix of node numbers
}
\description{
  There are many ways of numbering nodes in a tree.  This
  returns a matrix containing (column 1) R's native
  internal numbering scheme, and (column 2) the node
  numbers in the downpass numbering used by C++
  \code{LAGRANGE}, in particular in their .bgkey output
  file.  Note that this is different from
  \code{\link[ape]{ape}}'s \code{pruningwise} downpass
  ordering (see \code{\link{get_pruningwise_nodenums}}).
}
\details{
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
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
tmppath = paste(extdata_dir,
"/examples/Psychotria_M0/LGcpp/Psychotria_5.2.newick", sep="")
trfn = np(slashslash(tmppath))
tr = read.tree(trfn)
downpass_node_matrix = get_lagrange_nodenums(tr)
downpass_node_matrix


downpass_node_matrix = get_lagrange_nodenums(tr)
downpass_node_matrix = downpass_node_matrix[order(downpass_node_matrix[,2]), ]
plot(tr)
nodelabels(node=20:37, downpass_node_matrix[,1])
tiplabels(1:19)

plot(tr)
nodelabels(node=20:37, downpass_node_matrix[,2])
tiplabels(1:19)

downpass_node_matrix = get_lagrange_nodenums(tr)
downpass_node_matrix = downpass_node_matrix[order(downpass_node_matrix[,1]), ]
plot(tr)
nodelabels(node=20:37, downpass_node_matrix[,1])
tiplabels(1:19)

# THIS WORKS
plot(tr)
nodelabels(node=20:37, downpass_node_matrix[,2])
tiplabels(1:19)
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link{get_pruningwise_nodenums}},
  \code{\link{prt}},
  \code{\link{postorder_nodes_phylo4_return_table}}
}

