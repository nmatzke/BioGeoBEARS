\name{section_the_tree}
\alias{section_the_tree}
\title{Section a tree for stratified analysis}
\usage{
  section_the_tree(inputs, make_master_table = FALSE,
    plot_pieces = TRUE, cut_fossils = TRUE,
    fossils_older_than = 0.6)
}
\arguments{
  \item{inputs}{The list of inputs for stratified analysis}

  \item{make_master_table}{If desired, make an
  \code{inputs$master_table} containing the correspondance
  between the original tree and the sectioned pieces.}

  \item{plot_pieces}{If \code{TRUE}, plot the tree chunks
  (but not isolated branch segments) as they are created.}

  \item{cut_fossils}{If \code{TRUE} (default), the program
  is stopped if there are fossils, i.e. tips older than 0.6
  my (default).  Users should use
  code{\link[ape]{drop.tip}} or an external program to clip
  fossils out of the tree. PLEASE NOTE that several times I
  have experienced miserable long nights due, apparently,
  to \code{\link[ape]{drop.tip}} producing weird tree
  structures, resulting in weird Newick files, without me
  realizing it.  The solution is usually to open the Newick
  file in something like \code{FigTree}, resort the
  branches, and save to a new Newick file. Fossils have now
  been implemented in stratified analysis; this was
  complicated, as it involves inserting new branches in
  chopped trees.}

  \item{fossils_older_than}{Tips that are older than
  \code{fossils_older_than} will be marked as \code{TRUE}
  in a column called \code{fossil}. This is not currently
  set to 0, because Newick files can have slight precision
  issues etc. that mean not all tips quite come to zero.
  You can attempt to fix this with
  \code{\link{extend_tips_to_ultrametricize}} (but make
  sure you do not inappropriately average in fossils!!).}
}
\value{
  \code{inputs} with \code{inputs$tree_sections_list}
  added.
}
\description{
  A utility function for stratified analysis.  Sections the
  tree into a series of strata. Each stratum may have one
  or more subtrees (APE phylo3 objects, *WITH* root edges)
  and/or branch segments (which are just represented as
  numeric values, indicating the length of the sub-branch,
  i.e. the time-width of the stratum, if the branch crosses
  the whole stratum.
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

  Matzke_2012_IBS
}
\seealso{
  \code{\link{prt}}, \code{\link{chainsaw2}},
  \code{\link[ape]{drop.tip}}
}

