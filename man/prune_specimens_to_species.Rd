\name{prune_specimens_to_species}
\alias{prune_specimens_to_species}
\title{Take a tree and species names/geography table and produce a pruned tree and tipranges object}
\usage{
  prune_specimens_to_species(original_tr, xls,
    group_name = "default", titletxt = "",
    areas_abbr = NULL, plot_intermediate = TRUE)
}
\arguments{
  \item{original_tr}{The input tree (an
  \code{\link[ape]{ape}} \code{\link[ape]{phylo}} object).}

  \item{xls}{The input table (a
  \code{\link[base]{data.frame}})}

  \item{group_name}{The name of the clade in the tree.  For
  use in plots and output files. Default="default".}

  \item{titletxt}{Additional text for the plots. Default
  "".}

  \item{areas_abbr}{An optional table, containing the
  abbreviations (e.g. letters) corresponding to each region
  in xls$region. Default is NULL, in which case the program
  imposes A, B, C, D, etc. \code{areas_abbr} must have
  column headings \code{abbr} and \code{letter}.}

  \item{plot_intermediate}{If TRUE, the starting, ending,
  and intermediate stages of tree pruning are plotted.}
}
\value{
  The outputs are a \code{\link[base]{list}} with a pruned
  tree and, optionally, a tipranges object.
}
\description{
  This function takes a tree and species names/geography
  table and produces a pruned tree and (optionally) a
  tipranges object.
}
\details{
  Often, users will have an phylogeny where the tips/OTUs
  (operational taxonomic units) are specimens rather than
  species.  The analyses done by models like DEC, DEC+J,
  etc., in programs like LAGRANGE and BioGeoBEARS, assume
  as a core part of the model that species might occupy
  more than one areas.  A phylogeny of specimens, then,
  would not be an appropriate input to these programs, as
  each single specimen can only be found in one region. The
  exception would occur when the researcher is confident
  that each species lives in only one region; in that case,
  the specimen geography is representative of the species
  geography.

  This function requires a table containing

  (1) Column "OTUs": all tipnames in the input tree (often,
  original specimen/original OTU names) );

  (2) Column "species": the corresponding species names;

  (3) optionally, the geographic range inhabited by each
  specimen (column "region").  If an OTU has more than one
  geographic range in the original table, these should be
  split by "|".

  When the pruning occurs, all tips belonging to the same
  species are cut, except the first.

  NOTE: Tips that should be cut because they are outgroups,
  or because they are geographically outside of your domain
  of analysis, should be represented in xls$region by
  "out_group" or "Out". These will be cut from the final
  tree/geography table.
}
\note{
  Go BEARS!
}
\examples{
testval=1
tipranges_object = define_tipranges_object()
tipranges_object

areanames = getareas_from_tipranges_object(tipranges_object)
areanames
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\seealso{
  \code{\link[ape]{drop.tip}},
  \code{\link{define_tipranges_object}},
}

