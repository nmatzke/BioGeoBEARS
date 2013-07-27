\name{corner_coords}
\alias{corner_coords}
\title{Get the corner coordinates}
\usage{
  corner_coords(tr, coords_fun = "plot_phylo3_nodecoords",
    tmplocation = "manual")
}
\arguments{
  \item{tr}{A tree object in \code{\link[ape]{phylo}}
  format.}

  \item{coords_fun}{The name of the function to use to get
  node coordinates. Default: "plot_phylo3_nodecoords".}

  \item{tmplocation}{Default is "manual", which throws an
  error check unless your path structure matches the
  developer's.  Most users should probably use the
  \code{\link[base]{system.file}} command in the examples,
  below. The directory location containing the R script
  \code{plot_phylo3_nodecoords.R}. This function, modified
  from the \code{\link[ape]{ape}} function
  \code{\link[ape]{plot.phylo}}, cannot be included
  directly in the R package as it contains C code that does
  not pass CRAN's R CMD check. The default,
  cornercoords_loc="manual", will not allow split states to
  be plot.  The R script \code{plot_phylo3_nodecoords.R} is
  located in the BioGeoBEARS extension data directory,
  \code{extdata/a_scripts}.  You should be able to get the
  full path with
  \code{list.files(system.file("extdata/a_scripts",
  package="BioGeoBEARS"), full.names=TRUE)}.}
}
\value{
  \code{corners_list}
}
\description{
  Gets the coordinates of the corners when the tree is
  plotted.
}
\details{
  Because this function needs to use a modified version of
  the APE plot.phylo function, and for complex reasons
  APE's .C functions cannot be used elsewhere without
  causing problems with R CMD check, this function is left
  up to user specification.  Basically, the user puts in
  the name of the function, which is available in the
  extension data (\code{extdata/a_scripts}) directory of
  the package.  The defaults work on the developer's
  machine, other users may have to e.g. change "manual" to
  \code{tmplocation}, where \code{tmplocation} is specified
  as in the example.
}
\note{
  Go BEARS!
}
\examples{
# Set location like this if you don't have plot_phylo3_nodecoords
# hardcoded/sourced elsehwhere
# tmplocation = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
#
\dontrun{
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
trfn = np(paste(extdata_dir, "/Psychotria_5.2.newick", sep=""))
tr = read.tree(trfn)
tmplocation = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
corner_coords(tr, coords_fun="plot_phylo3_nodecoords", tmplocation=tmplocation)
}
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link[ape]{phylo}}, \code{\link{get_nodenums}}
}

