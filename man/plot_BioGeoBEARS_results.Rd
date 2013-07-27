\name{plot_BioGeoBEARS_results}
\alias{plot_BioGeoBEARS_results}
\title{Plot the results of a BioGeoBEARS run}
\usage{
  plot_BioGeoBEARS_results(results_object,
    analysis_titletxt = NULL, addl_params = list(),
    plotwhat = "text", label.offset = NULL, tipcex = 0.8,
    statecex = 0.7, splitcex = 0.6, titlecex = 0.8,
    plotsplits = TRUE, plotlegend = FALSE,
    legend_ncol = NULL, legend_cex = 1,
    cornercoords_loc = "manual", include_null_range = TRUE,
    tr = NULL, tipranges = NULL)
}
\arguments{
  \item{results_object}{The results object from
  \code{\link{bears_optim_run}} (with ancestral states
  on).}

  \item{analysis_titletxt}{The main title of the plot. If
  NULL, \code{results_object$inputs$description} is
  checked.}

  \item{addl_params}{The function will plot the
  log-likelihood (LnL) and the ML values of the free
  parameters. If you want additional parameters plotted,
  list them here.}

  \item{plotwhat}{To plot the ML discrete states, "text".
  To plot a piechart of the relative probability of all the
  states, "pie".}

  \item{label.offset}{Offset for the tree tip labels. If
  \code{NULL}, program chooses 0.05 x tree height.}

  \item{tipcex}{\code{cex} value for the tiplabels (scaling
  factor, i.e. 0.5 is half size)}

  \item{statecex}{\code{cex} value for the states (scaling
  factor, i.e. 0.5 is half size). Used on piecharts if
  plotwhat="pie".}

  \item{splitcex}{\code{cex} value for the splits (scaling
  factor, i.e. 0.5 is half size). Used on piecharts if
  plotwhat="pie".}

  \item{titlecex}{\code{cex} value for the title (scaling
  factor, i.e. 0.5 is half size).}

  \item{plotsplits}{If \code{TRUE}, plot states on the
  corners -- text or pie charts, depending on
  \code{plotwhat}.}

  \item{plotlegend}{If \code{TRUE}, make a (separate) plot
  with a legend giving the colors for each state/range,
  using \code{\link{colors_legend}}.}

  \item{legend_ncol}{The number of columns in the legend.
  If \code{NULL} (default), the function calculates
  \code{floor(sqrt(length(possible_ranges_list_txt) / 2))}
  when the number of states is <=64, and
  \code{sqrt(ceiling(length(possible_ranges_list_txt)))}
  when > 64. Note that when you have hundreds of states,
  there is probably no good way to have a readable legend,
  and it is easier to just rely upon printing the character
  codes for the ML states in the plots, with the colors,
  and users can then see and trace the common colors/states
  by eye.}

  \item{legend_cex}{The cex (character expansion size) for
  the legend.  Defaults to 1, which means the
  \code{\link[graphics]{legend}} function determines the
  size.  The value 2.5 works well for 15 or 16
  states/ranges.}

  \item{cornercoords_loc}{The directory location containing
  the R script \code{plot_phylo3_nodecoords.R}. This
  function, modified from the APE function
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

  \item{include_null_range}{If \code{TRUE} (default), the
  null range is included in calculation of colors. (Safest
  for now.)}

  \item{tr}{Tree to plot on. Default \code{NULL}, which
  means the tree will be read from the file at
  \code{results_object$inputs$trfn}.}

  \item{tipranges}{Tip geography data. Default \code{NULL},
  which means the tree will be read from the file at
  \code{results_object$inputs$geogfn}.}
}
\description{
  This function plots on a tree the highest-probability
  ancestral states (ranges), splits if desired (these are
  the ranges/states just after cladogenesis, and are
  plotted on the corners of a tree), and/or pie charts at
  nodes.  A legend tying the relationship between colors
  and states/ranges is also optionally plotted.
}
\details{
  The legend is plotted on a separate plot, as it is very
  difficult to predict whether or not there will be space
  on any given tree plot.  The utility of the legend is
  also debatable, as \code{plot_BioGeoBEARS_results} plots
  the colors and state/range names directly onto the plot.
  Any legend will get unwieldy above perhaps 32 states,
  which is just 5 areas with no constraints (see
  \code{\link[cladoRcpp]{numstates_from_numareas}}, or type
  \code{numstates_from_numareas(numareas5, maxareas5,
  include_null_range=TRUE)}.

  Note that this assumes that the ancestral states were
  calculated under the global optimum model (rather than
  the local optimum, with the model re-optimized for each
  possible state at each possible node, as done in e.g.
  \code{LAGRANGE}), and that these are marginal
  probabilities, i.e. this is not a joint reconstruction,
  instead it gives the probabilities of states at each
  node.  This will not always be readable as a joint
  reconstruction (it could depict split scenarios that are
  not possible, for instance.)
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
  \url{https://code.google.com/p/lagrange/}

  Matzke_2012_IBS

  ReeSmith2008
}
\seealso{
  \code{\link{get_leftright_nodes_matrix_from_results}},
  \code{\link{corner_coords}},
  \code{\link[ape]{plot.phylo}},
  \code{\link[ape]{plot.phylo}},
  \code{\link[ape]{tiplabels}},
  \code{\link[graphics]{legend}},
  \code{\link[base]{floor}}, \code{\link[base]{ceiling}},
  \code{\link[base]{floor}},
  \code{\link[cladoRcpp]{numstates_from_numareas}},
  \code{\link[base]{system.file}},
  \code{\link[base]{list.files}}
}

