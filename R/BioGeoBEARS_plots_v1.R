

#######################################################
# Functions to make nice plots
#######################################################


#######################################################
# add_statum_boundaries_to_phylo_plot
# Add stratum boundaries to a plot of a time-scaled phylogeny
# Correcting for the offset in the plot
#######################################################
#' Add stratum boundaries to a phylogeny plot
#' 
#' Adds vertical lines, representing stratum boundaries, 
#' to a plot of a phylogeny. 
#' 
#' \bold{Note:} This function asssumes the phylogeny is plotted with
#' tips to the right. 
#' 
#' \bold{Background:} This function may be useful, because plotting vertical lines
#' onto plots generated with \code{APE}'s \code{\link[ape]{plot.phylo}} plots at the 
#' time-points desired will not work. This is because, confusingly, APE's plotting
#' of phylogenies puts the x-axis minimum at 0, and the max at
#' the height of the tree above the root, plus a buffer for 
#' the tip labels.*
#' 
#' Therefore, \code{add_statum_boundaries_to_phylo_plot} calculates the height of 
#' the highest tip and subtracts the \code{timeperiods}, to get the plotted x-values, then 
#' uses abline to plot the lines. (This information may be helpful if you want to 
#' plot your arbitrary lines on top of phylogenies plotted in other orientations.
#' 
#' *This may change further if the \code{phylo} object
#' \code{tr} has a root edge (a branch below the root), so the function might not 
#' draw the lines in the correct place if there is a root edge.  You can tell if your 
#' tree has a root edge by running \code{names(tr)} and seeing if a \code{root.edge} 
#' is listed. Alternatively, look at \code{tr$root.edge}. 
#'
#' @param tr An \code{APE} \code{phylo} (tree) object.
#' @param timeperiods A vector of one or more timeperiods. The default is 1 timeperiod
#'                    at 1 million years (or 1 time unit).
#' @param lty Line type, e.g. "dashed".  See \code{\link{par}}
#' @param col Color of the line, e.g. "grey50" (default)
#' @param plotlines If TRUE (default), the lines are plotted on the current plot. 
#'                  If FALSE, the function just returns the x-axis positions
#'                  of the lines on the phylogeny plot. 

#' @return \code{line_positions_on_plot} The x-axis positions of the lines.
#' @export
#' @seealso \code{\link[base]{plot}}, \code{\link[base]{par}}, \code{\link[ape]{plot.phylo}}
#' @note Go (BioGeo)BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' 
#' @examples
#' 
#' # Loading the default tree
#' trfn = np(paste(addslash(extdata_dir), "Psychotria_5.2.newick", sep=""))
#' tr = read.tree(trfn)
#' 
#' # Get the tree coordinates (APE 5.0 or higher), i.e. the x and y of each 
#' node.
#' trcoords = plot_phylo3_nodecoords_APE5(tr, plot=FALSE, root.edge=TRUE)
#' 
#' # Set reasonable x-limits (unlike the defaults on APE5.0)
#' xlims = c(min(trcoords$xx), 1.42*max(trcoords$xx))
#' 
#' # Plot the tree
#' trplot = plot(tr, cex=1, x.lim=xlims); axisPhylo()
#' 
#' # Add the stratum boundaries
#' add_statum_boundaries_to_phylo_plot(tr, timeperiods=1, lty="dashed", col="gray50", plotlines=TRUE)
#' 
add_statum_boundaries_to_phylo_plot <- function(tr, timeperiods=1, lty="dashed", col="gray50", plotlines=TRUE)
	{
	# Example code for the programmer to run easily
	SETUP='
	# Loading the default tree
	trfn = np(paste(addslash(extdata_dir), "Psychotria_5.2.newick", sep=""))
	tr = read.tree(trfn)
	
	# Get the tree coordinates (APE 5.0 or higher), i.e. the x and y of each node.
	trcoords = plot_phylo3_nodecoords_APE5(tr, plot=FALSE, root.edge=TRUE)
	
	# Set reasonable x-limits (unlike the defaults on APE5.0)
	xlims = c(min(trcoords$xx), 1.42*max(trcoords$xx))
	
	# Plot the tree
	trplot = plot(tr, cex=1, x.lim=xlims); axisPhylo()
	
	# Add the stratum boundaries
	add_statum_boundaries_to_phylo_plot(tr, timeperiods=1, lty="dashed", col="gray50", plotlines=TRUE)
	' ## END SETUP
	
	ntips = length(tr$tip.label)
	tr_table = prt(tr, printflag=FALSE)
	tr_height = tr_table$time_bp[ntips+1]
	line_positions_on_plot = tr_height-timeperiods
	if (plotlines == TRUE)
		{
		abline(v=line_positions_on_plot, lty=lty, col=col)
		} # END if (plotlines == TRUE)
	return(line_positions_on_plot)
	} # END add_statum_boundaries_to_phylo_plot <- function(tr, lty="dashed")


#######################################################
# plot_BioGeoBEARS_results
#######################################################
#' Plot the results of a BioGeoBEARS run
#' 
#' This function plots on a tree the highest-probability ancestral states (ranges), 
#' splits if desired (these are the ranges/states just after cladogenesis, and are 
#' plotted on the corners of a tree), and/or pie charts at nodes/splits.  A legend 
#' tying the relationship between colors and states/ranges is also optionally plotted.
#'
#' The legend, if desired, is plotted on a separate plot, as it is very difficult 
#' to predict whether or not there will be appropriate space on any given tree plot. 
#' The utility of a classical legend is also debatable, as 
#' \code{plot_BioGeoBEARS_results} plots the colors and state/range names directly 
#' onto the plot.  Any legend will get unwieldy above perhaps 16
#' states/ranges, which is just 4 areas with no constraints (see \code{\link[cladoRcpp]{numstates_from_numareas}}, or type \code{numstates_from_numareas(numareas=4, maxareas=4, include_null_range=TRUE)}.
#'
#' In any case, the human eye can only easily read a few colors on a plot. The 
#' philosophy adopted in \code{BioGeoBEARS} plots is to use bright, primary colors for the 
#' single-area ranges, and then blend these colors for multi-areas ranges. A range 
#' all areas is colored white. Coupled with plotting the letter codes for the 
#' ranges ("A", "AB", "ABC"), this makes for reasonably readable plots. (Of course,
#' some researchers design their own custom plots in Adobe Illustrator or elsewhere.)
#' 
#' Note that \code{plot_BioGeoBEARS_results} assumes
#' that the ancestral states were calculated under the global optimum model (rather than the local optimum, with the model re-optimized for each 
#' possible state at each possible node, as done in e.g. \code{LAGRANGE}), and that these are marginal probabilities, i.e. this is not a joint reconstruction,
#' instead it gives the probabilities of states at each node.  This will not always be readable as a joint reconstruction (it could depict split scenarios
#' that are not possible, for instance.)
#' 
#' @param results_object The results object from \code{\link{bears_optim_run}} (with ancestral states on).
#' @param analysis_titletxt The main title of the plot. If NULL, \code{results_object$inputs$description} is checked.
#' @param addl_params The function will plot the log-likelihood (LnL) and the ML values of the free parameters. If you want additional parameters plotted, list them here.
#' @param plotwhat To plot the ML discrete states, "text".  To plot a piechart of the relative probability of all the states, "pie".
#' @param label.offset Offset for the tree tip labels. If \code{NULL}, program chooses 0.05 x tree height.
#' @param tipcex \code{cex} value for the tiplabels (scaling factor, i.e. 0.5 is half size)
#' @param statecex \code{cex} value for the states (scaling factor, i.e. 0.5 is half size). Used on piecharts if plotwhat="pie".
#' @param splitcex \code{cex} value for the splits (scaling factor, i.e. 0.5 is half size). Used on piecharts if plotwhat="pie".
#' @param titlecex \code{cex} value for the title (scaling factor, i.e. 0.5 is half size). 
#' @param plotsplits If \code{TRUE}, plot states on the corners -- text or pie charts, depending on \code{plotwhat}.
#' @param plotlegend If \code{TRUE}, make a (separate) plot with a legend giving the colors for each state/range, using \code{\link{colors_legend}}.
#' @param legend_ncol The number of columns in the legend.  If \code{NULL} (default), the function calculates \code{floor(sqrt(length(possible_ranges_list_txt) / 2))} 
#' when the number of states is <=64, and \code{sqrt(ceiling(length(possible_ranges_list_txt)))} when > 64. Note that when you have hundreds of states, there is probably 
#' no good way to have a readable legend, and it is easier to just rely upon printing the 
#' character codes for the ML states in the plots, with the colors, and users can then see and trace the common colors/states by eye.
#' @param legend_cex The cex (character expansion size) for the legend.  Defaults to 1, which means the \code{\link[graphics]{legend}} function determines the 
#' size.  The value 2.5 works well for 15 or 16 states/ranges.
#' @param cornercoords_loc The directory location containing the R script \code{plot_phylo3_nodecoords.R}. This function, modified from the APE function
#' \code{\link[ape]{plot.phylo}}, cannot be included directly in the R package as it contains C code that does not pass CRAN's R CMD check. The default, 
#' cornercoords_loc="manual", will not allow split states to be plot.  The R script \code{plot_phylo3_nodecoords.R} is located in the BioGeoBEARS extension data 
#' directory, \code{extdata/a_scripts}.  You should be able to get the full path with \code{list.files(system.file("extdata/a_scripts", package="BioGeoBEARS"), full.names=TRUE)}.
#' @param include_null_range If \code{TRUE} (default), the null range is included in calculation of colors. (Safest for now.)
#' @param tr Tree to plot on. Default \code{NULL}, which means the tree will be read from the file at \code{results_object$inputs$trfn}.
#' @param tipranges Tip geography data. Default \code{NULL}, which means the tree will be read from the file at \code{results_object$inputs$geogfn}.
#' @param if_ties What to do with ties in probability. Currently, 
#' the options are:
#' (1) "takefirst", which takes the first tied state in the 
#' probabilities column (The full probabilities of all states will be
#' shown in e.g. pie charts, of course); (2) "asterisk", which 
#' returns returns the first tied state, but marks it with an 
#' asterisk ("*").
#' @param juststats If \code{FALSE} (default), plots are plotted. If \code{TRUE}, 
#' no plots are done, 
#' the function just returns the summary statistics.
#' @param root.edge Should the root edge be plotted, if it exists in the tree?  Passed to
#' plot.phylo().  This can be handy if the root state display is getting cut off.
#' @param colors_list_for_states Default \code{NULL} auto-generates colors with 
#' \code{get_colors_for_states_list_0based}. Otherwise, users can specify colors for 
#' each state (remember that e.g. 4 areas can mean 2^4=16 states).
#' @param skiptree Skip the plotting of the tree -- useful for plotting the state labels
#' e.g. on top of a stochastic map. Default \code{FALSE}.
#' @show.tip.label Same as for APE's \code{plot.phylo}.
#' @tipcol The tip text color. Default "black".
#' @dej_params_row Parameters used to generate an SSE simulation. dej_params_row can be 
#' obtained from \code{get_dej_params_row_from_SSEsim_results_processed}.
#' @plot_max_age The maximum age of the plot (age below the tips). Default
#' is tree height
#' @skiplabels If \code{TRUE}, skip the nodelabels command, resulting in plotting just the 
#' tree. (Yes, you could have just used \code{FALSE}).  Default \code{FALSE}.
#' @plot_stratum_lines If \code{TRUE} (default), plot dotted lines at the stratum 
#' boundaries, *if* it's a time-stratified results object.
#' @simplify_piecharts If \code{TRUE}, just plot one slice for the maximum probability, 
#' and white for "other". This should help with large plots with many ranges, 
#' which can overwhelm some graphics programs (imagine 1000 slices per piechart,
#' times 1000+nodes). Default \code{FALSE}.
#' @tipboxes_TF Plot the tip boxes (tip states)?  Default \code{TRUE}.
#' @tiplabel_adj Justification for tiplabel boxes (same as "adj" parameter 
#' of text()). Default is c(0.5). Left justification: c(0). Right: c(1). 
#' Justification top left: c(0,0), etc.
#' @no.margin Same as in plot.phylo. Default is FALSE (meaning yes, 
#' there will be margins).
#' @xlims Same as in plot.phylo x.lim. Default is basically c(0,treeheight), so
#' if tiplabels are getting cut off, try e.g. xlims=c(0,1.5*treeheight).
#' @ylims Same as in plot.phylo y.lim. Default is basically c(0,numtips), so
#' if the spacing between the title and time axis and the tree are too big,
#' try e.g. ylims=c(0+10,numtips-10). Trial and error will get you there. 
#' See: https://stat.ethz.ch/pipermail/r-sig-phylo/2013-March/002540.html
#' @export
#' @seealso \code{\link{get_leftright_nodes_matrix_from_results}}, \code{\link{corner_coords}}, \code{\link[ape]{plot.phylo}}, \code{\link[ape]{plot.phylo}}, \code{\link[ape]{tiplabels}}, \code{\link[graphics]{legend}}, \code{\link[base]{floor}}, \code{\link[base]{ceiling}}, \code{\link[base]{floor}}, \code{\link[cladoRcpp]{numstates_from_numareas}}, \code{\link[base]{system.file}}, \code{\link[base]{list.files}}
#' @note Go (BioGeo)BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{https://code.google.com/p/lagrange/}
#' @examples
#' test=1
#' 
plot_BioGeoBEARS_results <- function(results_object, analysis_titletxt=NULL, addl_params=list(), plotwhat="text", label.offset=NULL, tipcex=0.8, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, plotlegend=FALSE, legend_ncol=NULL, legend_cex=1, cornercoords_loc="manual", tr=NULL, tipranges=NULL, if_ties="takefirst", pie_tip_statecex=0.7, juststats=FALSE, xlab="Millions of years ago", root.edge=TRUE, colors_list_for_states=NULL, skiptree=FALSE, show.tip.label=TRUE, tipcol="black", dej_params_row=NULL, plot_max_age=NULL, skiplabels=FALSE, plot_stratum_lines=TRUE, include_null_range=NULL, plot_null_range=FALSE, simplify_piecharts=FALSE, tipboxes_TF=TRUE, tiplabel_adj=c(0.5), no.margin=FALSE, xlims=NULL, ylims=NULL)
	{
	
	junk='
	# manual_ranges_txt=NULL, 
	# @manual_ranges_txt If you dont want to use the default text for each range, produced
	# by areas_list_to_states_list_new(), specify the list here.

	
	scriptdir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/a_scripts/"
	plot_BioGeoBEARS_results(results_object, analysis_titletxt=NULL, addl_params=list(), plotwhat="text", label.offset=NULL, tipcex=0.8, statecex=0.8, splitcex=0.8, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=NULL, tipranges=NULL)
	
	# Defaults
	addl_params=list("j"); plotwhat="text"; label.offset=0.45; tipcex=0.7; statecex=0.7; splitcex=0.6; titlecex=0.8; plotsplits=TRUE; cornercoords_loc=scriptdir; include_null_range=TRUE; tr=tr; tipranges=tipranges; juststats = FALSE; plotlegend=FALSE; 	xlab="Millions of years ago"; if_ties="takefirst"
	
	
	# Setup
results_object = resDEC
analysis_titletxt ="BioGeoBEARS DEC on Mariana M1v4_unconstrained"
addl_params=list("j"); plotwhat="text"; label.offset=0.45; tipcex=0.7; statecex=0.7; splitcex=0.6; titlecex=0.8; plotsplits=TRUE; cornercoords_loc=scriptdir; include_null_range=TRUE; tr=tr; tipranges=tipranges
juststats=FALSE; plotlegend=FALSE; 	xlab="Millions of years ago"; if_ties="takefirst"
show.tip.label=TRUE
tipcol="black"; dej_params_row=NULL; plot_max_age=NULL; skiplabels=FALSE; 
colors_list_for_states=NULL
skiptree=FALSE
include_null_range=NULL
plot_stratum_lines=TRUE
	plot_null_range = FALSE
	' # endjunk
	
	# Default; no longer used
	if (is.null(include_null_range) == TRUE)
		{
		include_null_range = results_object$inputs$include_null_range
		}
	
	# Force this in, if user-specified
	results_object$inputs$include_null_range = include_null_range
	
	#######################################################
	# User modifications to border color (externally, using
	# 'par(fg=NA)' or whatever
	#######################################################
	# border color (for modifying this)
	tmp_fg = par("fg")
	par(fg="black")	# set to default for most things

	
	#######################################################
	# Plot ancestral states - DEC
	#######################################################


	# Setup
	#results_object = resDEC_strat
	BioGeoBEARS_run_object = results_object$inputs

	# Read the tree from file, if needed
	if (is.null(tr))
		{
		#tr = read.tree(BioGeoBEARS_run_object$trfn)
		tr = check_trfn(trfn=BioGeoBEARS_run_object$trfn)
		}
	tr_pruningwise = reorder(tr, "pruningwise")
	
	# Basic tree info
	tips = 1:length(tr_pruningwise$tip.label)
	nodes = (length(tr_pruningwise$tip.label)+1):(length(tr_pruningwise$tip.label)+tr_pruningwise$Nnode)


	
	# Read the tipranges from file, if needed.
	if (is.null(tipranges))
		{
		# Get geographic ranges at tips
		if (BioGeoBEARS_run_object$use_detection_model == FALSE)
			{
			tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=np(BioGeoBEARS_run_object$geogfn))
			}
		if (BioGeoBEARS_run_object$use_detection_model == TRUE)
			{
			if (BioGeoBEARS_run_object$use_detection_model == TRUE)
				{
				tipranges = tipranges_from_detects_fn(detects_fn=BioGeoBEARS_run_object$detects_fn)
				} # END if (inputs$use_detection_model == TRUE)
			} # END if (BioGeoBEARS_run_object$use_detection_model == TRUE)
		} # END if (is.null(tipranges))
	
	# Basic areas info
	areas = getareas_from_tipranges_object(tipranges)
	areas

	numareas = length(areas)
	numareas
	
	if (!is.na(results_object$inputs$max_range_size))
		{
		max_range_size = results_object$inputs$max_range_size
		} else {
		max_range_size = length(areas)
		}
	max_range_size


	if (is.null(results_object$inputs$states_list))
		{
		numstates = numstates_from_numareas(numareas=length(areas), maxareas=max_range_size, include_null_range=results_object$inputs$include_null_range)
		numstates
		states_list_areaLetters = areas_list_to_states_list_new(areas, maxareas=max_range_size, include_null_range=results_object$inputs$include_null_range)
		#states_list
		states_list_0based_index = rcpp_areas_list_to_states_list(areas, maxareas=max_range_size, include_null_range=results_object$inputs$include_null_range)
		#states_list_0based_index
		} else {
		states_list_0based_index = results_object$inputs$states_list
		#states_list = states_list_indexes_to_areastxt(states_list=states_list_0based_index, areanames=areas, counting_base=0, concat=FALSE, sep="")
		}


	# calculate the ML marginal probs of states at the base of each branch
	# above each split (local, non-joint probs, global model)
	# 2014-05-15_NJM: Used to add:
	# results_object$ML_marginal_prob_each_split_at_branch_bottom_BELOW_node = ML_marginal_prob_each_split_at_branch_bottom_BELOW_node / rowSums(ML_marginal_prob_each_split_at_branch_bottom_BELOW_node)
	# ... but this is now totally pointless since this is done automatically
	#results_object = get_MLsplitprobs_from_results(results_object)
	#names(results_object)

	# Extract ML parameter values, and LnL
	# This will work with optim, optimx2012, or optimx2013
	
	# Handy summary outputs
	param_ests = extract_params_from_BioGeoBEARS_results_object(results_object, returnwhat="table", addl_params=addl_params, paramsstr_digits=4)

	
	# If you want to skip the plotting and just extract
	# the parameter values
	if (juststats == TRUE)
		{
		return(param_ests)		
		} else {
		paramstr = extract_params_from_BioGeoBEARS_results_object(results_object, returnwhat="string", addl_params=addl_params, paramsstr_digits=4)
		} # if (juststats == TRUE)

	# Get the parameter names
	param_names = extract_params_from_BioGeoBEARS_results_object(results_object, returnwhat="param_names", addl_params=addl_params, paramsstr_digits=4)

		
	# Major title (short description)
	if (is.null(analysis_titletxt))
		{
		tmptxt = results_object$inputs$description
		if (any(is.null(tmptxt), tmptxt=="", tmptxt=="defaults", tmptxt=="default"))
			{
			analysis_titletxt = ""
			} else {
			analysis_titletxt = results_object$inputs$description
			}
		}
	
	
	if (is.null(dej_params_row))
		{
		# Default text for an inference or stochastic mapping
		analysis_titletxt = paste(analysis_titletxt, "\n", "ancstates: global optim, ", max_range_size, " areas max. ", paramstr, sep="")
		analysis_titletxt
		} else {
		# Text for an SSE simulation
		dej_params_row

		brate_col_TF = names(dej_params_row) == "brate"
		brate_col = (1:length(dej_params_row))[brate_col_TF]
		biogeog_params = dej_params_row[1:(brate_col-1)]
		biogeog_param_names = names(dej_params_row)[1:(brate_col-1)]
		equals_col = "="
		
		tmpcols = cbind(biogeog_param_names, equals_col, unlist(biogeog_params))
		tmpcols
		txtrows = apply(X=tmpcols, MARGIN=1, FUN=paste, sep="", collapse="")
		txtrows
		biogeog_params_txt = paste(txtrows, sep="", collapse="; ")
		biogeog_params_txt
			
		titletxt2 = bquote(paste(.(max_range_size), " areas max., ", .(biogeog_params_txt), "; ", lambda, "=", .(dej_params_row$brate), "; ",  mu, "=", .(dej_params_row$drate), "; ", alpha, "=", .(dej_params_row$rangesize_b_exponent), "; ", omega, "=", .(dej_params_row$rangesize_d_exponent), "", sep=""))
		
		#print(titletxt2)
		
		} # END if (is.null(dej_params_row))




	#######################################################
	# Get the marginal probs of the splits (global ML probs, not local)
	# (These are marginal, rather than joint probs; but not local optima)
	#######################################################
	leftright_nodes_matrix = get_leftright_nodes_matrix_from_results(tr_pruningwise)

	# This gets you the prob. of each state at the left base above each node, and
	# at the right base above each node
	marprobs = results_object$ML_marginal_prob_each_state_at_branch_bottom_below_node
	left_ML_marginals_by_node = marprobs[leftright_nodes_matrix[, 2], ]
	right_ML_marginals_by_node = marprobs[leftright_nodes_matrix[, 1], ]
	right_ML_marginals_by_node

	# If they aren't matrices (because of a 2-species, 1-internal-node tree), fix that
	if (is.null(dim(left_ML_marginals_by_node)))
		{
		left_ML_marginals_by_node = matrix(data=left_ML_marginals_by_node, nrow=1)
		}
	if (is.null(dim(right_ML_marginals_by_node)))
		{
		right_ML_marginals_by_node = matrix(data=right_ML_marginals_by_node, nrow=1)
		}


	#######################################################
	# Extract the outputs ancestral states at nodes, and plot!
	#######################################################
	relprobs_matrix = results_object$ML_marginal_prob_each_state_at_branch_top_AT_node
	
	if (length(nodes) > 1)
		{
		relprobs_matrix_for_internal_states = relprobs_matrix[nodes,]	# subset to just internal nodes
		} else {
		relprobs_matrix_for_internal_states = relprobs_matrix[nodes,]	# subset to just internal nodes
		# Convert back to a matrix
		relprobs_matrix_for_internal_states = matrix(data=relprobs_matrix_for_internal_states, nrow=1, ncol=ncol(relprobs_matrix))
		}
	
	
	relprobs_matrix
	
	if (is.null(states_list_0based_index))
		{
		statenames = areas_list_to_states_list_new(areas, maxareas=max_range_size, include_null_range=results_object$inputs$include_null_range, split_ABC=FALSE)
		ranges_list = as.list(statenames)
		statenames
		} else {
		ranges_list = states_list_0based_to_ranges_txt_list(state_indices_0based=states_list_0based_index, areanames=areas)
		ranges_list
		statenames = unlist(ranges_list)
		statenames
		}


	MLprobs = get_ML_probs(relprobs_matrix)
	MLprobs
	MLstates = get_ML_states_from_relprobs(relprobs_matrix, statenames, returnwhat="states", if_ties=if_ties)


	# Set up colors for each state
	if (is.null(colors_list_for_states))
		{
		# Fix plot_null_range to FALSE (don't want to plot that color)
		colors_matrix = get_colors_for_numareas(length(areas))
		colors_list_for_states = mix_colors_for_states(colors_matrix, states_list_0based_index, plot_null_range=results_object$inputs$include_null_range)
		colors_list_for_states
		} # END if (is.null(colors_list_for_states))

	# Set up colors by possible ranges
	if (is.null(ranges_list))
		{
		possible_ranges_list_txt = areas_list_to_states_list_new(areas,  maxareas=max_range_size, split_ABC=FALSE, include_null_range=results_object$inputs$include_null_range)
		} else {
		possible_ranges_list_txt = ranges_list
		} # if (is.null(ranges_list))
	#possible_ranges_list_txt = areas_list_to_states_list_new(areas,  maxareas=max_range_size, split_ABC=FALSE, include_null_range=include_null_range)

# 	if (plot_null_range == FALSE)
# 		{
# 		possible_ranges_list_txt[possible_ranges_list_txt == "_"] = NULL
# 		possible_ranges_list_txt[possible_ranges_list_txt == ""] = NULL
# 		possible_ranges_list_txt[possible_ranges_list_txt == "NA"] = NULL
# 		possible_ranges_list_txt[is.na(possible_ranges_list_txt)] = NULL
# 		possible_ranges_list_txt[is.null(possible_ranges_list_txt)] = NULL
# 		}

	cols_byNode = rangestxt_to_colors(possible_ranges_list_txt, colors_list_for_states, MLstates)

	# Legend, if desired
	if (plotlegend == TRUE)
		{
		colors_legend(possible_ranges_list_txt, colors_list_for_states, legend_ncol=legend_ncol, legend_cex=legend_cex)
		}
	
	
	# Put in a 0 for root.edge
	if (root.edge == FALSE)
		{
		tr$root.edge = 0
		} # END if (root.edge == FALSE)
	if (root.edge == TRUE)
		{
		if (is.null(tr$root.edge) == TRUE)
			{
			tr$root.edge = 0
			} # END if (is.null(tr$root.edge) == TRUE)
		} # END if (root.edge == TRUE)

	# Default label offset
	if (is.null(label.offset))
		{
		label.offset = 0.007 * (get_max_height_tree(tr) + tr$root.edge)
		}
		
	# Manual changing of xlims for phylogeny plot, with plot_max_age
	if (show.tip.label == TRUE)
		{
		# If plot_max_age *IS NOT* specified
		if (is.null(plot_max_age))
			{
			max_x = 1.25*(get_max_height_tree(tr) + tr$root.edge)
			min_x = 0
			} else {
			# If plot_max_age *IS* specified
			nontree_part_of_x = plot_max_age - (get_max_height_tree(tr) + tr$root.edge)
			max_x = 1.25*(get_max_height_tree(tr) + tr$root.edge)
			min_x = -1 * nontree_part_of_x
			}
		} else {
		# NO tip labels
		if (is.null(plot_max_age))
			{
			# If plot_max_age *IS NOT* specified
			max_x = 1.05*(get_max_height_tree(tr) + tr$root.edge)
			min_x = 0
			} else {
			# If plot_max_age *IS* specified
			nontree_part_of_x = plot_max_age - (get_max_height_tree(tr) + tr$root.edge)
			max_x = 1.05*(get_max_height_tree(tr) + tr$root.edge)
			min_x = -1 * nontree_part_of_x
			}
		} # if (show.tip.label == TRUE)


	###################################################
	# Calculate x-axis ticks (axisPhylo alternative)
	###################################################	
	max_tree_x = 1.0 * (get_max_height_tree(tr) + tr$root.edge)
	
	# Plot min/max
	if (is.null(xlims))
		{
		xlims = c(min_x, max_x)
		} else {
		xlims = xlims
		}
	#print(xlims)

	# Tree min/max
	nodecoords = node_coords(tr, tmplocation=cornercoords_loc, root.edge=root.edge)
	max_tree_x = max(nodecoords$x)
	
	# AxisPhylo() min/max
	if (is.null(plot_max_age))
		{
		xticks_desired_lims = c(0, max_tree_x)
		} else {
		xticks_desired_lims = c(0, plot_max_age)
		}

	xticks_desired = pretty(xticks_desired_lims)
	
	# Translate into plot coordinates
	xaxis_ticks_locs = max_tree_x - xticks_desired
	
	#print(xticks_desired)
	#print(xaxis_ticks_locs)
	###################################################	
		
	# Skip tree plotting if it has already been done:
	if (skiptree != TRUE)
		{
		#######################################################
		# Otherwise, plot the tree!!
		#######################################################
		plot(tr_pruningwise, x.lim=xlims, y.lim=ylims, show.tip.label=FALSE, label.offset=label.offset, cex=tipcex, no.margin=no.margin, root.edge=root.edge)
		
		if (show.tip.label == TRUE)
			{
			tiplabels_to_plot = sapply(X=tr_pruningwise$tip.label, FUN=substr, start=1, stop=30)
			if (skiplabels == FALSE)
				{
				tiplabels(text=tiplabels_to_plot, tip=tips, cex=tipcex, adj=0, bg="white", frame="n", pos=4, offset=label.offset, col=tipcol)	# pos=4 means labels go to the right of coords
				} # END if (skiplabels == FALSE)
			} # END if (show.tip.label == TRUE)
		
		#axisPhylo(cex.axis=titlecex)
		axis(side=1,  at=xaxis_ticks_locs, label=xticks_desired)
		# Plot the title
		mtext(text=xlab, side=1, line=2, cex=titlecex)
		}
	
	# Add states / piecharts
	if (plotwhat == "text")
		{
		# Use statecex for pie chart size at nodes AND states at tips
		par(fg=tmp_fg)	# so that user can change border externally
		
		if (skiplabels == FALSE)
			{
			nodelabels(text=MLstates[nodes], node=nodes, bg=cols_byNode[nodes], cex=statecex)		
			tiplabels(text=MLstates[tips], tip=tips, bg=cols_byNode[tips], cex=statecex, adj=tiplabel_adj)
			} # END if (skiplabels == FALSE)
		
		par(fg="black")	# set to default for most things
		}
	if (plotwhat == "pie")
		{
		# Use statecex for pie chart size at nodes BUT for states at tips,
		# use "pie_tip_statecex"
		par(fg=tmp_fg)	# so that user can change border externally

		if (skiplabels == FALSE)
			{
			# DOSIMPLIFY PIE CHARTS
			if (simplify_piecharts == TRUE)
				{
				# columns to keep in the final
				colnums_to_keep_in_probs = NULL
				
				# Probs table
				probs = results_object$ML_marginal_prob_each_state_at_branch_top_AT_node
				probs2 = probs
				maxprob = rep(0, nrow(probs))
				other = rep(0, nrow(probs))
				num_to_keep = 1
				cat("\nSince simplify_piecharts==TRUE, reducing prob pie charts to (most probable, other)...\n")
				for (i in 1:nrow(probs))
					{
					cat(i, " ", sep="")
					tmprow = probs[i,]
					positions_highest_prob_to_lowest = rev(order(tmprow))
					# If there are ties, we take the first ones
					positions_to_keep = positions_highest_prob_to_lowest[1:num_to_keep]
					colnums_to_keep_in_probs = c(colnums_to_keep_in_probs, positions_to_keep)
					keepTF = rep(FALSE, length(tmprow))
					keepTF[positions_to_keep] = TRUE
	
					# Sum the others
					otherTF = keepTF == FALSE
					other[i] = sum(tmprow[otherTF])
					tmprow[otherTF] = 0
					probs2[i,] = tmprow
					} # END for (i in 1:nrow(probs))
				cat("\n")
				
				colnums_to_keep_in_probs_in_order = sort(unique(colnums_to_keep_in_probs))
				probs3 = cbind(probs2[,colnums_to_keep_in_probs_in_order], other)
				# Subset to just internal nodes
				probs3 = probs3[nodes,]
				
				newcols = c(colors_list_for_states[colnums_to_keep_in_probs_in_order], "white")

				# DO SIMPLIFY PIE CHARTS
				nodelabels(pie=probs3, node=nodes, piecol=newcols, cex=statecex)
				} else {
				# DON'T SIMPLIFY PIE CHARTS
				nodelabels(pie=relprobs_matrix_for_internal_states, node=nodes, piecol=colors_list_for_states, cex=statecex)
				} # END if (simplify_piecharts == TRUE)

			# Plot the tiplabels, if desired
			if (tipboxes_TF == TRUE)
				{
				tiplabels(text=MLstates[tips], tip=tips, bg=cols_byNode[tips], cex=pie_tip_statecex, adj=tiplabel_adj)
				} # END if (tipboxes_TF = TRUE)
			} # END if (skiplabels == FALSE)

		par(fg="black")	# set to default for most things
		}
	
	if (skiptree != TRUE)
		{
		if (titlecex > 0)
			{
			#par(ps = 12, cex = titlecex, cex.main = titlecex)
			par(cex.main = titlecex)
			title(analysis_titletxt)
			if (!is.null(dej_params_row))
				{
				# Subtitle for SSE simulation plots
				title(titletxt2, line=1)
				#print(titletxt2)
				}
			#par("font.main") = 2
			}
		}

	if (plotsplits == TRUE)
		{
		# Error check; users must specify the location of the function "plot_phylo3_nodecoords"
		if (cornercoords_loc == "manual")
			{
			stoptxt = cat("\nNOTE: To plot splits, this function needs to access the function 'plot_phylo3_nodecoords'.\n",
							"The function is modified from an APE function, and cannot be directly included in the package,\n",
							"due to some C code that does not meet CRAN standards. To solve this, give plot_BioGeoBEARS_results\n",
							"a 'cornercoords_loc' string that gives the directory of plot_phylo3_nodecoords.R.  Typically this\n",
							"can be found via: ", 'tmp=np(system.file("extdata/a_scripts", package="BioGeoBEARS"))\n',
							"then: list.files(tmp); print(tmp)\n", sep="")
			plotsplits = FALSE
			}
		}

	if (plotsplits == TRUE)
		{
		#######################################################
		# Also add the splits to the plot
		#######################################################
		# First, get the corner coordinates
		coords = corner_coords(tr, tmplocation=cornercoords_loc, root.edge=root.edge)
		coords

		# LEFT SPLITS
		relprobs_matrix = left_ML_marginals_by_node
		
		if (plotwhat == "text")
			{
			MLprobs = get_ML_probs(relprobs_matrix)
			MLprobs
			MLstates = get_ML_states_from_relprobs(relprobs_matrix, statenames, returnwhat="states", if_ties=if_ties)
			MLstates
			length(MLstates)
		
			# Set up colors
			#possible_ranges_list_txt = areas_list_to_states_list_new(areas,  maxareas=max_range_size, split_ABC=FALSE, include_null_range=include_null_range)
			cols_byNode = rangestxt_to_colors(possible_ranges_list_txt, colors_list_for_states, MLstates)
			
			par(fg=tmp_fg)	# so that user can change border externally
			if (skiplabels == FALSE)
				{			
				cornerlabels(text=MLstates, coords=coords$leftcorns, bg=cols_byNode, cex=splitcex)
				} # END if (skiplabels == FALSE)

			par(fg="black")	# set to default for most things
			}
		
		if (plotwhat == "pie")
			{
			par(fg=tmp_fg)	# so that user can change border externally
			cornerpies(pievals=relprobs_matrix, coords$leftcorns, piecol=colors_list_for_states, cex=splitcex)
			par(fg="black")	# set to default for most things
			}



		# RIGHT SPLITS
		relprobs_matrix = right_ML_marginals_by_node

		if (plotwhat == "text")
			{
			MLprobs = get_ML_probs(relprobs_matrix)
			MLprobs
			MLstates = get_ML_states_from_relprobs(relprobs_matrix, statenames, returnwhat="states", if_ties=if_ties)
			MLstates
			length(MLstates)

			# Set up colors
			#possible_ranges_list_txt = areas_list_to_states_list_new(areas,  maxareas=max_range_size, split_ABC=FALSE, include_null_range=include_null_range)
			cols_byNode = rangestxt_to_colors(possible_ranges_list_txt, colors_list_for_states, MLstates)

			par(fg=tmp_fg)	# so that user can change border externally
			if (skiplabels == FALSE)
				{
				cornerlabels(text=MLstates, coords=coords$rightcorns, bg=cols_byNode, cex=splitcex)
				} # END if (skiplabels == FALSE)
			par(fg="black")	# set to default for most things
			}

		if (plotwhat == "pie")
			{
			par(fg=tmp_fg)	# so that user can change border externally
			cornerpies(pievals=relprobs_matrix, coords$rightcorns, piecol=colors_list_for_states, cex=splitcex)			
			par(fg="black")	# set to default for most things
			}
		}

	# Plot vertical dashed lines for timeperiods
	# Is it time-stratified? Plot lines if desired
	if ( ((is.null(BioGeoBEARS_run_object$timeperiods) == FALSE)) && (plot_stratum_lines==TRUE) )
		{
		timeperiods = BioGeoBEARS_run_object$timeperiods
		line_positions_on_plot = add_statum_boundaries_to_phylo_plot(tr, timeperiods=timeperiods, lty="dashed", col="gray50", plotlines=TRUE)
		} # END plot vertical dashed lines for timeperiods



	# Handy summary outputs
	param_ests = matrix(data=param_ests, nrow=1)
	param_ests = adf2(param_ests)
	
	param_ests = dfnums_to_numeric(param_ests)
	names(param_ests) = c("LnL", "nparams", param_names)
	
	return(param_ests)
	} # END plot_BioGeoBEARS_results






#######################################################
# Colors and legends
#######################################################


#######################################################
# get_colors_for_numareas
#######################################################
#' Get colors for a certain number of single areas
#' 
#' Like it says.
#' 
#' @param numareas The number of areas
#' @param use_rainbow If TRUE, force use of \code{rainbow()}
#' @return \code{colors_matrix} The colors for the single areas, 1 column per area
#' @export
#' @seealso \code{\link[stats]{optim}}
#' @note Go (BioGeo)BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' 
get_colors_for_numareas = function(numareas, use_rainbow=FALSE)
	{
	if (use_rainbow == TRUE)
		{
		area_colors = rainbow(numareas)
		colors_matrix = col2rgb(area_colors)
		return(colors_matrix)
		}
	
	if (numareas == 2)
		{
		area_colors = c("blue", "green")
		}
	if (numareas == 3)
		{
		area_colors = c("blue", "green", "red")
		}
	if (numareas == 4)
		{
		area_colors = c("blue", "green3", "yellow", "red")
		}
	if (numareas == 5)
		{
		area_colors = c("blue", "cyan", "green3", "yellow", "red")
		}
	if (numareas == 6)
		{
		area_colors = c("blue", "cyan", "green3", "yellow", "red", "violet")
		}
	if (numareas == 7)
		{
		area_colors = c("blue", "cyan", "green3", "yellow", "orange", "red", "violet")
		}
	if (numareas == 7)
		{
		area_colors = c("blue", "cyan", "green3", "yellow", "orange", "red", "violet")
		}
	if (numareas > 7)
		{
		area_colors = col2rgb(rainbow(numareas))
		return(area_colors)
		}
	
	colors_matrix = col2rgb(area_colors)
	return(colors_matrix)	
	}





#######################################################
# mix_colors_for_states
#######################################################
#' Mix colors logically to produce colors for multi-area ranges
#' 
#' Like it says.
#' 
#' @param colors_matrix A column with a color for each single area
#' @param states_list_0based_index States list giving areas, 0-based
#' @param plot_null_range If FALSE, null ranges are excluded (however coded). 
#' Default FALSE, and should be false unless you *really* want to plot null range.
#' @return \code{colors_list_for_states} The colors for the ML states
#' @export
#' @seealso \code{\link[stats]{optim}}
#' @note Go (BioGeo)BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' 
mix_colors_for_states <- function(colors_matrix, states_list_0based_index, plot_null_range=FALSE)
	{
	# Eliminate the null range, if present
	if (plot_null_range == FALSE)
		{
		states_list_0based_index[states_list_0based_index == "_"] = NULL
		states_list_0based_index[states_list_0based_index == ""] = NULL
		states_list_0based_index[states_list_0based_index == "NA"] = NULL
		states_list_0based_index[is.na(states_list_0based_index)] = NULL
		states_list_0based_index[is.null(states_list_0based_index)] = NULL
		}
	
	
	colors_list_for_states = rep(NA, length(states_list_0based_index))   #matrix(data=NA, nrow=3, ncol=length(states_list_0based_index))
	
	# There is a column for each single area
	numareas = ncol(colors_matrix)
	
	
	# Combine the colors of the input areas, and, if 2 areas or more, divide by (numareas-0.5)
	# (assures no duplication of colors)
	for (i in 1:length(states_list_0based_index))
		{
		tmpareas_in_state_0based_index = states_list_0based_index[[i]]
		tmp_numareas = length(tmpareas_in_state_0based_index)
		
		# Check for null etc.
		if (tmpareas_in_state_0based_index == "_" || tmpareas_in_state_0based_index == "" || is.na(tmpareas_in_state_0based_index) || is.null(tmpareas_in_state_0based_index))
			{
			# NULL/empty range gets black
			colors_list_for_states[i] = rgb(red=0, green=0, blue=0, maxColorValue=255)
			next()
			}
		
		
		
		# Fill in the colors for this area
		
		# Make it white, if its all areas
		if (tmp_numareas == numareas)
			{
			# input white
			colors_list_for_states[i] = rgb(red=255, green=255, blue=255, maxColorValue=255)
			next()
			}
		
		# Otherwise
		sum_color_nums = rep(0, 3)
		for (j in 1:tmp_numareas)
			{
			tmparea_1based_index = 1 + tmpareas_in_state_0based_index[j]
			
			color_nums = colors_matrix[ ,tmparea_1based_index]
			sum_color_nums = sum_color_nums + color_nums
			}
		
		if (tmp_numareas >= 2)
			{
			divide_by = (tmp_numareas - 0.0)
			} else {
			divide_by = 1
			}
		
		sum_color_nums = round(sum_color_nums / divide_by)
		
		# Save the colors as hex format
		colors_list_for_states[i] = rgb(
			red=sum_color_nums[1], 
			green=sum_color_nums[2], 
			blue=sum_color_nums[3], 
			maxColorValue=255)
		}
	
	return(colors_list_for_states)
	}



# Convert a list of ranges text (KOM, MH, KOMIH, etc.) 
#######################################################
# rangestxt_to_colors
#######################################################
#' Convert a list of ranges text (KOM, MH, KOMIH, etc.) 
#' 
#' Like it says.
#' 
#' @param possible_ranges_list_txt A list of the allowed ranges/states
#' @param colors_list_for_states The corresponding colors
#' @param MLstates The ML states for the internal nodes
#' @return \code{MLcolors} The colors for the ML states
#' @export
#' @seealso \code{\link[stats]{optim}}
#' @note Go (BioGeo)BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' 
rangestxt_to_colors <- function(possible_ranges_list_txt, colors_list_for_states, MLstates)
	{
	if (length(possible_ranges_list_txt) != length(colors_list_for_states))
		{
		cat("\nlength(possible_ranges_list_txt) = ", length(possible_ranges_list_txt))
		cat("\nlength(colors_list_for_states) = ", length(colors_list_for_states))
		cat("\nlength(MLstates) = ", length(MLstates))
		
		
		stop("Error: possible_ranges_list_txt and colors_list_for_states must be the same length")
		}
	
	ranges_colors = 1:length(MLstates)
	
	# Nums
	#nums = 1:length(possible_ranges_list_txt)
	

	match_indices = get_indices_where_list1_occurs_in_list2(list1=MLstates, list2=possible_ranges_list_txt)
	MLcolors = colors_list_for_states[match_indices]
	
	return(MLcolors)
	}




#######################################################
# colors_legend
#######################################################
#' Plot a colors legend for geographic ranges
#' 
#' Like it says.
#' 
#' @param possible_ranges_list_txt A list of the allowed ranges/states
#' @param colors_list_for_states The corresponding colors
#' @param legend_ncol The number of columns in the legend.  If \code{NULL} (default), the function calculates \code{floor(sqrt(length(possible_ranges_list_txt) / 2))}.
#' Note that when you have hundreds of states, there is probably no good way to have a coherent legend, and it is easier to just rely upon printing the 
#' character codes for the ML states in the plots, with the colors, and users can then see and trace the common colors/states by eye.
#' @param legend_cex The cex (character expansion size) for the legend.  Defaults to 1, which means the \code{\link[graphics]{legend}} function determines the 
#' size.  The value 2.5 works well for 15 or 16 states/ranges.
#' @location The "location" input goes to the legend() command. Default is "top".
#' @make_blank_plot_first If TRUE (default), colors_legend() will first plot a blank plot
#' with the x- and y-axes ranging from 1-10. make_blank_plot_first=FALSE may be useful for
#' adding the legend on top of other plots.
#' @return Nothing
#' @export
#' @seealso \code{\link[graphics]{legend}}, \code{\link[base]{floor}}, \code{\link[base]{ceiling}}, \code{\link[base]{floor}}
#' @note Go (BioGeo)BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' 
colors_legend <- function(possible_ranges_list_txt, colors_list_for_states, legend_ncol=NULL, legend_cex=1, location="top", make_blank_plot_first=TRUE, ...)
	{

	if (is.null(legend_ncol))
		{
		tmp_numstates = length(possible_ranges_list_txt)
		if (tmp_numstates <= 64)
			{
			legend_ncol = floor(sqrt(tmp_numstates / 2))
			legend_ncol
			} else {
			# For huge numbers of states, just do a square
			legend_ncol = ceiling(sqrt(tmp_numstates))
			}
		}
	
	# Plot, no borders (bty="n"), no labels (xlab, ylab), no tick marks (xaxt, yaxt)
	if (make_blank_plot_first == TRUE)
		{
		plot(1:10, 1:10, pch=".", col="white", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
		}
	#lines(1:4,4:1, col="blue") 
	#legend("top", leg=c("a","b"),col=c("black","blue"), fill=TRUE) 
	legend(x=location, legend=possible_ranges_list_txt, fill=colors_list_for_states, ncol=legend_ncol, title="Legend", cex=legend_cex, ...)#, fill=TRUE) 
	
	} # END colors_legend <- function(possible_ranges_list_txt, colors_list_for_states, legend_ncol=NULL, legend_cex=1, ...)










#######################################################
# map_LGpy_MLsplits_to_tree
#######################################################
#' Take the table of ML splits and node number and map on tree (Python version)
#' 
#' Given a table of splits probabilities from \code{\link{LGpy_splits_fn_to_table}}, map the splits on the tree.
#' 
#' See \code{\link{get_lagrange_nodenums}} for connecting these node numbers to APE node numbers.
#' 
#' @param MLsplits_LGpy A data.frame containing the node numbers, splits, and split probabilities.
#' @param tr An ape phylo object
#' @param tiprange_names The geographic ranges at the tips (i.e. the input data)
#' @return \code{MLsplits_LGpy} A data.frame containing the node numbers, ML splits, and split probabilities; reordered for this plot
#' @export
#' @seealso \code{\link{get_lagrange_nodenums}}, \code{\link{LGpy_splits_fn_to_table}}, \code{\link{LGcpp_splits_fn_to_table}}
#' @note Go (BioGeo)BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @examples
#' test=1
#' 
map_LGpy_MLsplits_to_tree <- function(MLsplits_LGpy, tr, tiprange_names)
	{
	# Order in LAGRANGE order
	MLsplits_LGpy = MLsplits_LGpy[order(MLsplits_LGpy$nodenum_LGpy),]
	MLsplits_LGpy
	
	# Get the names of the columns in the MLsplits table
	tmpnames = names(MLsplits_LGpy)

	# Add R node numbering, if not already done
	if (("Rnodes" %in% tmpnames) == FALSE)
		{
		downpass_node_matrix = get_lagrange_nodenums(tr)
		#downpass_node_matrix[,2] = 1:18
		#downpass_node_matrix = downpass_node_matrix[order(downpass_node_matrix[,1]), ]
	
		MLsplits_LGpy = cbind(MLsplits_LGpy, downpass_node_matrix)
		names(MLsplits_LGpy) = c(tmpnames, "Rnodes", "LGnodes")
		MLsplits_LGpy
		}
	
	MLsplits_LGpy = MLsplits_LGpy[order(MLsplits_LGpy$Rnodes), ]
	MLsplits_LGpy

	# Plot them	
	par(mfrow=c(2,1))
	plot(tr, label.offset=0.15)
	nodelabels(text=MLsplits_LGpy$splits, node=20:37)
	tiplabels(tiprange_names)
	title("LAGRANGE (python) ML splits")
	axisPhylo()
	mtext(text="million years ago", side=1, line=2)
	
	plot(tr, label.offset=0.15)
	pievals = as.matrix(MLsplits_LGpy[,c("relprob","relprob2")])
	nodelabels(pie=pievals, piecol=c("blue", "white"))
	tiplabels(tiprange_names)
	title("LAGRANGE (python) split probs")
	axisPhylo()
	mtext(text="million years ago", side=1, line=2)

	return(MLsplits_LGpy)
	}



#######################################################
# map_LG_MLsplits_to_tree
#######################################################
#' Take the table of ML splits and node number and map on tree (C++ LAGRANGE version)
#' 
#' Given a table of splits probabilities from \code{\link{LGcpp_splits_fn_to_table}}, map the splits on the tree.
#' 
#' See \code{\link{get_lagrange_nodenums}} for connecting these node numbers to APE node numbers.
#' 
#' @param MLsplits_LGcpp A data.frame containing the node numbers, splits, and split probabilities.
#' @param tr An ape phylo object
#' @param tiprange_names The geographic ranges at the tips (i.e. the input data)
#' @param removechar The character to remove, if needed.
#' @param type The type of LAGRANGE input (default C++)
#' @return \code{MLsplits_LGcpp} A data.frame containing the node numbers, ML splits, and split probabilities; reordered for this plot.
#' @export
#' @seealso \code{\link{get_lagrange_nodenums}}, \code{\link{LGpy_splits_fn_to_table}}, \code{\link{LGcpp_splits_fn_to_table}}
#' @note Go (BioGeo)BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @examples
#' test=1
#' 
map_LG_MLsplits_to_tree <- function(MLsplits_LGcpp, tr, tiprange_names, removechar=NULL, type="C++")
	{
	defaults='
	MLsplits_LGpy, tr=tr, tiprange_names=tiprange_names, type="python"
	'
	# Order in LAGRANGE order
	MLsplits_LGcpp = order_LGnodes(MLsplits_LGcpp, tr, removechar, type)
	
	
	# Plot them	
	par(mfrow=c(2,1))
	plot(tr, label.offset=0.15)
	nodelabels(text=MLsplits_LGcpp$splits, node=20:37)
	tiplabels(tiprange_names)
	title(paste("LAGRANGE (", type, ") ML splits", sep=""))
	axisPhylo()
	mtext(text="million years ago", side=1, line=2)
	
	plot(tr, label.offset=0.15)
	pievals = as.matrix(MLsplits_LGcpp[,c("relprob","relprob2")])
	nodelabels(pie=pievals, piecol=c("blue", "white"))
	tiplabels(tiprange_names)
	title(paste("LAGRANGE (", type, ") split probs", sep=""))
	axisPhylo()
	mtext(text="million years ago", side=1, line=2)

	return(MLsplits_LGcpp)
	}


#######################################################
# get_statesColors_table
#######################################################
#' Make a color table for each area and their combinations
#' 
#' Given a list of areas, make a color table for the various combinations.
#' 
#' @param areanames A list of the area names.
#' @param plot_null_range If FALSE, the null range is excluded from the state space. 
#' Note: The null range is never plotted, regardless of the setting, but the function 
#' needs to know which, for proper counting/coloring.
#' @return \code{statesColors_table} A table giving the colors for each state.
#' @export
#' @seealso \code{\link{get_lagrange_nodenums}}, \code{\link{LGpy_splits_fn_to_table}}, \code{\link{LGcpp_splits_fn_to_table}}
#' @note Go (BioGeo)BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @examples
#' test=1
#' 
get_statesColors_table <- function(areanames=c("K","O","M","H"), plot_null_range=FALSE)
	{
	# Make the color matrix for the individual areas
	colors_matrix = get_colors_for_numareas(numareas=length(areanames), use_rainbow=FALSE)
	colors_matrix
	
	# Get the states
	# Here, include_null_range is fixed to FALSE, so that this state is not displayed
	states_list_0based_index = rcpp_areas_list_to_states_list(areas=areanames, include_null_range=plot_null_range)
	
	colors_list_for_states = mix_colors_for_states(colors_matrix, states_list_0based_index, plot_null_range=plot_null_range)
	colors_list_for_states

	possible_ranges_list_txt = states_list_indexes_to_areastxt(states_list=states_list_0based_index, areanames, counting_base=0, sep="")
	possible_ranges_list_txt
	
	statesColors_table = adf2(cbind(possible_ranges_list_txt, colors_list_for_states))
	names(statesColors_table) = c("range", "color")
	return(statesColors_table)
	}


#######################################################
# map_LG_MLsplits_to_tree_corners
#######################################################
#' Map splits to the corners on a phylogeny
#' 
#' What it says.
#' 
#' @param MLsplits A data.frame containing the node numbers, splits, and split probabilities.
#' @param tr An ape phylo object
#' @param tipranges Tipranges object
#' @param removechar The character to remove, if needed.
#' @param type The type of LAGRANGE input (default C++)
#' @param statesColors_table If not default, a table with a color for each area combination.
#' @param bgcol The background color
#' @param areanames The area names, if different from those in the tipranges object
#' @param newplot Default TRUE; should there be a new plot, or should the splits be added to another plot?
#' @param ... Additional arguments to standard functions
#' @return \code{MLsplits} The splits table, ordered appropriately.
#' @export
#' @seealso \code{\link{get_lagrange_nodenums}}, \code{\link{LGpy_splits_fn_to_table}}, \code{\link{LGcpp_splits_fn_to_table}}
#' @note Go (BioGeo)BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @examples
#' test=1
#' 
map_LG_MLsplits_to_tree_corners <- function(MLsplits, tr, tipranges, removechar=NULL, type="C++", statesColors_table="default", bgcol="green3", areanames="default", newplot=TRUE, root.edge, ...)
	{
	defaults='
	MLsplits=MLsplits_LGpy
	type="python"
	'
	# Get corner coordinates
	corners_list = corner_coords(tr, root.edge=root.edge)
	leftcorns = corners_list$leftcorns
	rightcorns = corners_list$rightcorns
	
	# Plot splits on corners
	# Ensure correct order
	MLsplits = order_LGnodes(MLsplits, tr, removechar=removechar, type=type)
	MLsplits
	
	
	# Get the ranges at the tips (in the right order)
	tipranges = order_tipranges_by_tree_tips(tipranges, tr)
	tiprange_names = tipranges_to_area_strings(tipranges=tipranges)
	
	if (areanames == "default")
		{
		areanames = getareas_from_tipranges_object(tipranges=tipranges)
		}
	
	# Colors 
	if (statesColors_table == "default")
		{
		statesColors_table = get_statesColors_table(areanames)
		tmp_index = get_indices_where_list1_occurs_in_list2(list1=MLsplits$leftBB, list2=statesColors_table$range)
		left_colors = statesColors_table$color[tmp_index]
		tmp_index = get_indices_where_list1_occurs_in_list2(list1=MLsplits$rightBB, list2=statesColors_table$range)
		right_colors = statesColors_table$color[tmp_index]
		
		tmp_index = get_indices_where_list1_occurs_in_list2(list1=tiprange_names, list2=statesColors_table$range)
		tip_colors = statesColors_table$color[tmp_index]
		} else {
		left_colors = bgcol
		right_colors = bgcol
		tip_colors = bgcol
		}
	
	# Plot them
	if (newplot == TRUE)
		{
		plot(tr, label.offset=0.15)
		axisPhylo()
		mtext(text="million years ago", side=1, line=2)
		}
	cornerlabels(text=MLsplits$leftBB, coords=leftcorns, bg=left_colors, ...)
	cornerlabels(text=MLsplits$rightBB, coords=rightcorns, bg=right_colors, ...)
	tiplabels(text=tiprange_names, bg=tip_colors, ...)

	return(MLsplits)
	}



#######################################################
# map_LG_MLstates_to_tree
#######################################################
#' Map states to the nodes on a phylogeny
#' 
#' What it says.
#' 
#' @param MLstates_LGcpp A data.frame containing the node numbers, states, and states probabilities.
#' @param tr An ape phylo object
#' @param tipranges Tipranges object
#' @param removechar The character to remove, if needed.
#' @param type The type of LAGRANGE input (default C++)
#' @param statesColors_table If not default, a table with a color for each area combination.
#' @param bgcol The background color
#' @param areanames The area names, if different from those in the tipranges object
#' @param newplot Default TRUE; should there be a new plot, or should the splits be added to another plot?
#' @param ... Additional arguments to standard functions
#' @return \code{MLstates_LGcpp} The states table, ordered appropriately.
#' @export
#' @seealso \code{\link{get_lagrange_nodenums}}, \code{\link{LGpy_splits_fn_to_table}}, \code{\link{LGcpp_splits_fn_to_table}}
#' @note Go (BioGeo)BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @examples
#' test=1
#' 
map_LG_MLstates_to_tree <- function(MLstates_LGcpp, tr, tipranges, removechar=NULL, type="C++", statesColors_table="default", bgcol="green3", areanames="default", newplot=TRUE, ...)
	{
	# Order in LAGRANGE order
	MLstates_LGcpp = order_LGnodes(MLstates_LGcpp, tr, removechar, type=type, type2="states")
	
	# Get the ranges at the tips (in the right order)
	tipranges = order_tipranges_by_tree_tips(tipranges, tr)
	tiprange_names = tipranges_to_area_strings(tipranges=tipranges)
	
	if (areanames == "default")
		{
		areanames = getareas_from_tipranges_object(tipranges=tipranges)
		}
	
	# Colors 
	if (statesColors_table == "default")
		{
		statesColors_table = get_statesColors_table(areanames)
		tmp_index = get_indices_where_list1_occurs_in_list2(list1=MLstates_LGcpp$states, list2=statesColors_table$range)
		state_colors = statesColors_table$color[tmp_index]
		
		tmp_index = get_indices_where_list1_occurs_in_list2(list1=tiprange_names, list2=statesColors_table$range)
		tip_colors = statesColors_table$color[tmp_index]
		} else {
		left_colors = bgcol
		right_colors = bgcol
		tip_colors = bgcol
		}

	# Plot them
	if (newplot == TRUE)
		{
		plot(tr, label.offset=0.15)
		axisPhylo()
		mtext(text="million years ago", side=1, line=2)
		}

	nodelabels(text=MLstates_LGcpp$states, node=20:37, bg=state_colors, ...)
	tiplabels(tiprange_names, bg=tip_colors, ...)
	#title(paste("LAGRANGE (", type, ") ML states", sep=""))
	
# 	plot(tr, label.offset=0.15)
# 	pievals = as.matrix(MLstates_LGcpp[,c("relprob","relprob2")])
# 	nodelabels(pie=pievals, piecol=c("blue", "white"))
# 	tiplabels(tiprange_names)
# 	title(paste("LAGRANGE (", type, ") state probs", sep=""))
# 	axisPhylo()
# 	mtext(text="million years ago", side=1, line=2)

	return(MLstates_LGcpp)
	}






#######################################################
# order_LGnodes
#######################################################
#' Order LAGRANGE-numbered nodes so that they can be plotted in R
#' 
#' What it says.
#' 
#' @param MLsplits_LGcpp A data.frame containing the node numbers, splits, and split probabilities.
#' @param tr An ape phylo object
#' @param removechar The character to remove, if needed.
#' @param type The type of LAGRANGE input (default C++)
#' @param type2 "splits" or "states"
#' @return \code{MLsplits} The splits table, ordered appropriately.
#' @export
#' @seealso \code{\link{get_lagrange_nodenums}}, \code{\link{LGpy_splits_fn_to_table}}, \code{\link{LGcpp_splits_fn_to_table}}
#' @note Go (BioGeo)BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @examples
#' test=1
#' 
order_LGnodes <- function(MLsplits_LGcpp, tr=NULL, removechar=NULL, type="C++", type2="splits")
	{
	# Order in LAGRANGE order
	if (type == "C++")
		{
		MLsplits_LGcpp = MLsplits_LGcpp[order(MLsplits_LGcpp$nodenum_LGcpp),]
		} else {
		MLsplits_LGcpp = MLsplits_LGcpp[order(MLsplits_LGcpp$nodenum_LGpy),]
		}
	MLsplits_LGcpp
	
	# Get the names of the columns in the MLsplits table
	tmpnames = names(MLsplits_LGcpp)

	# Add R node numbering, if not already done
	if (("Rnodes" %in% tmpnames) == FALSE)
		{
		downpass_node_matrix = get_lagrange_nodenums(tr)
		#downpass_node_matrix[,2] = 1:18
		#downpass_node_matrix = downpass_node_matrix[order(downpass_node_matrix[,1]), ]
	
		MLsplits_LGcpp = cbind(MLsplits_LGcpp, downpass_node_matrix)
		names(MLsplits_LGcpp) = c(tmpnames, "Rnodes", "LGnodes")
		MLsplits_LGcpp
		}
	
	MLsplits_LGcpp = MLsplits_LGcpp[order(MLsplits_LGcpp$Rnodes), ]
	MLsplits_LGcpp
	
	# Change e.g. O_M -> OM, K_O_M_H -> KOMH
	if (!is.null(removechar))	
		{
		if (type2 == "splits")
			{
			MLsplits_LGcpp$splits = gsub(pattern=removechar, replacement="", x=MLsplits_LGcpp$splits)
			MLsplits_LGcpp$leftBB = gsub(pattern=removechar, replacement="", x=MLsplits_LGcpp$leftBB)
			MLsplits_LGcpp$rightBB = gsub(pattern=removechar, replacement="", x=MLsplits_LGcpp$rightBB)
			} else {
			MLsplits_LGcpp$states = gsub(pattern=removechar, replacement="", x=MLsplits_LGcpp$states)
			}
		}
	
	return(MLsplits_LGcpp)
	}



#######################################################
# cornerlabels
#######################################################
#' Make labels for plotting ranges on corners
#' 
#' This function makes labels for plotting ranges on corners.
#' 
#' @param text The text to put at the corners.
#' @param coords The coordinates at which to plot the labels
#' @param bg The background color
#' @param col The text color
#' @param adj Position adjustment; default \code{adj=c(0.5,0.5)}
#' @param ... Additional arguments to standard functions
#' @return nothing
#' @export
#' @seealso \code{\link{cornerpies}}, \code{\link{corner_coords}}, \code{\link{get_lagrange_nodenums}}, 
#' \code{\link{LGpy_splits_fn_to_table}}, \code{\link{LGcpp_splits_fn_to_table}}
#' @note Go (BioGeo)BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @examples
#' test=1
#' 
cornerlabels <- function(text, coords, bg="green3", col="black", adj=c(0.5,0.5), ...)
	{
	args <- list(...)
    CEX <- if ("cex" %in% names(args)) 
    	{
	    args$cex
	    } else {
	    par("cex")
	    }
	   
	# Draw the rectangles	   
	width <- strwidth(text, units = "inches", cex = CEX)
	height <- strheight(text, units = "inches", cex = CEX)

	width <- xinch(width)
	height <- yinch(height)

	XX = coords$x
	YY = coords$y	
	xl <- XX - width * adj[1] - xinch(0.03)
	xr <- xl + width + xinch(0.03)
	yb <- YY - height * adj[2] - yinch(0.02)
	yt <- yb + height + yinch(0.05)
	rect(xl, yb, xr, yt, col = bg)
	
	# Write the text
	text(XX, YY, text, adj = adj, col = col, ...)
	
	return()
	}



#######################################################
# cornerpies
#######################################################
#' Make pie charts for plotting ranges on corners
#' 
#' This function makes pie charts for plotting ranges on corners.  It makes use of 
#' \code{ape:::floating.pie.asp} to plot the pie charts on the corners.
#' 
#' To get the corner coordinates, use \code{\link{corner_coords}}.  Please note the 
#' special input required in that function to get it to access a corner-coordinates 
#' function in the extensions data (\code{extdata}) directory.
#' 
#' @param pievals The matrix (numnodes x numstates) of probabilities to plot.
#' @param coords The coordinates at which to plot the labels.
#' @param piecol The color for each possible state.
#' @param adj Position adjustment; default \code{adj=c(0.5,0.5)}
#' @param ... Additional arguments to standard functions
#' @return nothing
#' @export
#' @seealso \code{\link{cornerlabels}}, \code{\link{corner_coords}}, \code{\link{get_lagrange_nodenums}}, 
#' \code{\link{LGpy_splits_fn_to_table}}, \code{\link{LGcpp_splits_fn_to_table}}
#' @note Go (BioGeo)BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @examples
#' test=1
#' 
cornerpies <- function(pievals, coords, piecol, adj=c(0.5,0.5), ...)
	{
	require(ape)	# for ape:::floating.pie.asp
	
	args <- list(...)
    CEX <- if ("cex" %in% names(args)) 
    	{
	    args$cex
	    } else {
	    par("cex")
	    }
	   
# 	# Draw the rectangles	   
# 	width <- strwidth(text, units = "inches", cex = CEX)
# 	height <- strheight(text, units = "inches", cex = CEX)
# 
# 	width <- xinch(width)
# 	height <- yinch(height)

	XX = coords$x
	YY = coords$y	
# 	xl <- XX - width * adj[1] - xinch(0.03)
# 	xr <- xl + width + xinch(0.03)
# 	yb <- YY - height * adj[2] - yinch(0.02)
# 	yt <- yb + height + yinch(0.05)
# 	rect(xl, yb, xr, yt, col = bg)
# 	
# 	# Write the text
# 	text(XX, YY, text, adj = adj, col = col, ...)
	
	if (is.vector(pie))
		{
		pie <- cbind(pie, 1 - pie)
		}
	xrad <- CEX * diff(par("usr")[1:2])/50
	xrad <- rep(xrad, nrow(pievals))
	XX <- XX + adj[1] - 0.5
	YY <- YY + adj[2] - 0.5
	
	# Loop through the nodes
	for (i in 1:nrow(pievals))
		{
		if (any(is.na(pievals[i, ]))) 
			{
			next()
			}
		ape:::floating.pie.asp(XX[i], YY[i], pievals[i, ], radius = xrad[i], col = piecol, ...)
		}
	
	
	return()
	}





#######################################################
# add_corners
#######################################################
#' Iterate up through a plotted tree, getting the coordinates of the corners 
#' 
#' What it says.
#'
#' @param startnode The node to start at (this is a recursive function)
#' @param tr A tree object in \code{\link[ape]{phylo}} format.
#' @param nodecoords The accumulating list of node coordinates
#' @param corners_list The accumulating list of corners
#' @return \code{corners_list} 
#' @export
#' @seealso \code{\link[ape]{phylo}}, \code{\link{get_nodenums}}
#' @note Go (BioGeo)BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' blah=1
#' 
add_corners <- function(startnode, tr, nodecoords, corners_list)
	{
	# Get the daughters of this node 
	daughtersTF = tr$edge[,1] == startnode
	
	# Run through them if they exist
	if (sum(daughtersTF) > 0)
		{
		# Get the daughter nodenums
		daughters = tr$edge[daughtersTF,2]

		# Get and store node number
		nums = 1:nrow(corners_list$nodevals)
		filled_TF = !is.na(corners_list$nodevals)
		
		if (sum(filled_TF) > 0)
			{
			row_we_are_at = 1 + max(nums[filled_TF])
			} else {
			row_we_are_at = 1
			}
		
		# Store the location
		corners_list$nodevals[row_we_are_at, 1] = startnode
		
		# Get left corner
		# x coordinate
		corners_list$corner_coords_left[row_we_are_at,1] = nodecoords[startnode,1]
		# y coordinate
		corners_list$corner_coords_left[row_we_are_at,2] = nodecoords[daughters[2],2]

		# Get right corner
		# x coordinate
		corners_list$corner_coords_right[row_we_are_at,1] = nodecoords[startnode,1]
		# y coordinate
		corners_list$corner_coords_right[row_we_are_at,2] = nodecoords[daughters[1],2]
		
		# Then propagate up
		for (d in daughters)
			{
			#startnode = d
			corners_list = add_corners(startnode=d, tr, nodecoords, corners_list)
			}
		
		return(corners_list)
		} else {
		return(corners_list)		
		}
	return(corners_list)
	}

# Get the corner coordinates
#######################################################
# corner_coords
#######################################################
#' Get the corner coordinates
#' 
#' Gets the coordinates of the corners when the tree is plotted.
#'
#' Because this function needs to use a modified version of the APE plot.phylo
#' function, and for complex reasons APE's .C functions cannot be used 
#' elsewhere without causing problems with R CMD check, this function is
#' left up to user specification.  Basically, the user puts in
#' the name of the function, which is available in the extension data
#' (\code{extdata/a_scripts}) directory of the package.  The defaults work on the 
#' developer's machine, other users may have to e.g. change "manual" to \code{tmplocation},
#' where \code{tmplocation} is specified as in the example.
#'
#' @param tr A tree object in \code{\link[ape]{phylo}} format.
#' @param coords_fun The name of the function to use to get node coordinates. Default: 
#' "plot_phylo3_nodecoords". 
#' @param tmplocation Default is "manual", which throws an error check unless your path structure matches the developer's.  Most users
#' should probably use the \code{\link[base]{system.file}} command in the examples, below. The directory location containing the R 
#' script \code{plot_phylo3_nodecoords.R}. This function, modified from the \code{\link[ape]{ape}} function
#' \code{\link[ape]{plot.phylo}}, cannot be included directly in the R package as it contains C code that does not pass CRAN's R CMD check. The default, 
#' cornercoords_loc="manual", will not allow split states to be plot.  The R script \code{plot_phylo3_nodecoords.R} is located in the BioGeoBEARS extension data 
#' directory, \code{extdata/a_scripts}.  You should be able to get the full path with 
#' \code{list.files(system.file("extdata/a_scripts", package="BioGeoBEARS"), full.names=TRUE)}.
#' @param root.edge Plot the root.edge, if it exists? Default TRUE.
#' @return \code{corners_list} 
#' @export
#' @seealso \code{\link[ape]{phylo}}, \code{\link{get_nodenums}}
#' @note Go (BioGeo)BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' 
#' # Set location like this if you don't have plot_phylo3_nodecoords
#' # hardcoded/sourced elsehwhere
#' # tmplocation = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
#' # 
#' \dontrun{
#' extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
#' trfn = np(paste(extdata_dir, "/Psychotria_5.2.newick", sep=""))
#' tr = read.tree(trfn)
#' tmplocation = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
#' corner_coords(tr, coords_fun="plot_phylo3_nodecoords", tmplocation=tmplocation)
#' }
#' 
#'
#' 
corner_coords <- function(tr, coords_fun="plot_phylo3_nodecoords", tmplocation="manual", root.edge=TRUE)
	{
	defaults='
	coords_fun="plot_phylo3_nodecoords"
	tmplocation="manual"
	'
	
	# The plot_phylo3_nodecoords causes problems for CRAN, since APE's .C functions
	# aren't exportable, or something.
	# So, we'll source this script from extdata

	# coords_fun="plot_phylo3_nodecoords"; location="manual"

	# trfn = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/
	# inst/extdata/examples/Psychotria_M0/LGcpp/Psychotria_5.2.newick"
	
	# tr = read.tree(trfn)
	
	if (tmplocation != "manual")
		{
		scriptdir = tmplocation
		}

	if (tmplocation == "manual")
		{
		scriptdir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/a_scripts/"
		}
	file_to_source = np(paste(addslash(scriptdir), coords_fun, ".R", sep=""))
	source(file=file_to_source)
	
	# Get the coordinates of the vertices in the tree
	
	# Initialize
	trcoords = matrix(data=c(1,1), nrow=1, ncol=2)
	
	# Set up the command as a string
	# trcoords = plot_phylo3_nodecoords(tr, plot=FALSE)

	# 2017-11-02_edit for new versino of APE
	if (packageVersion("ape") < 5.0)
		{
		# Set up the command as a string
		# trcoords = plot_phylo3_nodecoords(tr, plot=FALSE)
		cmdstr = paste("trcoords = ", coords_fun, "(tr, plot=FALSE, root.edge=root.edge)", sep="")
		eval(parse(text=cmdstr))
	
		# X and Y coords for nodes, 1-37
		nodecoords = cbind(trcoords$xx, trcoords$yy)
		} else {
		# Temporary plot, not to screen (hopefully)
		trcoords = plot_phylo3_nodecoords_APE5(tr, plot=FALSE, root.edge=root.edge)
		
		# X and Y coords for nodes, 1-37
		nodecoords = cbind(trcoords$xx, trcoords$yy)
		} # END if (packageVersion("ape") < 5.0)
	
	# Go through the edge matrix from the root, take the x coord of the node,
	# and the y coord of the descendant
	rootnode = get_nodenum_structural_root(tr)
	
	# Get corner coords
	corners_list = NULL
	corners_list$corner_coords_left = matrix(NA, ncol=2, nrow=tr$Nnode)
	corners_list$corner_coords_right = matrix(NA, ncol=2, nrow=tr$Nnode)
	corners_list$nodevals = matrix(NA, ncol=1, nrow=tr$Nnode)
	
	corners_list = add_corners(startnode=rootnode, tr, nodecoords, corners_list)
	
	
	left = corners_list$corner_coords_left
	right = corners_list$corner_coords_right
	node = corners_list$nodevals
	
	leftcorns = adf(cbind(node, left))
	names(leftcorns) = c("node", "x", "y")
	row.names(leftcorns) = node

	rightcorns = adf(cbind(node, right))
	names(rightcorns) = c("node", "x", "y")
	row.names(rightcorns) = node
	
	corners_list = NULL
	corners_list$leftcorns = leftcorns
	corners_list$rightcorns = rightcorns
	
	return(corners_list)
	}



# Get the node coordinates
#######################################################
# node_coords
#######################################################
#' Get the node coordinates
#' 
#' Gets the coordinates of the nodes when the tree is plotted.
#'
#' Because this function needs to use a modified version of the APE plot.phylo
#' function, and for complex reasons APE's .C functions cannot be used 
#' elsewhere without causing problems with R CMD check, this function is
#' left up to user specification.  Basically, the user puts in
#' the name of the function, which is available in the extension data
#' (\code{extdata/a_scripts}) directory of the package.  The defaults work on the 
#' developer's machine, other users may have to e.g. change "manual" to \code{tmplocation},
#' where \code{tmplocation} is specified as in the example.
#'
#' @param tr A tree object in \code{\link[ape]{phylo}} format.
#' @param coords_fun The name of the function to use to get node coordinates. Default: 
#' "plot_phylo3_nodecoords". 
#' @param tmplocation Default is "manual", which throws an error check unless your path structure matches the developer's.  Most users
#' should probably use the \code{\link[base]{system.file}} command in the examples, below. The directory location containing the R 
#' script \code{plot_phylo3_nodecoords.R}. This function, modified from the \code{\link[ape]{ape}} function
#' \code{\link[ape]{plot.phylo}}, cannot be included directly in the R package as it contains C code that does not pass CRAN's R CMD check. The default, 
#' cornercoords_loc="manual", will not allow split states to be plot.  The R script \code{plot_phylo3_nodecoords.R} is located in the BioGeoBEARS extension data 
#' directory, \code{extdata/a_scripts}.  You should be able to get the full path with 
#' \code{list.files(system.file("extdata/a_scripts", package="BioGeoBEARS"), full.names=TRUE)}.
#' @param root.edge Plot the root.edge, if it exists? Default TRUE.
#' @return \code{nodecoords}, a data.frame with nodecoords$x and nodecoords$y specifying
#' the plot coordinates of each node.
#' @export
#' @seealso \code{\link[ape]{phylo}}, \code{\link{get_nodenums}}
#' @note Go (BioGeo)BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' 
#' # Set location like this if you don't have plot_phylo3_nodecoords
#' # hardcoded/sourced elsehwhere
#' # tmplocation = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
#' # 
#' \dontrun{
#' extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
#' trfn = np(paste(extdata_dir, "/Psychotria_5.2.newick", sep=""))
#' tr = read.tree(trfn)
#' tmplocation = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
#' nodes_xy = node_coords(tr, coords_fun="plot_phylo3_nodecoords", tmplocation=tmplocation)
#' nodes_xy
#' }
#' 
#'
#' 
node_coords <- function(tr, coords_fun="plot_phylo3_nodecoords", tmplocation="manual", root.edge=TRUE)
	{
	defaults='
	coords_fun="plot_phylo3_nodecoords"
	tmplocation="manual"
	'
	
	# The plot_phylo3_nodecoords causes problems for CRAN, since APE's .C functions
	# aren't exportable, or something.
	# So, we'll source this script from extdata

	# coords_fun="plot_phylo3_nodecoords"; location="manual"

	# trfn = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/
	# inst/extdata/examples/Psychotria_M0/LGcpp/Psychotria_5.2.newick"
	
	# tr = read.tree(trfn)
	
	if (tmplocation != "manual")
		{
		scriptdir = tmplocation
		}

	if (tmplocation == "manual")
		{
		scriptdir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/a_scripts/"
		}
	file_to_source = np(paste(addslash(scriptdir), coords_fun, ".R", sep=""))
	source(file=file_to_source)
	
	# Get the coordinates of the vertices in the tree
	
	# Initialize
	trcoords = matrix(data=c(1,1), nrow=1, ncol=2)
	
	# Set up the command as a string
	# trcoords = plot_phylo3_nodecoords(tr, plot=FALSE)

	# 2017-11-02_edit for new versino of APE
	if (packageVersion("ape") < 5.0)
		{
		# Set up the command as a string
		# trcoords = plot_phylo3_nodecoords(tr, plot=FALSE)
		cmdstr = paste("trcoords = ", coords_fun, "(tr, plot=FALSE, root.edge=root.edge)", sep="")
		eval(parse(text=cmdstr))
	
		# X and Y coords for nodes, 1-37
		nodecoords = cbind(trcoords$xx, trcoords$yy)
		} else {
		# Temporary plot, not to screen (hopefully)
		trcoords = plot_phylo3_nodecoords_APE5(tr, plot=FALSE, root.edge=root.edge)
		
		# X and Y coords for nodes, 1-37
		nodecoords = cbind(trcoords$xx, trcoords$yy)
		} # END if (packageVersion("ape") < 5.0)


	nodes_xy = adf(nodecoords)
	names(nodes_xy) = c("x", "y")
	row.names(nodes_xy) = 1:nrow(nodecoords)
	
	return(nodes_xy)
	}






#######################################################
# PLOT PROBABILITIES ON A PER-AREA BASIS
#######################################################
#' @param border The color of the border of the boxes holding areas. Default is
#' the string "default", which means that the borders will be "gray50" if 
#' there is no color specified (that is, cols_each_area=NULL, which means bars 
#' will be "gray70"), and borders will be "black" if color is 
#' specified.  By changing the box borders and the tree color (see 
#' parameter \code{plot_per_area_probs} in \code{plot_per_area_probs}, or \code{edge.color} in \code{ape::plot.phylo}).
#' the user can emphasize or de-emphasize one or the other.
#' @param trcol The color of the lines in the phylogeny when plotted.  
#' Default is "black", but "gray60" looks good if you want to 
#' de-emphasize the tree with respect to e.g. node areas.

plot_per_area_probs <- function(tr, res, areas, states_list_0based, titletxt="", cols_each_area=NULL, barwidth_proportion=0.02, barheight_proportion=0.025, offset_nodenums=NULL, offset_xvals=NULL, offset_yvals=NULL, root.edge=TRUE, border="default", trcol="black", plot_rangesizes=FALSE)
	{
	defaults='
	areas = getareas_from_tipranges_object(tipranges)
	states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
	states_list_0based

	cols_each_area=NULL
	offset_nodenums=NULL
	offset_xvals=NULL
	offset_yvals=NULL
	barwidth_proportion=0.02
	barheight_proportion=0.025
	
	border="default"
	trcol="black"
	plot_rangesizes=FALSE
	' # END defaults
	
	
	# Error check
	if ((length(cols_each_area) == 1) && (cols_each_area == FALSE))
		{
		errortxt = "\n\nERROR IN plot_per_area_probs():\ncols_each_area=FALSE, but it should be either:\nNULL (which gives grey boxes)\nTRUE (which makes colors auto-generated), or \na list of colors the same length as 'areas'.\n\n"
		cat(errortxt)
		
		stop("\n\nStopping on error.\n\n")
		
		} # END if ((length(cols_each_area) == 1) && (cols_each_area == FALSE))
	
	
	# Get the relative probabilities of each state/range
	relprobs_matrix = res$ML_marginal_prob_each_state_at_branch_top_AT_node
	
	# Collapse to the probabilities of each area
	if (plot_rangesizes == FALSE)
		{
		probs_each_area = infprobs_to_probs_of_each_area(relprobs_matrix, states_list=states_list_0based)
		}
	
	# Collapse to the probabilities of each range size
	if (plot_rangesizes == TRUE)
		{
		rangesizes_by_state = sapply(FUN=length, X=states_list_0based)
		rangesizes = sort(unique(rangesizes_by_state))
		rangesizes
		
		probs_each_area = infprobs_to_probs_of_each_rangesize(relprobs_matrix, states_list=states_list_0based)
		}
	dim(probs_each_area)


	# To get offset_tiplabels:
	ntips = length(tr$tip.label)
	numnodes = tr$Nnode
	tipnums = 1:ntips
	nodenums = (ntips+1):(ntips+numnodes)

	extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
	tmplocation=paste(extdata_dir, "a_scripts" , sep="/")
	nodes_xy = node_coords(tr, root.edge=root.edge, tmplocation=tmplocation)
	nodes_xy

	# Make bar width a proportion of the width of the plot in x
	#barwidth_proportion = 0.02
	barwidth = (max(nodes_xy$x) - min(nodes_xy$x)) * barwidth_proportion
	barwidth

	#barheight_proportion = 0.025
	barheight = (max(nodes_xy$y) - min(nodes_xy$y)) * barheight_proportion
	barheight

	numareas = ncol(probs_each_area)
	areanums = 1:numareas
	middle = median(areanums)
	middle
	offsets_nodes = (areanums - middle)*barwidth
	offsets_tips = (areanums - 0)*barwidth
	offset_tiplabels = max(offsets_tips) + barwidth/1



	# Plot the tree
	plot(tr, label.offset=offset_tiplabels, root.edge=root.edge, edge.color=trcol)
	# plot(tr, label.offset=offset_tiplabels)
	axisPhylo()
	title(titletxt)

	# Add the areas boxes
	add_per_area_probs_to_nodes(tr, probs_each_area, cols_each_area=cols_each_area, barwidth_proportion=barwidth_proportion, barheight_proportion=barheight_proportion, offset_nodenums=offset_nodenums, offset_xvals=offset_xvals, offset_yvals=offset_yvals, border=border)
	
	return(probs_each_area)
	}





#######################################################
# add_per_area_probs_to_nodes
#######################################################

#' @param border The color of the border of the boxes holding areas. Default is
#' the string "default", which means that the borders will be "gray50" if 
#' there is no color specified (that is, cols_each_area=NULL, which means bars 
#' will be "gray70"), and borders will be "black" if color is 
#' specified.  By changing the box borders and the tree color (see 
#' parameter \code{plot_per_area_probs} in \code{plot_per_area_probs}, 
#' or \code{edge.color} in \code{ape::plot.phylo}).

add_per_area_probs_to_nodes <- function(tr, probs_each_area, cols_each_area=NULL, barwidth_proportion=0.02, barheight_proportion=0.025, offset_nodenums=NULL, offset_xvals=NULL, offset_yvals=NULL, root.edge=TRUE, border="default")
	{
	defaults='
	cols_each_area=NULL
	offset_nodenums=NULL
	offset_xvals=NULL
	offset_yvals=NULL
	barwidth_proportion=0.02
	barheight_proportion=0.025
	' # END defaults


	# Error check
	if ((length(cols_each_area) == 1) && (cols_each_area == FALSE))
		{
		errortxt = "\n\nERROR IN add_per_area_probs_to_nodes():\n\ncols_each_area=FALSE, but it should be either:\nNULL (which gives grey boxes)\nTRUE (which makes colors auto-generated), or \na list of colors the same length as 'areas'.\n\n"
		cat(errortxt)
		
		stop("\n\nStopping on error.\n\n")
		
		} # END if ((length(cols_each_area) == 1) && (cols_each_area == FALSE))
	

	ntips = length(tr$tip.label)
	numnodes = tr$Nnode
	tipnums = 1:ntips
	nodenums = (ntips+1):(ntips+numnodes)

	extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
	tmplocation=paste(extdata_dir, "a_scripts" , sep="/")
	nodes_xy = node_coords(tr, root.edge=root.edge, tmplocation=tmplocation)
	nodes_xy

	# Make bar width a proportion of the width of the plot in x
	#barwidth_proportion = 0.02
	barwidth = (max(nodes_xy$x) - min(nodes_xy$x)) * barwidth_proportion
	barwidth

	#barheight_proportion = 0.025
	barheight = (max(nodes_xy$y) - min(nodes_xy$y)) * barheight_proportion
	barheight

	numareas = ncol(probs_each_area)
	areanums = 1:numareas
	middle = median(areanums)
	middle
	offsets_nodes = (areanums - middle)*barwidth
	offsets_tips = (areanums - 0)*barwidth
	offset_tiplabels = max(offsets_tips) + barwidth/1

	# Draw the empty boxes

	# xcoords for tips
	xlefts_tips = sapply(X=nodes_xy$x[tipnums], FUN="+", (offsets_tips - barwidth/2))
	xlefts_nodes = sapply(X=nodes_xy$x[nodenums], FUN="+", (offsets_nodes - barwidth/2))
	xrights_tips = sapply(X=nodes_xy$x[tipnums], FUN="+", (offsets_tips + barwidth/2))
	xrights_nodes = sapply(X=nodes_xy$x[nodenums], FUN="+", (offsets_nodes + barwidth/2))
	
	xlefts_matrix = t(cbind(xlefts_tips, xlefts_nodes))
	xrights_matrix = t(cbind(xrights_tips, xrights_nodes))

	ybottoms_per_node = sapply(X=nodes_xy$y, FUN="-", (barheight/2))
	ytops_per_node = sapply(X=nodes_xy$y, FUN="+", (barheight/2))
	
	# Manually modify some positions
	if (is.null(offset_nodenums) == FALSE)
		{
		xlefts_matrix[offset_nodenums,] = xlefts_matrix[offset_nodenums,] + (offset_xvals*barwidth)
		xrights_matrix[offset_nodenums,] = xrights_matrix[offset_nodenums,] + (offset_xvals*barwidth)
		ybottoms_per_node[offset_nodenums] = ybottoms_per_node[offset_nodenums] + (offset_yvals*barheight)
		ytops_per_node[offset_nodenums] = ytops_per_node[offset_nodenums] + (offset_yvals*barheight)
		}
	
	# Convert these into plain lists
	xlefts = c(xlefts_matrix)
	xrights = c(xrights_matrix)

	ybottoms = rep(ybottoms_per_node, times=numareas)
	ytops = rep(ytops_per_node, times=numareas)


	# Plot the box outlines
	#nulls = mapply(FUN=rect, xleft=xlefts, ybottom=ybottoms, xright=xrights, ytop=ytops, MoreArgs=list(col="white", border=border))

	# Then draw black boxes inside these
	# You just have to adjust the top of the black bar, based on prob
	yranges_per_node = ytops_per_node - ybottoms_per_node
	yadd_above_ybottom = yranges_per_node * probs_each_area
	ytops_probs_node = yadd_above_ybottom + ybottoms
	ytops_probs = c(ytops_probs_node)
	

	# IF cols_each_area == TRUE, auto-generate colors
	if (length(cols_each_area) == 1 && (cols_each_area == TRUE))
		{
		tmp_colors_matrix = get_colors_for_numareas(numareas, use_rainbow=FALSE)
		#cols_each_area = c("blue", "green", "yellow", "red")
		cols_each_area = mapply(FUN=rgb, red=tmp_colors_matrix[1,], green=tmp_colors_matrix[2,], blue=tmp_colors_matrix[3,], MoreArgs=list(maxColorValue=255))

		}
	

	
	# Default color is darkgray
	if (is.null(cols_each_area))
		{
		# Auto-generate border, also -- gray50 looks black against lighter gray70
		if (border == "default")
			{
			border = "gray50"
			}

		# Plot the box outlines
		nulls = mapply(FUN=rect, xleft=xlefts, ybottom=ybottoms, xright=xrights, ytop=ytops, MoreArgs=list(col="white", border=border))

		# Draw bars
		nulls = mapply(FUN=rect, xleft=xlefts, ybottom=ybottoms, xright=xrights, ytop=ytops_probs, MoreArgs=list(col="gray70", border=border))
		} else {
		
		
		if ( length(cols_each_area) != length(areas))
			{
			errortxt = paste("\n\nERROR in add_per_area_probs_to_nodes():\n\nif cols_each_area is specified, length(cols_each_area) must equal length(areas).\n\n", sep="")
			
			cat(errortxt)
			
			cat("Your 'areas':\n\n", sep="")
			print(areas)
			
			cat("Your 'cols_each_area':\n\n", sep="")
			print(cols_each_area)
			
			stop("Stopping on error.")
			} # END if ( length(cols_each_area) != length(areas))
		
		# Otherwise, expand colors to each box
		cols_each_area_matrix = matrix(data=cols_each_area, nrow=(ntips+numnodes), ncol=length(cols_each_area), byrow=TRUE)
		
		
		# Plot with different internal box colors
		cols_each_area = c(cols_each_area_matrix)

		# Auto-generate border with colored boxes (black looks best), also
		if (border == "default")
			{
			border = "black"
			}
		
		# Plot the box outlines
		nulls = mapply(FUN=rect, xleft=xlefts, ybottom=ybottoms, xright=xrights, ytop=ytops, MoreArgs=list(col="white", border=border))
		
		# Plot the colored bars
		nulls = mapply(FUN=rect, xleft=xlefts, ybottom=ybottoms, xright=xrights, ytop=ytops_probs, col=cols_each_area, MoreArgs=list(border=border))
		} # END if (is.null(cols_each_area))

	return(NULL)
	}





#######################################################
# add_per_area_probs_to_corners
#######################################################

#' @param border The color of the border of the boxes holding areas. Default is
#' the string "default", which means that the borders will be "gray50" if 
#' there is no color specified (that is, cols_each_area=NULL, which means bars 
#' will be "gray70"), and borders will be "black" if color is 
#' specified.  By changing the box borders and the tree color (see 
#' parameter \code{plot_per_area_probs} in \code{plot_per_area_probs}, 
#' or \code{edge.color} in \code{ape::plot.phylo}).

add_per_area_probs_to_corners <- function(tr, probs_each_area, left_or_right, cols_each_area=NULL, barwidth_proportion=0.02, barheight_proportion=0.025, offset_nodenums=NULL, offset_xvals=NULL, offset_yvals=NULL, root.edge=TRUE, border="default", trcol)
	{
	defaults='
	cols_each_area=NULL
	offset_nodenums=NULL
	offset_xvals=NULL
	offset_yvals=NULL
	barwidth_proportion=0.02
	barheight_proportion=0.025
	' # END defaults
	
	# Error check
	if ((length(cols_each_area) == 1) && (cols_each_area == FALSE))
		{
		errortxt = "\n\nERROR IN add_per_area_probs_to_corners():\n\ncols_each_area=FALSE, but it should be either:\nNULL (which gives grey boxes)\nTRUE (which makes colors auto-generated), or \na list of colors the same length as 'areas'.\n\n"
		cat(errortxt)
		
		stop("\n\nStopping on error.\n\n")
		
		} # END if ((length(cols_each_area) == 1) && (cols_each_area == FALSE))

	
	ntips = length(tr$tip.label)
	numnodes = tr$Nnode
	tipnums = 1:ntips
	#nodenums = (ntips+1):(ntips+numnodes)
	nodenums = (ntips+1):(ntips+numnodes)
	
	# Get the plot coordinates of the corners ABOVE each node
	corners_list = corner_coords(tr, root.edge=root.edge)
	corners_list
	
	# Get the node numbers of the nodes ABOVE each corner above each node
	leftright_nodes_matrix = get_leftright_nodes_matrix_from_results(tr)
	leftright_nodes_matrix
	
	if (left_or_right == "left")
		{
		corners_xy = corners_list$leftcorns
		nodenums_above_LorR_corner = leftright_nodes_matrix$left
		}
		
	if (left_or_right == "right")
		{
		corners_xy = corners_list$rightcorns
		nodenums_above_LorR_corner = leftright_nodes_matrix$right
		}
	
	
	
	
	# LEFT OR RIGHT SPLITS
	# Probs of each area BELOW the node (nodenum = rownum)
	
	# Make bar width a proportion of the width of the plot in x
	#barwidth_proportion = 0.02
	barwidth = (max(corners_xy$x) - min(corners_xy$x)) * barwidth_proportion
	barwidth

	#barheight_proportion = 0.025
	barheight = (max(corners_xy$y) - min(corners_xy$y)) * barheight_proportion
	barheight

	numareas = ncol(probs_each_area)
	areanums = 1:numareas
	middle = median(areanums)
	middle
	offsets_nodes = (areanums - middle)*barwidth
	offsets_tips = (areanums - 0)*barwidth
	offset_tiplabels = max(offsets_tips) + barwidth/1

	# Draw the empty boxes

	# xcoords for tips
	nodenums_above_LorR_corner
	
	# We just need to plot at the corners above internal nodes, not tips
	#xlefts_tips = sapply(X=corners_xy$x[nodenums], FUN="+", (offsets_tips - barwidth/2))
	rownums = nodenums - ntips
	xlefts_nodes = sapply(X=corners_xy$x[rownums], FUN="+", (offsets_nodes - barwidth/2))
	#xrights_tips = sapply(X=corners_xy$x[tipnums], FUN="+", (offsets_tips + barwidth/2))
	xrights_nodes = sapply(X=corners_xy$x[rownums], FUN="+", (offsets_nodes + barwidth/2))
	
	xlefts_matrix = t(xlefts_nodes)
	xrights_matrix = t(xrights_nodes)

	ybottoms_per_node = sapply(X=corners_xy$y, FUN="-", (barheight/2))
	ytops_per_node = sapply(X=corners_xy$y, FUN="+", (barheight/2))
	
	# Manually modify some positions
	if (is.null(offset_nodenums) == FALSE)
		{
		rownums = match(x=offset_nodenums, table=nodenums)
		xlefts_matrix[rownums,] = xlefts_matrix[rownums,] + (offset_xvals*barwidth)
		xrights_matrix[rownums,] = xrights_matrix[rownums,] + (offset_xvals*barwidth)
		ybottoms_per_node[rownums] = ybottoms_per_node[rownums] + (offset_yvals*barheight)
		ytops_per_node[rownums] = ytops_per_node[rownums] + (offset_yvals*barheight)
		}
	
	# Convert these into plain lists
	xlefts = c(xlefts_matrix)
	xrights = c(xrights_matrix)

	ybottoms = rep(ybottoms_per_node, times=numareas)
	ytops = rep(ytops_per_node, times=numareas)


	# Plot the box outlines
	#nulls = mapply(FUN=rect, xleft=xlefts, ybottom=ybottoms, xright=xrights, ytop=ytops, MoreArgs=list(col="white", border=border))

	# Then draw black boxes inside these
	# You just have to adjust the top of the black bar, based on prob
	yranges_per_node = ytops_per_node - ybottoms_per_node
	yadd_above_ybottom = yranges_per_node * probs_each_area[nodenums_above_LorR_corner,]
	ytops_probs_node = yadd_above_ybottom + ybottoms
	ytops_probs = c(ytops_probs_node)
	

	# IF cols_each_area == TRUE, auto-generate colors
	if (length(cols_each_area) == 1 && (cols_each_area == TRUE))
		{
		tmp_colors_matrix = get_colors_for_numareas(numareas, use_rainbow=FALSE)
		#cols_each_area = c("blue", "green", "yellow", "red")
		cols_each_area = mapply(FUN=rgb, red=tmp_colors_matrix[1,], green=tmp_colors_matrix[2,], blue=tmp_colors_matrix[3,], MoreArgs=list(maxColorValue=255))
		}
	
	
	# Default color is darkgray
	if (is.null(cols_each_area))
		{
		# Auto-generate border, also -- gray50 looks black against lighter gray70
		if (border == "default")
			{
			border = "gray50"
			}
		
		# Plot the box outlines
		nulls = mapply(FUN=rect, xleft=xlefts, ybottom=ybottoms, xright=xrights, ytop=ytops, MoreArgs=list(col="white", border=border))
		
		# Draw bars
		nulls = mapply(FUN=rect, xleft=xlefts, ybottom=ybottoms, xright=xrights, ytop=ytops_probs, MoreArgs=list(col="gray70", border=border))
		} else {
		
		
		if ( length(cols_each_area) != length(areas))
			{
			errortxt = paste("\n\nERROR in add_per_area_probs_to_corners():\n\nif cols_each_area is specified, length(cols_each_area) must equal length(areas).\n\n", sep="")
			cat(errortxt)
			
			cat("Your 'areas':\n\n", sep="")
			print(areas)
			
			cat("Your 'cols_each_area':\n\n", sep="")
			print(cols_each_area)
			
			stop("Stopping on error.")
			} # END if ( length(cols_each_area) != length(area))
		
		# Otherwise, expand colors to each box
		cols_each_area_matrix = matrix(data=cols_each_area, nrow=numnodes, ncol=length(cols_each_area), byrow=TRUE)
		
		
		# Plot with different internal box colors
		cols_each_area = c(cols_each_area_matrix)

		# Auto-generate border with colored boxes (black looks best), also
		if (border == "default")
			{
			border = "black"
			}

		# Plot the box outlines
		nulls = mapply(FUN=rect, xleft=xlefts, ybottom=ybottoms, xright=xrights, ytop=ytops, MoreArgs=list(col="white", border=border))
		
		nulls = mapply(FUN=rect, xleft=xlefts, ybottom=ybottoms, xright=xrights, ytop=ytops_probs, col=cols_each_area, MoreArgs=list(border=border))
		} # END if (is.null(cols_each_area))

	return(NULL)
	}










#######################################################
# plot_BioGeoBEARS_model
#######################################################
#' Graphical display of your anagenetic and cladogenetic biogeography models
#' 
#' This function produces a graphical summary of the model stored in a \code{BioGeoBEARS_run_object}.
#' This could be either an input model, or the result of the ML parameter search.
#'
#' Understanding of phylogenetic methoods in historical biogeography methods is hampered by the
#' difficulty of displaying the models the computer is using.  This function is one attempt to 
#' improve the situation, by plotting the relative weights of the various parameters.
#'
#' @param obj The input object, either a \code{BioGeoBEARS_run_object} (if so, set 
#' \code{obj_is_run_or_results="run"} or an output object from \code{\link{bears_optim_run}} 
#' (if so, specify \code{obj_is_run_or_results="results"}.
#' @param obj_is_run_or_results Specify \code{"run"} or \code{"results"}, as described above 
#' for parameter \code{obj}.
#' @param plotwhat Default is \code{"init"}, which means plotting the starting model parameters. 
#' \code{"est"} plots the estimated model parameters.
#' @param titletxt Additional text for the title of the plot
#' @param statenames State names to pass to \code{\link{plot_cladogenesis_size_probabilities}}. 
#' If \code{NULL} (default), these are auto-generated assuming all areas up to the maximum number 
#' are allowed.
#' @return nada
#' @export
#' @seealso \code{\link{plot_cladogenesis_size_probabilities}}, \code{\link{define_BioGeoBEARS_run}}, \code{\link{define_BioGeoBEARS_model_object}}
#' @note Go (BioGeo)BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' blah=1
#' 
plot_BioGeoBEARS_model <- function(obj, obj_is_run_or_results=NULL, plotwhat="init", titletxt="", statenames=NULL)
	{
	runjunk='
	obj = BioGeoBEARS_run_object
	obj_is_run_or_results = "run"
	plotwhat="init"
	titletxt=""
	'
	
	if (is.null(obj_is_run_or_results) == TRUE)
		{
		stoptxt = "\n\nWith plot_BioGeoBEARS_model(), you must specify whether obj_is_run_or_results='run' or 'results'.\n"
		cat(stoptxt)
		stop(stoptxt)
		}
		
	if (obj_is_run_or_results == "run")
		{
		BioGeoBEARS_run_object = obj
		BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
		} else {
		if (obj_is_run_or_results == "results")
			{
			# Extract the correct inputs and outputs
			BioGeoBEARS_run_object = obj$inputs
			BioGeoBEARS_run_object$BioGeoBEARS_model_object = obj$outputs
			BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
			}		
		else {
			stoptxt = "\n\nWith plot_BioGeoBEARS_model(), you must specify whether obj_is_run_or_results='run' or 'results'.\n"
			cat(stoptxt)
			stop(stoptxt)
			}
		}
	
	# Get areas/states/tips
	# Get geographic ranges at tips
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=np(BioGeoBEARS_run_object$geogfn))
	
	# Get the list of geographic areas
	areas = getareas_from_tipranges_object(tipranges)
	areas_list = seq(0, length(areas)-1, 1)		# 0-base indexes

	# Change the names to tipranges@df:
	# this doesn't make sense if areas_list is 0-based indexes
	names(tipranges@df) = areas_list
	
	if (is.na(BioGeoBEARS_run_object$max_range_size))
		{
		if (is.null(BioGeoBEARS_run_object$states_list))
			{
			# Maximum range size is all areas
			max_range_size = length(areas)
			} else {
			# If not NA
			# Get max rangesize from states list
			max_range_size = max(sapply(X=BioGeoBEARS_run_object$states_list, FUN=length), na.rm=TRUE)
			}
		} else {
		# Maximum range size hard-coded
		max_range_size = BioGeoBEARS_run_object$max_range_size
		}
	max_numareas = max_range_size

	#######################################################
	# Check that no tips have larger ranges than you allowed
	#######################################################
	tipranges_df_tmp = tipranges@df
	tipranges_df_tmp[tipranges_df_tmp=="?"] = 0
	TF = (rowSums(dfnums_to_numeric(tipranges_df_tmp))) > max_range_size
	if (sum(TF, na.rm=TRUE) > 0)
		{
		cat("\n\nERROR: Tips with ranges too big:\n", sep="")
		print(dfnums_to_numeric(tipranges@df)[TF, ])
		cat("\n\nCheck your input geography file!\n", sep="")
		txt = paste("ERROR: Some tips (listed above) have range sizes larger than ", max_range_size, sep="")
		stop(txt)
		}
	
	# Take the list of areas, and get list of possible states
	# (the user can manually input states if they like)
	if (is.null(BioGeoBEARS_run_object$states_list))
		{
		states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=BioGeoBEARS_run_object$include_null_range)
		states_list
		#BioGeoBEARS_run_object$states_list = states_list
		#inputs$states_list = states_list
		} else {
		states_list = BioGeoBEARS_run_object$states_list
		#BioGeoBEARS_run_object$states_list = states_list
		#inputs$states_list = states_list
		}





	
	
	# Get everything
	#BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
	#BioGeoBEARS_model_object = params_into_BioGeoBEARS_model_object(BioGeoBEARS_model_object=BioGeoBEARS_model_object, params=params)
	
	# Update linked parameters
	BioGeoBEARS_model_object = calc_linked_params_BioGeoBEARS_model_object(BioGeoBEARS_model_object)

	# Set the dispersal and extinction rate
	d = BioGeoBEARS_model_object@params_table["d",plotwhat]
	e = BioGeoBEARS_model_object@params_table["e",plotwhat]
	a = BioGeoBEARS_model_object@params_table["a",plotwhat]


	#######################################################
	#######################################################
	# Do branch-length exponentiation if desired
	#######################################################
	#######################################################
	b = BioGeoBEARS_model_object@params_table["b",plotwhat]
	# Modify the edge.lengths
	#phy$edge.length = phy$edge.length ^ b


	#######################################################
	#######################################################
	# Do distance-dependence and dispersal multipliers matrix
	#######################################################
	#######################################################
	# Equal dispersal in all directions (unconstrained)
	# Equal extinction probability for all areas
	areas = areas_list
	
	# If there is a distance matrix, use the first one 
	# (non-stratified analysis, here)
	if ( (is.null(BioGeoBEARS_run_object$list_of_distances_mats) == FALSE))
		{
		distances_mat = BioGeoBEARS_run_object$list_of_distances_mats[[1]]
		} else {
		# Default is all areas effectively equidistant
		distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))
		}

	# Get the exponent on distance, apply to distances matrix
	x = BioGeoBEARS_model_object@params_table["x",plotwhat]
	dispersal_multipliers_matrix = distances_mat ^ x

	# Environmental distances
	if ( (is.null(BioGeoBEARS_run_object$list_of_envdistances_mats) == FALSE))
		{
		envdistances_mat = BioGeoBEARS_run_object$list_of_envdistances_mats[[1]]
		} else {
		# Default is all areas effectively equidistant
		envdistances_mat = matrix(1, nrow=length(areas), ncol=length(areas))
		}

	# Get the exponent on environmental distance, apply to distances matrix
	n = BioGeoBEARS_model_object@params_table["n","est"]
	dispersal_multipliers_matrix = dispersal_multipliers_matrix * envdistances_mat ^ n


	# Apply manual dispersal multipliers, if any
	# If there is a manual dispersal multipliers matrix, use the first one 
	# (non-stratified analysis, here)
	if ( (is.null(BioGeoBEARS_run_object$list_of_dispersal_multipliers_mats) == FALSE))
		{
		manual_dispersal_multipliers_matrix = BioGeoBEARS_run_object$list_of_dispersal_multipliers_mats[[1]]
		} else {
		# Default is all areas effectively equidistant
		manual_dispersal_multipliers_matrix = matrix(1, nrow=length(areas), ncol=length(areas))
		}
	
	# Get the exponent on manual dispersal multipliers
	w = BioGeoBEARS_model_object@params_table["w","est"]

	# Apply element-wise
	dispersal_multipliers_matrix = dispersal_multipliers_matrix * manual_dispersal_multipliers_matrix ^ w

	#######################################################
	# multiply parameter d by dispersal_multipliers_matrix
	#######################################################
	dmat_times_d = dispersal_multipliers_matrix * matrix(d, nrow=length(areas), ncol=length(areas))
	amat = dispersal_multipliers_matrix * matrix(a, nrow=length(areas), ncol=length(areas))




	#######################################################
	#######################################################
	# Do area-dependence and extinction multipliers list
	#######################################################
	#######################################################
	if ( (is.null(BioGeoBEARS_run_object$list_of_area_of_areas) == FALSE))
		{
		area_of_areas = BioGeoBEARS_run_object$list_of_area_of_areas[[1]]
		} else {
		# Default is all areas effectively equidistant
		area_of_areas = rep(1, length(areas))
		}
		
	# Get the exponent on extinction, apply to extinction modifiers	
	u = BioGeoBEARS_model_object@params_table["u",plotwhat]
	extinction_modifier_list = area_of_areas ^ (1 * u)
	
	# Apply to extinction rate
	elist = extinction_modifier_list * rep(e, length(areas))



	
	#######################################################
	# Start a big plot with layout()
	#######################################################
	
	#require(plotrix)
	#par(mfrow=c(1,4))
	#Pdecsize_given_ancsize = cbind(maxent01y, maxent01s, maxent01v, maxent01j)
	#color2D.matplot(x=Pdecsize_given_ancsize, c(1,0), c(1,0), c(1,0), show.values=TRUE, show.legend=TRUE, xlab="Ancestral range size", ylab="P(range size) in smaller descendant", axes=FALSE, nslices=100)

	# Layout so plot has 4 columns, 2 main rows and header/footer, and 2 side columns
	numcols = 6
	plotgroup_header = matrix(data=1, nrow=1, ncol=numcols)
	plotgroup_row2 = matrix(data=c(2,3,4,5,6,7), nrow=1, ncol=numcols)

	spmodel_header = matrix(data=8, nrow=1, ncol=numcols)
	plotgroup_row3 = matrix(data=c(9,10,11,12,13,14), nrow=1, ncol=numcols)
	#plotgroup_footer = matrix(data=15, nrow=1, ncol=numcols)
	plotgroup_footer = matrix(data=c(15,15,16,16,16,16), nrow=1, ncol=numcols)
	plotgroups = rbind(plotgroup_header, plotgroup_row2, spmodel_header, plotgroup_row3, plotgroup_footer)
	plotgroups
	
	layout(mat=plotgroups, widths=c(0.2,1,1,1,1,0.3), heights=c(0.2,1,0.1,1,1))
	
	
	#######################################################
	# Header (cell #1)
	#######################################################
	# set the inside box ("i") ahead of time!!  -- removes the 4% extension
	# this friggin' solution took an hour to find,
	# solution here: http://tolstoy.newcastle.edu.au/R/help/06/08/32529.html
	par(mar=c(0,0,0,0), xaxs = "i", yaxs = "i") 
	plot(x=c(0,1), y=c(0,1), pch=".", col="white", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
	tmptxt = paste("BioGeoBEARS model: ", titletxt, sep="")
	text(x=0.5, y=0.5, tmptxt, cex=1.5, cex.main=2)


	#################################
	# Row 2
	#################################
	# Cell #2 (leftmost, thin)
	plot(0,0, pch=".", col="white", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")

	# cell #3:
	# Table of free vs. fixed params
	plot(x=c(0,1),y=c(0,1), pch=".", col="white", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
	val_txt = c("d", "e", "a", "x", "n", "w", "u", "b")
	Param = val_txt
	Type = BioGeoBEARS_model_object@params_table[val_txt, "type"]
	Init = BioGeoBEARS_model_object@params_table[val_txt, "init"]
	Est = BioGeoBEARS_model_object@params_table[val_txt, "est"]
	dtf = adf2(cbind(Param, Type, Init, Est))
	dtf$Init = round(as.numeric(dtf$Init), 3)
	dtf$Est = round(as.numeric(dtf$Est), 3)
	row.names(dtf) = NULL
	dtf

	
	# Make the table into a plot
	#require(plotrix)
	#background_colors = matrix(data="#FFEFDB", nrow=1+nrow(dtf), ncol=ncol(dtf))
	addtable2plot(x=0.1, y=0.80, table=dtf, xjust=0, yjust=0)#, bty="o", box.col="black")#, bg=background_colors) # bg doesn't work
	title("Anagenetic\nparameters", font.main=2, cex.main=1, line=-2)
	

	# cell #4
	# Barplot of anagenetic parameters
	par(mar=c(5,3,3,1), xaxs = "r", yaxs = "r") 
	#plot(0,0, pch=".", col="white", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
	# plot(0,0, pch=".", col="white", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
	val_txt = c("d", "e", "a", "x", "n", "w", "u")
	plotvals =  BioGeoBEARS_model_object@params_table[val_txt, plotwhat]
	
	# xaxt="s" means plot (anything other than "n")
	# Get x-axis bar centers
	x_axis_bar_centers = barplot(height=plotvals, width=1, space=NULL, names.arg="", legend.text=FALSE, horiz=FALSE, col="white", border="black", xpd=TRUE, xaxt="s", yaxt="s", tck=0, plot=FALSE) 

	# Now plot for realsies
	barplot(height=plotvals, width=1, space=NULL, names.arg="", legend.text=FALSE, horiz=FALSE, col="white", border="black", xpd=TRUE, xaxt="s", yaxt="n", tck=0, yaxp=c(0,1,1)) 
	axis(side=1, at=x_axis_bar_centers, labels=val_txt, tick=FALSE, line=-1)
	title("Anagenetic\nparameters", font.main=2, cex.main=1, line=1)
	axis(side=2, at=NULL, labels=NULL, tick=TRUE, padj=1, las=0, cex.axis=1, tcl=-0.25, line=0, hadj=0.5)
	mtext(text="param. est.", side=2, line=2.1, padj=0.5, adj=0.5, las=3, cex=0.8)#, adj=0.65)
	

	
	# Cell #5 (table of cladogensis params)	
	# Cells #6 (cladogenesis params barplot), #7 (rightmost thin spacer)
	# Cells #8 (cladogenesis header) and #9-14 (cladogenesis rangesize model)
	# Plot the cladogeneis model P(size|event,ancsize) (requires 6 slots)
	plot_cladogenesis_size_probabilities(BioGeoBEARS_run_object, plotwhat=plotwhat)

	
	
	}







#######################################################
# plot_cladogenesis_size_probabilities
#######################################################
#' Graphical display of P(daughter rangesize) for your input or inferred speciation model
#' 
#' This function produces a graphical summary of the daughter rangesize aspect of the 
#' cladogenesis model stored in a \code{BioGeoBEARS_run_object}. This could be either an 
#' input model, or the result of the ML parameter search.
#'
#' The \code{LAGRANGE} DEC model assumes that at cladogenesis events, one daughter species has a 
#' range size of 1 area, and the other daughter either inherits the full ancestral range 
#' (sympatric-subset speciation), inherits the remainder of the ancestral range (vicariance),
#' or as the same range (sympatric-range copying, which is the only option when the ancestor
#' range is of size 1 area.  
#'
#' BioGeoBEARS enables numerous additional models. To see how these are similar or different from
#' the LAGRANGE DEC cladogenesis model, this function can be used.  E.g., comparison of
#' \code{LAGRANGE} DEC to a \code{DIVA}-like model is instructive: see examples. DIVA disallows
#' sympatric-subset speciation (probability 0 under this model), but allows classic vicariance 
#' (a species with 4 areas splitting into 2 daughters, each occupying 2 areas).  LAGRANGE DEC
#' gives 0 probability to a \code{4->(2,2)} history, allowing only \code{4->(3,1)} or 
#' \code{4->(1,3)} histories.
#'
#' Several additional plots relating to the cladogenesis model are also produced.  Best used via
#' \code{\link{plot_BioGeoBEARS_model}}.
#'
#' @param BioGeoBEARS_run_object The input run object.
#' @param plotwhat Default is "input", which means plotting the starting model.
#' @param statenames State names to pass to \code{\link{plot_cladogenesis_size_probabilities}}. 
#' If \code{NULL} (default), these are auto-generated assuming all areas up to the maximum number 
#' are allowed.
#' @return \code{Nothing} 
#' @export
#' @seealso \code{\link{plot_BioGeoBEARS_model}}, \code{\link{define_BioGeoBEARS_run}}, \code{\link{define_BioGeoBEARS_model_object}}
#' @note Go (BioGeo)BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' blah=1
#' 
plot_cladogenesis_size_probabilities <- function(BioGeoBEARS_run_object, plotwhat="est", statenames=NULL)
	{
	# Get the model
	BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
	
	# Get areas/states/tips
	# Get geographic ranges at tips
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=np(BioGeoBEARS_run_object$geogfn))
	
	# Get the list of geographic areas
	areas = getareas_from_tipranges_object(tipranges)
	areas_list = seq(0, length(areas)-1, 1)		# 0-base indexes
	areanames = areas
	areanames
	
	
	if (is.na(BioGeoBEARS_run_object$max_range_size))
		{
		if (is.null(BioGeoBEARS_run_object$states_list))
			{
			# Maximum range size is all areas
			max_range_size = length(areas)
			} else {
			# If not NA
			# Get max rangesize from states list
			max_range_size = max(sapply(X=BioGeoBEARS_run_object$states_list, FUN=length), na.rm=TRUE)
			}
		} else {
		# Maximum range size hard-coded
		max_range_size = BioGeoBEARS_run_object$max_range_size
		}
	max_numareas = max_range_size




	# Take the list of areas, and get list of possible states
	# (the user can manually input states if they like)
	if (is.null(BioGeoBEARS_run_object$states_list))
		{
		states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=BioGeoBEARS_run_object$include_null_range)
		states_list
		#BioGeoBEARS_run_object$states_list = states_list
		#inputs$states_list = states_list
		} else {
		states_list = BioGeoBEARS_run_object$states_list
		#BioGeoBEARS_run_object$states_list = states_list
		#inputs$states_list = states_list
		}







	# If there is a distance matrix, use the first one 
	# (non-stratified analysis, here)
	if ( (is.null(BioGeoBEARS_run_object$list_of_distances_mats) == FALSE))
		{
		distances_mat = BioGeoBEARS_run_object$list_of_distances_mats[[1]]
		} else {
		# Default is all areas effectively equidistant
		distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))
		}

	# Get the exponent on distance, apply to distances matrix
	x = BioGeoBEARS_model_object@params_table["x",plotwhat]
	dispersal_multipliers_matrix = distances_mat ^ x

	# Environmental distances
	if ( (is.null(BioGeoBEARS_run_object$list_of_envdistances_mats) == FALSE))
		{
		envdistances_mat = BioGeoBEARS_run_object$list_of_envdistances_mats[[1]]
		} else {
		# Default is all areas effectively equidistant
		envdistances_mat = matrix(1, nrow=length(areas), ncol=length(areas))
		}

	# Get the exponent on environmental distance, apply to distances matrix
	n = BioGeoBEARS_model_object@params_table["n","est"]
	dispersal_multipliers_matrix = dispersal_multipliers_matrix * envdistances_mat ^ n


	# Apply manual dispersal multipliers, if any
	# If there is a manual dispersal multipliers matrix, use the first one 
	# (non-stratified analysis, here)
	if ( (is.null(BioGeoBEARS_run_object$list_of_dispersal_multipliers_mats) == FALSE))
		{
		manual_dispersal_multipliers_matrix = BioGeoBEARS_run_object$list_of_dispersal_multipliers_mats[[1]]
		} else {
		# Default is all areas effectively equidistant
		manual_dispersal_multipliers_matrix = matrix(1, nrow=length(areas), ncol=length(areas))
		}
	
	# Get the exponent on manual dispersal multipliers
	w = BioGeoBEARS_model_object@params_table["w","est"]

	# Apply element-wise
	dispersal_multipliers_matrix = dispersal_multipliers_matrix * manual_dispersal_multipliers_matrix ^ w


	
	# Set up the instantaneous rate matrix (Q matrix)
	# someday we'll have to put "a" (anagenic range-switching) in here...
	#Qmat = rcpp_states_list_to_DEmat(areas_list=areas_list, states_list=states_list, dmat=dmat_times_d, elist=elist, amat=amat, include_null_range=BioGeoBEARS_run_object$include_null_range, normalize_TF=TRUE, makeCOO_TF=force_sparse)

	#######################################################
	# Cladogenic model
	#######################################################
	j = BioGeoBEARS_model_object@params_table["j",plotwhat]
	ysv = BioGeoBEARS_model_object@params_table["ysv",plotwhat]
	ys = BioGeoBEARS_model_object@params_table["ys",plotwhat]
	v = BioGeoBEARS_model_object@params_table["v",plotwhat]
	y = BioGeoBEARS_model_object@params_table["y",plotwhat]
	s = BioGeoBEARS_model_object@params_table["s",plotwhat]
	sum_SPweights = y + s + j + v


	maxent_constraint_01 = BioGeoBEARS_model_object@params_table["mx01",plotwhat]
	
	# Text version of speciation matrix	
	maxent_constraint_01v = BioGeoBEARS_model_object@params_table["mx01v",plotwhat]
	#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
		
	# Set the parameter controlling the size distribution of 
	# the smaller descendant species
	maxent01s_param = BioGeoBEARS_model_object@params_table["mx01s",plotwhat]
	maxent01v_param = BioGeoBEARS_model_object@params_table["mx01v",plotwhat]
	maxent01j_param = BioGeoBEARS_model_object@params_table["mx01j",plotwhat]
	maxent01y_param = BioGeoBEARS_model_object@params_table["mx01y",plotwhat]


	# Cladogenesis model inputs
	spPmat_inputs = NULL

	# Note that this gets the dispersal multipliers matrix, which is applied to 
	# e.g. the j events, NOT the dmat_times_d above which is d*dispersal_multipliers_matrix
	dmat = dispersal_multipliers_matrix
	spPmat_inputs$dmat = dmat

	states_indices = states_list
	# shorten the states_indices by 1 (cutting the 
	# null range state from the speciation matrix)
	if (include_null_range == TRUE)
		{
		states_indices[1] = NULL
		} # END if (include_null_range == TRUE)
	spPmat_inputs$l = states_indices
	spPmat_inputs$s = s
	spPmat_inputs$v = v
	spPmat_inputs$j = j
	spPmat_inputs$y = y
	spPmat_inputs$maxent01s_param = maxent01s_param
	spPmat_inputs$maxent01v_param = maxent01v_param
	spPmat_inputs$maxent01j_param = maxent01j_param
	spPmat_inputs$maxent01y_param = maxent01y_param

	
	
	
	# Calculate the ancsize/decsize speciation model...
	maxent01s = relative_probabilities_of_subsets(max_numareas=max_numareas, maxent_constraint_01=maxent01s_param, NA_val=0)
	maxent01v = relative_probabilities_of_vicariants(max_numareas=max_numareas, maxent_constraint_01v=maxent01v_param, NA_val=0)
	maxent01j = relative_probabilities_of_subsets(max_numareas=max_numareas, maxent_constraint_01=maxent01j_param, NA_val=0)
	maxent01y = relative_probabilities_of_subsets(max_numareas=max_numareas, maxent_constraint_01=maxent01y_param, NA_val=0)

	# Matrix of probs for each ancsize
	maxprob_as_function_of_ancsize_and_decsize = mapply(FUN=max, maxent01s, maxent01s, maxent01s, maxent01s, MoreArgs=list(na.rm=TRUE))
	maxprob_as_function_of_ancsize_and_decsize = matrix(data=maxprob_as_function_of_ancsize_and_decsize, nrow=nrow(maxent01s), ncol=ncol(maxent01s))
	maxprob_as_function_of_ancsize_and_decsize[maxprob_as_function_of_ancsize_and_decsize > 0] = 1
	maxprob_as_function_of_ancsize_and_decsize[maxprob_as_function_of_ancsize_and_decsize <= 0] = 0

	# Now, go through, and make a list of the max minsize for each decsize
	max_minsize_as_function_of_ancsize = apply(X=maxprob_as_function_of_ancsize_and_decsize, MARGIN=1, FUN=maxsize)
	
	rangesizes = 1:nrow(maxent01s)
	rangesizes_at_y = rev(rangesizes - 0.5)
	rangesizes_at_x = rangesizes - 0.5

	maxent01y = adf2(round(maxent01y,2))
	#rownames(maxent01y) = paste("ancsize=", rangesizes, sep="")
	rownames(maxent01y) = rangesizes
	colnames(maxent01y) = rangesizes
	maxent01s = adf2(round(maxent01s,2))
	#rownames(maxent01s) = paste("ancsize=", rangesizes, sep="")
	rownames(maxent01s) = rangesizes
	colnames(maxent01s) = rangesizes
	maxent01v = adf2(round(maxent01v,2))
	#rownames(maxent01v) = paste("ancsize=", rangesizes, sep="")
	rownames(maxent01v) = rangesizes
	colnames(maxent01v) = rangesizes
	maxent01j = adf2(round(maxent01j,2))
	#rownames(maxent01j) = paste("ancsize=", rangesizes, sep="")
	rownames(maxent01j) = rangesizes
	colnames(maxent01j) = rangesizes

	# Set to 0, if the parameter is 0
	y_cellcolors = NA
	s_cellcolors = NA
	v_cellcolors = NA
	j_cellcolors = NA
	if (y == 0)
		{
		maxent01y[maxent01y != 0] = 0
		y_cellcolors = matrix(data="white", nrow=nrow(maxent01y), ncol=ncol(maxent01y))
		}
	if (s == 0)
		{
		maxent01s[maxent01s != 0] = 0
		s_cellcolors = matrix(data="white", nrow=nrow(maxent01s), ncol=ncol(maxent01s))
		}
	if (v == 0)
		{
		maxent01v[maxent01v != 0] = 0
		v_cellcolors = matrix(data="white", nrow=nrow(maxent01v), ncol=ncol(maxent01v))
		}
	if (j == 0)
		{
		maxent01j[maxent01j != 0] = 0
		j_cellcolors = matrix(data="white", nrow=nrow(maxent01j), ncol=ncol(maxent01j))
		}



	# Cell #5
	# Return to no margins
	par(mar=c(0,0,0,0), xaxs = "i", yaxs = "i") 
	
	# cells # 5, 6
	#plot(0,0, pch=".", col="white", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
	# Table of free vs. fixed params
	plot(x=c(0,1),y=c(0,1), pch=".", col="white", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
	val_txt = c("ysv", "ys", "y", "s", "v", "mx01", "mx01y", "mx01s", "mx01v", "mx01j")
	Param = val_txt
	Type = BioGeoBEARS_model_object@params_table[val_txt, "type"]
	Init = BioGeoBEARS_model_object@params_table[val_txt, "init"]
	Est = BioGeoBEARS_model_object@params_table[val_txt, "est"]
	dtf = adf2(cbind(Param, Type, Init, Est))
	dtf = dfnums_to_numeric(dtf)
	row.names(dtf) = NULL
	dtf$Init = round(dtf$Init, 2)
	dtf$Est = round(dtf$Est, 2)
	dtf
	
	# Make the table into a plot
	#require(plotrix)
	addtable2plot(x=0.1, y=0.80, table=dtf, xjust=0, yjust=0, cex=1)
	title("Cladogenetic\nparameters", font.main=2, cex.main=1, line=-2)


	# Row #2, right side: speciation model per-event probabilities
	# Rightmost square cell
	#plot(0,0, pch=".", col="white", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
	par(mar=c(5,3,3,1), xaxs = "r", yaxs = "r") 
	
	
	# If you just want to plot the default (unscaled) y,s,v,j values:
	plot_unscaled_ysvj = FALSE
	if (plot_unscaled_ysvj)
		{
		barplot_yaxis = c(0, 0.5, 1)
		barplot_yaxis_txt = c(0, 0.5, 1)
		ylims = c(0,1.35)
		
		val_txt = c("y", "s", "v", "j")
		plotvals =  BioGeoBEARS_model_object@params_table[val_txt, plotwhat]
		#plotvals = c(y,s,v,j)
		
		# xaxt="s" means plot (anything other than "n")
	
		# Get x-axis bar centers
		x_axis_bar_centers = barplot(height=plotvals, width=1, space=NULL, names.arg="", legend.text=FALSE, horiz=FALSE, col="white", border="black", ylim=ylims, xpd=TRUE, xaxt="s", yaxt="n", tck=0, yaxp=c(0,1,1), plot=FALSE) 

		# Now plot for realsies
		barplot(height=plotvals, width=1, space=NULL, names.arg="", legend.text=FALSE, horiz=FALSE, col="white", border="black", ylim=ylims, xpd=TRUE, xaxt="s", yaxt="n", tck=0, yaxp=c(0,1,1)) 
		axis(side=1, at=x_axis_bar_centers, labels=c("y", "s", "v", "j"), tick=FALSE, line=-1)
		title("Relative prob. of\neach type of\ncladogenesis event", font.main=2, cex.main=1)
		axis(side=2, at=barplot_yaxis, labels=barplot_yaxis_txt, tick=TRUE, padj=0.5, las=1, cex.axis=1, tcl=-0.25, line=0, hadj=0.5)
		mtext(text="relative prob.", side=2, line=2.1, padj=0.5, adj=0.5, las=3, cex=0.8)#, adj=0.65)
		
		} else {
		# If you ARE rescaling, e.g. so ysvj add up to 1
		barplot_yaxis = c(0, 0.5, 1)
		barplot_yaxis_txt = c(0, 0.5, 1)
		ylims = c(0, 1.15)

		params_table = BioGeoBEARS_model_object@params_table
		tmpvals = get_clado_perEvent_weights(params_table, plotwhat=plotwhat)
		
		plotvals = c(tmpvals$y, tmpvals$s, tmpvals$v, tmpvals$j)
		# xaxt="s" means plot (anything other than "n")
	
		# Get x-axis bar centers
		x_axis_bar_centers = barplot(height=plotvals, width=1, space=NULL, names.arg="", legend.text=FALSE, horiz=FALSE, col="white", border="black", ylim=ylims, xpd=TRUE, xaxt="s", yaxt="n", tck=0, yaxp=c(0,1,1), plot=FALSE) 

		# Now plot for realsies
		barplot(height=plotvals, width=1, space=NULL, names.arg="", legend.text=FALSE, horiz=FALSE, col="white", border="black", ylim=ylims, xpd=TRUE, xaxt="s", yaxt="n", tck=0, yaxp=c(0,1,1)) 
		axis(side=1, at=x_axis_bar_centers, labels=c("y", "s", "v", "j"), tick=FALSE, line=-1)
		title("Relative prob. of\neach type of\ncladogenesis event", font.main=2, cex.main=1)
		axis(side=2, at=barplot_yaxis, labels=barplot_yaxis_txt, tick=TRUE, padj=0.5, las=1, cex.axis=1, tcl=-0.25, line=0, hadj=0.5)
		mtext(text="relative prob.", side=2, line=2.1, padj=0.5, adj=0.5, las=3, cex=0.8)#, adj=0.65)
		
		} # end parplot in right square cell of row #2.

	# Plot text of the numbers, above the bars, if desired
	plot_text = TRUE
	if (plot_text == TRUE)
		{
		# pos=3 means above
		text(x=x_axis_bar_centers, y=plotvals, labels=round(plotvals,2), pos=3, offset=0.2, cex=0.9)
		}
	
	
	# Rightmost thin cell
	par(mar=c(0,0,0,0), xaxs = "i", yaxs = "i") 
	plot(0,0, pch=".", col="white", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")


	# Row 3: title for speciation model
	par(mar=c(0,0,0,0), xaxs = "i", yaxs = "i") 
	plot(x=c(0,1), y=c(0,1), pch=".", col="white", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
	tmptxt = paste("Prob(smaller daughter range size, given ancestor range size)", sep="")
	#text(x=0.5, y=0.5, tmptxt, cex=1.2)
	title(tmptxt, cex=1.5, line=-1)

	# Row 4 left thin column
	plot(x=c(0,1), y=c(0,1), pch=".", col="white", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
	#mtext(text="ancsize", side=2, line=-3, hadj=0.5, padj=1, las=3, cex=1)
	mtext(text="ancsize", side=2, line=-3, adj=0.65, padj=1, las=3, cex=0.8)
	
	# Plot maxent
	par(mar=c(5,3,3,1), xaxs = "r", yaxs = "r") 
	
	# for color2D.matplot
	#require(plotrix)	# for color2D.matplot
	
	color2D.matplot(x=maxent01y, c(1,0), c(1,0), c(1,0), cellcolors=y_cellcolors, show.values=TRUE, show.legend=FALSE, xlab="", ylab="", axes=FALSE, nslices=100, xaxt="n", yaxt="n")
	axis(side=2, at=rangesizes_at_y, labels=rangesizes, tick=FALSE, padj=0.5, las=1, cex.axis=1.5, line=-0.8)
	axis(side=3, at=rangesizes_at_x, labels=rangesizes, tick=FALSE, hadj=0.5, cex.axis=1.5, line=-0.8)			
	title("y:Sympatric (copying)", line=2, font.main=2, cex.main=1.1)
	mtext("decsize", side=1, cex=0.8, line=0.2, font.main=1, cex.main=0.9)
	
	color2D.matplot(x=maxent01s, c(1,0), c(1,0), c(1,0), cellcolors=s_cellcolors, show.values=TRUE, show.legend=FALSE, xlab="", ylab="", axes=FALSE, nslices=100, xaxt="n", yaxt="n")
	axis(side=2, at=rangesizes_at_y, labels=rangesizes, tick=FALSE, padj=0.5, las=1, cex.axis=1.5, line=-0.8)
	axis(side=3, at=rangesizes_at_x, labels=rangesizes, tick=FALSE, hadj=0.5, cex.axis=1.5, line=-0.8)			
	title("s:Sympatric (subset)", line=2, font.main=2, cex.main=1.1)
	mtext("decsize", side=1, cex=0.8, line=0.2, font.main=1, cex.main=0.9)
	
	color2D.matplot(x=maxent01v, c(1,0), c(1,0), c(1,0), cellcolors=v_cellcolors, show.values=TRUE, show.legend=FALSE, xlab="", ylab="", axes=FALSE, nslices=100, xaxt="n", yaxt="n")
	axis(side=2, at=rangesizes_at_y, labels=rangesizes, tick=FALSE, padj=0.5, las=1, cex.axis=1.5, line=-0.8)
	axis(side=3, at=rangesizes_at_x, labels=rangesizes, tick=FALSE, hadj=0.5, cex.axis=1.5, line=-0.8)			
	title("v:Vicariance", line=2, font.main=2, cex.main=1.1)
	mtext("decsize", side=1, cex=0.8, line=0.2, font.main=1, cex.main=0.9)
	
	color2D.matplot(x=maxent01j, c(1,0), c(1,0), c(1,0), cellcolors=j_cellcolors, show.values=TRUE, show.legend=FALSE, xlab="", ylab="", axes=FALSE, nslices=100, xaxt="n", yaxt="n")
	axis(side=2, at=rangesizes_at_y, labels=rangesizes, tick=FALSE, padj=0.5, las=1, cex.axis=1.5, line=-0.8)
	axis(side=3, at=rangesizes_at_x, labels=rangesizes, tick=FALSE, hadj=0.5, cex.axis=1.5, line=-0.8)			
	title("j:Founder-event (jump)", line=2, font.main=2, cex.main=1.1)
	mtext("decsize", side=1, cex=0.8, line=0.2, font.main=1, cex.main=0.9)
	
	
	par(mar=c(0,0,0,0), xaxs = "i", yaxs = "i") 
	# Row 2 right thin column
	plot(x=c(0,1), y=c(0,1), pch=".", col="white", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
	testcol = rev(color.gradient(c(0,1),c(0,1),c(0,1),nslices=100))
	color.legend(xl=0.15, yb=0.3, xr=0.4, yt=0.8, legend=c(0,1), align="rb", rect.col=testcol, gradient="y")	# gradient in x-axis
	


	
	#######################################################
	# Bottom row: depict the cladogenesis process in some fashion.
	# Let's just do the rowSums of the cladogenesis matrix
	#######################################################
	# Rangesize of each state
	tmpstates = states_list
	if ((length(tmpstates[[1]]) == 1) && (is.na(tmpstates[[1]]) == TRUE))
		{
		tmpstates[[1]] = NULL
		}
	rangesizes = sapply(X=tmpstates, FUN=length, simplify=TRUE)
	rangesizes

	# Footer
	#par(mar=c(0,0,0,0), xaxs = "i", yaxs = "i") 
	par(mar=c(5,3,3,1), xaxs = "r", yaxs = "r") 

	tmpca_1 = rep(1, (ncol(tipranges@df)-1))
	tmpcb_1 = rep(1, (ncol(tipranges@df)-1))

	COO_weights_columnar = rcpp_calc_anclikes_sp_COOweights_faster(Rcpp_leftprobs=tmpca_1, Rcpp_rightprobs=tmpcb_1, l=tmpstates, s=s, v=v, j=j, y=y, dmat=dmat, maxent01s=maxent01s, maxent01v=maxent01v, maxent01j=maxent01j, maxent01y=maxent01y, max_minsize_as_function_of_ancsize=max_minsize_as_function_of_ancsize, printmat=FALSE)
	COO_weights_columnar

	# This gives 15 states
	Rsp_rowsums = rcpp_calc_rowsums_for_COOweights_columnar(COO_weights_columnar=COO_weights_columnar)
	Rsp_rowsums

	# Do a count of nonzeros
	COO_weights_columnar_count = COO_weights_columnar
	COO_weights_columnar_count[[4]][COO_weights_columnar_count[[4]] > 0] = 1
	Rsp_rowsums_count = rcpp_calc_rowsums_for_COOweights_columnar(COO_weights_columnar=COO_weights_columnar_count)
	Rsp_rowsums_count
	
	
	
	
	# Make a plot of ancestral range size, versus number of possible descendent pairs
	ylims = c(0, max(max(Rsp_rowsums_count), max(rangesizes)))
	xlims = c(0, max(rangesizes))
	#plot(x=rangesizes, y=Rsp_rowsums, xlim=xlims, ylim=ylims, xlab="anc. range size", ylab="# desc. pairs with prob. > 0", main="Bias towards widespread ancestors")
	
	plot(x=rangesizes, y=Rsp_rowsums_count, xlim=xlims, ylim=ylims, main="", xaxt="n", yaxt="n", tck=0, ylab="", xlab="")
	axis(side=1, at=NULL, labels=NULL, tick=TRUE, padj=-1, las=0, cex.axis=1, tcl=-0.25, line=0, hadj=0.5)
	axis(side=2, at=NULL, labels=NULL, tick=TRUE, padj=1, las=3, cex.axis=1, tcl=-0.25, line=0, hadj=0.5)

	#mtext(text="param. est.", side=2, line=2.1, padj=0.5, adj=0.5, las=3, cex=0.8)#, adj=0.65)

	mtext("anc. range size", side=1, cex=0.8, line=2, font.main=1, cex.main=0.8)
	mtext("# desc. pairs with prob>0", side=2, line=2.1, padj=0.5, adj=0.5, las=3, cex=0.8)#, adj=0.65)
	title("Bias towards\nwidespread ancestors", font.main=2, cex.main=1)




	#######################################################
	# In the last plot, depict a bit of the cladogenesis matrix
	#######################################################
	plot(x=c(0,1), y=c(0,1), pch=".", col="white", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")

	maxrows = 6
	if (length(rangesizes) < maxrows)
		{
		maxrows = length(rangesizes)
		}
	
	
	
	
	# Just take the 1st 8 rows of the COO_weights_columnar, and the last 8
	
	# First part of speciation matrix to print
	ancstates_1st = COO_weights_columnar[[1]][1:maxrows] + 1
	leftstates_1st = COO_weights_columnar[[2]][1:maxrows] + 1
	rightstates_1st = COO_weights_columnar[[3]][1:maxrows] + 1
	relprobs_1st = round(as.numeric(COO_weights_columnar[[4]][1:maxrows] / Rsp_rowsums[ancstates_1st]), 3)
	
	# Last part of speciation matrix to print
	tmpend = length(COO_weights_columnar[[1]])
	tmpstart = tmpend - maxrows
	ancstates_2nd = COO_weights_columnar[[1]][tmpstart:tmpend] + 1
	leftstates_2nd = COO_weights_columnar[[2]][tmpstart:tmpend] + 1
	rightstates_2nd = COO_weights_columnar[[3]][tmpstart:tmpend] + 1
	relprobs_2nd = round((COO_weights_columnar[[4]][tmpstart:tmpend] / Rsp_rowsums[ancstates_2nd]), 3)

	

	
	# Now make a big matrix, with a gap in the middle
	if (is.null(statenames) == TRUE)
		{
		# Fix include_null_range here to false, for plotting
		# cladogenesis probabilities
		statenames = areas_list_to_states_list_new(areas=areanames, maxareas=max_range_size, include_null_range=FALSE, split_ABC=FALSE)
		statenames
		}
	first_mat = rbind(statenames[ancstates_1st], statenames[leftstates_1st], statenames[rightstates_1st], relprobs_1st)
	second_mat = rbind(statenames[ancstates_2nd], statenames[leftstates_2nd], statenames[rightstates_2nd], relprobs_2nd)
	spacer = matrix(" ", nrow=4, ncol=1)
	printmat = cbind(first_mat, spacer, second_mat)
	printmat = adf2(printmat)
	names(printmat) = c(paste("", 1:maxrows, sep=""), "_", paste("", tmpstart:tmpend, sep=""))
	rownames(printmat) = c("Ancestor", "Left desc.", "Right desc.", "Rel. prob.")
	printmat
	
	
	addtable2plot(printmat, x=0, y=0.8, table=printmat, display.rownames=TRUE, display.colnames=TRUE, bty="o", hlines=TRUE, vlines=TRUE, xjust=0, yjust=0, cex=1.1)
	title("Conditional probabilities of example cladogenesis events", font.main=2, cex.main=1)
	
	} # END plot_cladogenesis_size_probabilities














#######################################################
# Convert stochastic map to state probabilities for 
# plotting with plot_BioGeoBEARS_results()
#######################################################
stochastic_map_states_into_res <- function(res, master_table_cladogenetic_events, stratified=FALSE)
	{
	defaults='
	resmod = res
	master_table_cladogenetic_events = stochastic_mapping_results2$master_table_cladogenetic_events
	' # 
	
	# Load the tree
	tr = ape::read.tree(file=res$inputs$trfn)
	
	# Error checks -- is this a stratified analysis?
	strat_TF = ("SUBnode.type" %in% names(master_table_cladogenetic_events))
	if ( (strat_TF == TRUE) && (stratified == FALSE) )
		{
		errortxt = paste("\n\nError in stochastic_map_states_into_res():\n\nYou said 'stratified=FALSE' in your inputs (perhaps implicitly), but\nyour master_table_cladogenetic_events has e.g. 'SUBnode.type' in the column names,\nindicating a stratified analysis.\n\n", sep="")
		cat(errortxt)
		
		stop("Stopping on error.")
		}
		
	if ( (strat_TF == FALSE) && (stratified == TRUE) )
		{
		errortxt = paste("\n\nError in stochastic_map_states_into_res():\n\nYou said 'stratified=TRUE' in your inputs , but\nyour master_table_cladogenetic_events LACKS e.g. 'SUBnode.type' in the column names,\nindicating a non-stratified analysis.\n\n", sep="")
		cat(errortxt)
		
		stop("Stopping on error.")
		}	
	
	
	#######################################################
	# NOTE: YOU CANNOT REALLY PLOT THE 'BRANCH BOTTOM' STATES
	# AT BRANCH BOTTOMS, SINCE THEY ARE THE BOTTOM OF THE 
	# SUBTREE ROOT BRANCHES, NOT THE BOTTOM OF THE BRANCHES IN THE MASTER TREE
	# (although they will often, not always, be the same)
	#######################################################
	resmod = res

	# Input the sampled node states into the state probabilities
	resmod$ML_marginal_prob_each_state_at_branch_top_AT_node
	
	# Get the main nodes on the master tree from the master table
	# This has to be done somewhat differently
	if (stratified == TRUE)
		{
		# Get the main nodes on the master tree
		TF1 = master_table_cladogenetic_events$node.type != "tip"
		TF2 = master_table_cladogenetic_events$SUBnode.type != "tip"
		TF3 = master_table_cladogenetic_events$piececlass == "subtree"
		TF = ((TF1 + TF2 + TF3) == 3)
		sum(TF)
		} else {
		# Get the main nodes on the master tree
		TF = master_table_cladogenetic_events$node.type != "tip"
		sum(TF)
		} # END if (stratified == TRUE)
		
	# Look at the cladogenetic events table
	lastcol = ncol(master_table_cladogenetic_events)
	rownums_for_nodes_in_master_tree = (1:nrow(master_table_cladogenetic_events))[TF]
	internal_nodes = master_table_cladogenetic_events$node[rownums_for_nodes_in_master_tree]
	order_rows = order(internal_nodes)
	internal_nodes = sort(internal_nodes)
	master_table_cladogenetic_events[rownums_for_nodes_in_master_tree[order_rows],-lastcol]

	# Look at the sampled states
	master_table_cladogenetic_events$sampled_states_AT_nodes[rownums_for_nodes_in_master_tree[order_rows]]
	master_table_cladogenetic_events$samp_LEFT_dcorner[rownums_for_nodes_in_master_tree[order_rows]]
	master_table_cladogenetic_events$samp_RIGHT_dcorner[rownums_for_nodes_in_master_tree[order_rows]]

	# Zero out the old probabilities, and put in 1s for the sampled states

	# Nodes
	resmod$ML_marginal_prob_each_state_at_branch_top_AT_node[internal_nodes,] = 0
	states_1based_to_change = master_table_cladogenetic_events$sampled_states_AT_nodes[rownums_for_nodes_in_master_tree[order_rows]]
	resmod$ML_marginal_prob_each_state_at_branch_top_AT_node[internal_nodes,states_1based_to_change] = 1

	# Corners
#	if (strat_TF == FALSE)
#		{
		leftright_nodes_matrix = matrix(data=unlist(master_table_cladogenetic_events[rownums_for_nodes_in_master_tree[order_rows],]$daughter_nds), ncol=2, byrow=TRUE)
		Lcorners_above_nodenums = leftright_nodes_matrix[,2]
		Rcorners_above_nodenums = leftright_nodes_matrix[,1]
#		} else {
#		tr = read.tree(res$inputs$trfn)
#		tr_pruningwise = reorder(tr, "pruningwise")
#		leftright_nodes_matrix = get_leftright_nodes_matrix_from_results(tr_pruningwise)
#		Lcorners_above_nodenums = leftright_nodes_matrix[,2]
#		Rcorners_above_nodenums = leftright_nodes_matrix[,1]
#		}

	# Internal node states
	resmod$ML_marginal_prob_each_state_at_branch_top_AT_node[internal_nodes,] = 0
	states_1based_to_change = master_table_cladogenetic_events$sampled_states_AT_nodes[rownums_for_nodes_in_master_tree[order_rows]]
	for (i in 1:length(internal_nodes))
		{
		resmod$ML_marginal_prob_each_state_at_branch_top_AT_node[internal_nodes[i],states_1based_to_change[i]] = 1
		} # END for (i in 1:length(internal_nodes))


	# Get the row numbers in the master table corresponding to the 
	# standard tips and nodes in the unsliced master tree
	ntips = length(tr$tip.label)
	rownums_for_allnodes_in_master_tree = c(1:ntips, rownums_for_nodes_in_master_tree[order_rows])

	resmod$ML_marginal_prob_each_state_at_branch_bottom_below_node[,] = NA
	
	
	#######################################################
	# JUNK TO CUT
	#######################################################
	junk='
	# Branch bottoms OF SUBTREES 
	# (useful, but not for plotting states)
	for (i in 1:length(Lcorners_above_nodenums))
		{

		# Change left corners
		states_1based_to_change = master_table_cladogenetic_events$sampled_states_AT_brbots[rownums_for_allnodes_in_master_tree][Lcorners_above_nodenums[i]]
		resmod$ML_marginal_prob_each_state_at_branch_bottom_below_node[Lcorners_above_nodenums[i],states_1based_to_change] = 1
		resmod$ML_marginal_prob_each_state_at_branch_bottom_below_node[Lcorners_above_nodenums[i],-states_1based_to_change] = 0

		# Right corners
		states_1based_to_change = master_table_cladogenetic_events$sampled_states_AT_brbots[rownums_for_allnodes_in_master_tree][Rcorners_above_nodenums[i]]
		resmod$ML_marginal_prob_each_state_at_branch_bottom_below_node[Rcorners_above_nodenums[i],states_1based_to_change] = 1
		resmod$ML_marginal_prob_each_state_at_branch_bottom_below_node[Rcorners_above_nodenums[i],-states_1based_to_change] = 0
		}
	' #END junk
	#######################################################
	# END JUNK TO CUT
	#######################################################


	for (i in 1:length(internal_nodes))
		{

		# Change left corners
		states_1based_to_change = master_table_cladogenetic_events$samp_LEFT_dcorner[rownums_for_allnodes_in_master_tree][internal_nodes[i]]
		resmod$ML_marginal_prob_each_state_at_branch_bottom_below_node[Rcorners_above_nodenums[i],states_1based_to_change] = 1
		resmod$ML_marginal_prob_each_state_at_branch_bottom_below_node[Rcorners_above_nodenums[i],-states_1based_to_change] = 0

		# Right corners
		states_1based_to_change = master_table_cladogenetic_events$samp_RIGHT_dcorner[rownums_for_allnodes_in_master_tree][internal_nodes[i]]
		resmod$ML_marginal_prob_each_state_at_branch_bottom_below_node[Lcorners_above_nodenums[i],states_1based_to_change] = 1
		resmod$ML_marginal_prob_each_state_at_branch_bottom_below_node[Lcorners_above_nodenums[i],-states_1based_to_change] = 0
		}



	resmod$ML_marginal_prob_each_state_at_branch_bottom_below_node


	# States
	#analysis_titletxt = "Stochastic map"
	#scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
	#res2 = plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=plotsplits, cornercoords_loc=scriptdir, include_null_range=BioGeoBEARS_run_object$include_null_range, tr=tr, tipranges=tipranges)

	return(resmod)
	} # END plot_stochastic_map_states



get_colors_for_states_list_0based <- function(areanames, states_list_0based=NULL, max_range_size=NA, plot_null_range=FALSE)
	{
	defaults='
	areanames=c("K", "O", "M", "H")
	states_list_0based=NULL
	max_range_size=NA
	include_null_range=TRUE
	'
	
	numareas = length(areanames)
	numareas
	
	if (is.na(max_range_size))
		{
		max_range_size = length(areanames)
		}
	max_range_size

	if (is.null(states_list_0based))
		{
		states_list_0based = rcpp_areas_list_to_states_list(areas=areanames, maxareas=max_range_size, include_null_range=plot_null_range)
		}

	# Set up colors
	colors_matrix = get_colors_for_numareas(length(areanames))
	states_list_0based_index = states_list_0based
	
	# Run it 
	# Fix plot_null_range
	colors_list_for_states = mix_colors_for_states(colors_matrix, states_list_0based_index=states_list_0based_index, plot_null_range=plot_null_range)
	colors_list_for_states
	
	return(colors_list_for_states)
	}




# tr is required input (since users could change it in a variety of ways)...

paint_stochastic_map_branches <- function(res, master_table_cladogenetic_events, colors_list_for_states, tr=NULL, lwd=5, lty=par("lty"), root.edge=TRUE, cornercoords_loc=np(system.file("extdata/a_scripts", package="BioGeoBEARS")), stratified=FALSE, plot_clado_1desc=FALSE, plot_clado_1desc_points=FALSE, thickness_by_area=FALSE, states_list_0based=NULL, include_null_range=TRUE)
	{
	defaults='
	res = resDEC
	master_table_cladogenetic_events = stochastic_mapping_results$master_table_cladogenetic_events

	plot_stochastic_map_states_stratified(res, tr, master_table_cladogenetic_events, plotsplits=TRUE)

	scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
	#cornercoords_loc="manual"
	cornercoords_loc=scriptdir
	root.edge = TRUE
	lwd = 5
	lty=par("lty")
	
	
	# Get colors_list_for_states
	# Setup 
	#include_null_range = TRUE
	areanames = c("K", "O", "M", "H")
	areas = areanames
	max_range_size = 4
	states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)

	plot_clado_1desc=FALSE
	plot_clado_1desc_points=FALSE
	thickness_by_area=FALSE

	# Get colors
	plot_null_range = FALSE
	colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=max_range_size, plot_null_range=plot_null_range)
	' # END defaults

	# Keep a 2-column list of the lengths in each state
	list_of_each_state = NULL
	lengths_in_each_state = NULL
	
	
	# Read a tree, if needed (but bad idea!)
	if (is.null(tr))
		{
		#tr = read.tree(res$inputs$trfn)
		tr = check_trfn(trfn=res$inputs$trfn)
		} # END if (is.null(tr))

	# If we want to paint by thickness, make a list of range areas
	if (thickness_by_area == FALSE)
		{
		rangesizes = rep(lwd, times=length(colors_list_for_states))
		} else {
		if (is.null(states_list_0based) == TRUE)
			{
			errortxt = paste("\n\nERROR in paint_stochastic_map_branches(): 'thickness_by_area' is TRUE,\nbut if you want to paint branches with their thickness, you need to supply 'states_list_0based', which you did not.\nHave a hoopy froody day.\n\n", sep="")
			cat(errortxt)
			rangesizes = rep(lwd, times=length(colors_list_for_states))
			#stop(errortxt)
			} else {
			# Calculate the number of areas
			
			# Plotted line with is lwd * numareas
			rangesizes = lwd*sapply(FUN=length, X=states_list_0based)
			rangesizes
			
			# Set the null range to size 0 areas
			if (is.na(states_list_0based[[1]]))
				{
				rangesizes[[1]] = 0
				} # if (is.na(states_list_0based[[1]]))
			} # if (is.null(states_list_0based))
		} # if (is.null(thickness_by_area) == FALSE)



	# Error checks -- is this a stratified analysis?
	strat_TF = ("SUBnode.type" %in% names(master_table_cladogenetic_events))
	if ( (strat_TF == TRUE) && (stratified == FALSE) )
		{
		errortxt = paste("\n\nError in paint_stochastic_map_branches():\n\nYou said 'stratified=FALSE' in your inputs (perhaps implicitly), but\nyour master_table_cladogenetic_events has e.g. 'SUBnode.type' in the column names,\nindicating a stratified analysis.\n\n", sep="")
		cat(errortxt)
		
		stop("Stopping on error.")
		}
		
	if ( (strat_TF == FALSE) && (stratified == TRUE) )
		{
		errortxt = paste("\n\nError in paint_stochastic_map_branches():\n\nYou said 'stratified=TRUE' in your inputs , but\nyour master_table_cladogenetic_events LACKS e.g. 'SUBnode.type' in the column names,\nindicating a non-stratified analysis.\n\n", sep="")
		cat(errortxt)
		
		stop("Stopping on error.")
		}	


	# Get the coordinates of the nodes
	nodes_xy = node_coords(tr, coords_fun="plot_phylo3_nodecoords", tmplocation=cornercoords_loc, root.edge=root.edge)
	nodes_xy

	# Get the maximum height / x value of the tree
	max_x = max(nodes_xy$x)

	# Plot branches with no change
	naTF = is.na(master_table_cladogenetic_events$anagenetic_events_txt_below_node) == TRUE
	master_table_cladogenetic_events$anagenetic_events_txt_below_node[naTF] = "none"
	nochange_TF = master_table_cladogenetic_events$anagenetic_events_txt_below_node == "none"
	sum(nochange_TF)

	# Plot branch histories with no change
	for (rownum in 1:nrow(master_table_cladogenetic_events))
		{
		
		# Check for cladogenetic events at this node
		# (all nodes except tips)
		# Do the test differently depending on stratified or not
		if (stratified == TRUE)
			{
			test_TF = ((master_table_cladogenetic_events$SUBnode.type[rownum] == "internal") || (master_table_cladogenetic_events$SUBnode.type[rownum] == "root") )
			} else {
			test_TF = ((master_table_cladogenetic_events$node.type[rownum] == "internal") || (master_table_cladogenetic_events$node.type[rownum] == "root") )
			}
		
		
		# Plot the sampled cladogenesis events
		if (test_TF == TRUE )
			{
			# Get the nodenums of the descendants
			master_nodenum = master_table_cladogenetic_events$node[rownum]
			daughter_nds = master_table_cladogenetic_events$daughter_nds[rownum][[1]]
			Ldesc_nodenum = daughter_nds[1]
			Rdesc_nodenum = daughter_nds[2]
		
			# Draw the left descendant
			# Get the corresponding plot coordinates
			time_bp = master_table_cladogenetic_events$time_bp[rownum]
			start_x = max_x - time_bp
			end_x = max_x - time_bp
		
			# Y-coord at node
			start_y = nodes_xy$y[master_nodenum]
		
			# Left Node
			# Y-coord at corner (is just the y-coord of the desc. node)
			end_y = nodes_xy$y[Ldesc_nodenum]
		
			#######################################################
			# PAINT THE LEFT BRANCH (NODE TO CORNER)
			#######################################################
			# Plot the segment
			statenum_1based = as.numeric(master_table_cladogenetic_events$samp_LEFT_dcorner[rownum])
			tmpcolor = colors_list_for_states[statenum_1based]
			lwd_tmp = rangesizes[statenum_1based]
			segments(x0=start_x, x1=end_x, y0=start_y, y1=end_y, col=tmpcolor, lwd=lwd_tmp, lty=lty, lend="round")
			
			# Right Node
			# Y-coord at corner (is just the y-coord of the desc. node)
			end_y = nodes_xy$y[Rdesc_nodenum]
		
			#######################################################
			# PAINT THE RIGHT BRANCH (NODE TO CORNER)
			#######################################################
			# Plot the segment
			statenum_1based = as.numeric(master_table_cladogenetic_events$samp_RIGHT_dcorner[rownum])
			tmpcolor = colors_list_for_states[statenum_1based]
			lwd_tmp = rangesizes[statenum_1based]
			segments(x0=start_x, x1=end_x, y0=start_y, y1=end_y, col=tmpcolor, lwd=lwd_tmp, lty=lty, lend="round")
			} # END plot the sampled cladogenesis events
	
		# If there is no change on this piece of branch
		if (master_table_cladogenetic_events$anagenetic_events_txt_below_node[rownum] == "none")
			{
			starting_state_1based = master_table_cladogenetic_events$sampled_states_AT_brbots[rownum]
			ending_state_1based = master_table_cladogenetic_events$sampled_states_AT_nodes[rownum]
		
			# Get the relative start and stop of this "event" (non-event)
		
			# Get the master nodenum
			master_nodenum = master_table_cladogenetic_events$node[rownum]
			start_y = nodes_xy$y[master_nodenum]
			end_y = nodes_xy$y[master_nodenum]
		
			# Get the x of the top of the branch ( think...
			# The SUB times of this piece are relative to 
			# reltimept
			# This is absolute time before present
			if (stratified == TRUE)
				{
				time_top = master_table_cladogenetic_events$time_top[rownum]
				# This is absolute time before present for the branch piece we are
				# looking at
				branch_piece_TOP_time_bp = time_top + master_table_cladogenetic_events$SUBtime_bp[rownum]
				} else {
				time_top = master_table_cladogenetic_events$time_bp[rownum]
				branch_piece_TOP_time_bp = time_top
				} # END if (stratified == TRUE)

		

			# 2014-05-25_NJM: check if 
			# Is reltimept smaller?
			# Setup
			trtable = master_table_cladogenetic_events
			
			if (stratified == TRUE)
				{
				# Check for NA on branch length below node (e.g. at root)
				if (is.na(trtable$SUBedge.length[rownum]))
					{
					trtable$SUBedge.length[rownum] = 0
					}
		
				# Check if sub-edge length and edge length are the same:
				subedge_length_equals_edge_length_WORRY = FALSE
				if (trtable$SUBedge.length[rownum] > trtable$reltimept[rownum])
					{
					subedge_length_equals_edge_length_WORRY = TRUE

					# We can fix sub-branches easily:
					if (trtable$piececlass[rownum] == "subbranch")
						{
						# Fix to reltimept IF it's *NOT* a fossil:
						if ( (is.na(trtable$fossils[rownum])) || (trtable$fossils[rownum] == FALSE) )
							{
							brlen_in_section = trtable$reltimept[rownum]
							subedge_length_equals_edge_length_WORRY = FALSE
							} else {
							# It *IS* a fossil, it's brlen is based on time_bp
							brlen_in_section = trtable$time_bot[rownum] - trtable$time_bp[rownum]
							subedge_length_equals_edge_length_WORRY = FALSE
							}
						} # END check for fossils
					} else {
					brlen_in_section = trtable$SUBedge.length[rownum]
					}

				if (subedge_length_equals_edge_length_WORRY == TRUE)
					{
					errortxt = paste("\n\nError in plotting_stochastic_maps() (NO CHANGE ON BRANCH): your master_tree table, at row 'rownum'=", rownum, "\nhas an SUBedge.length > reltimept and was not corrected, as it's not a subbranch.\n\n", sep="")
					cat(errortxt)

					cat("\n\n")
					print("rownum:")
					cat("\n\n")
					print(rownum)
					cat("\n\n")
					print("trtable[rownum,]:")
					cat("\n\n")
					print(trtable[rownum,])
					cat("\n\n")
					}
				} else {
				# Stratified == FALSE
				# Check for NA on branch length below node (e.g. at root)
				if (is.na(trtable$edge.length[rownum]))
					{
					trtable$edge.length[rownum] = 0
					}
				# If stratified == FALSE, brlen is just the edge.length
				brlen_in_section = trtable$edge.length[rownum]
				} # END if (stratified == TRUE)

		
			branch_piece_BOT_time_bp = branch_piece_TOP_time_bp + brlen_in_section
		
			# Get the corresponding plot coordinates
			start_x = max_x - branch_piece_BOT_time_bp
			end_x = max_x - branch_piece_TOP_time_bp
			
			# 2017-04-07
			if (start_x > end_x)
				{
				warning("WARNING: start_x > end_x")
				start_x = end_x
				}
			
			
			#######################################################
			# PAINT THE BRANCH WITH NO ANAGENETIC CHANGE
			#######################################################
			# Plot the segment
			statenum_1based = master_table_cladogenetic_events$sampled_states_AT_nodes[rownum]
			tmpcolor = colors_list_for_states[statenum_1based]
			lwd_tmp = rangesizes[statenum_1based]
			# For a shift to NULL, don't plot line
			if (include_null_range==TRUE && statenum_1based==1)
				{
				segments(x0=end_x, x1=end_x, y0=start_y, y1=end_y, col=tmpcolor, lwd=lwd_tmp, lty=lty, lend="butt")
				} else {
				segments(x0=start_x, x1=end_x, y0=start_y, y1=end_y, col=tmpcolor, lwd=lwd_tmp, lty=lty, lend="butt")
				}
			
			list_of_each_state = c(list_of_each_state, statenum_1based)
			lengths_in_each_state = c(lengths_in_each_state, end_x - start_x)

			
			# END if no anagenetic events
			} else {
			# ANAGENETIC EVENTS ON BRANCH
			# Get the master nodenum
			master_nodenum = master_table_cladogenetic_events$node[rownum]
			start_y = nodes_xy$y[master_nodenum]
			end_y = nodes_xy$y[master_nodenum]

			# Get the branch text and extract the history
			branch_events_txt = master_table_cladogenetic_events$anagenetic_events_txt_below_node[rownum]
			events_table_for_branch = events_txt_into_events_table(branch_events_txt)
			
			
			if (plot_clado_1desc == FALSE)
				{
				events_table_for_branch_orig = events_table_for_branch
				
				TF1 = events_table_for_branch$event_type == "d"
				TF2 = events_table_for_branch$event_type == "e"
				TF3 = events_table_for_branch$event_type == "a"
				anagenetic_TF = (TF1 + TF2 + TF3) == 1
				
				events_table_for_branch = events_table_for_branch[anagenetic_TF, ]
				
				# Skip the loop if no anagenetic events
				if (sum(anagenetic_TF) == 0)
					{
					next()
					}
				} # END if (plot_clado_1desc == FALSE)
			
		
			if (stratified == TRUE)
				{
				time_top = master_table_cladogenetic_events$time_top[rownum]
				# This is absolute time before present for the branch piece we are
				# looking at
				branch_piece_TOP_time_bp = time_top + master_table_cladogenetic_events$SUBtime_bp[rownum]
				} else {
				time_top = master_table_cladogenetic_events$time_bp[rownum]
				branch_piece_TOP_time_bp = time_top
				} # END if (stratified == TRUE)

			# Just store the branch length and recall it here
			
			if (stratified == TRUE)
				{
				brlen_in_section = as.numeric(events_table_for_branch$brlen)
				branch_piece_BOT_time_bp = branch_piece_TOP_time_bp + brlen_in_section[1]
				} else {
				brlen_in_section = master_table_cladogenetic_events$edge.length[rownum]
				branch_piece_BOT_time_bp = branch_piece_TOP_time_bp + brlen_in_section
				}
		
		
			# Go through the events along the branch
			current_time_in_subbranch = 0
			
			# I think everything is off by 1 event
			old_painting_method = FALSE
			
			# Sort the branch events by time!!
			events_table_for_branch = events_table_for_branch[order(events_table_for_branch$event_time),]
			
			if (old_painting_method == FALSE)
				{
				# The first "event" is actually the corner state
				# up to the first event
				branch_event_BOT_time_bp = branch_piece_BOT_time_bp - current_time_in_subbranch
				branch_event_TOP_time_bp = as.numeric(events_table_for_branch$abs_event_time[1])
				
				# Need special edit if the ending state is "_" (null) - 2018-01-05
				if (events_table_for_branch$new_rangetxt[1] == "_")
					{
					branch_event_TOP_time_bp = branch_piece_TOP_time_bp
					}
				
				# Get the state at the bottom of the branch
				statenum_1based = master_table_cladogenetic_events$sampled_states_AT_brbots[rownum]
				tmpcolor = colors_list_for_states[statenum_1based]
				lwd_tmp = rangesizes[statenum_1based]
			
				#######################################################
				# PAINT THE FIRST EVENT ON THE BRANCH (new painting method)
				#######################################################
				# Get the corresponding plot coordinates
				start_x = max_x - branch_event_BOT_time_bp
				end_x = max_x - branch_event_TOP_time_bp
				# For a shift to NULL, don't plot line
				if (include_null_range==TRUE && statenum_1based==1)
					{
					segments(x0=end_x, x1=end_x, y0=start_y, y1=end_y, col=tmpcolor, lwd=lwd_tmp, lty=lty, lend="butt")
					} else {
					segments(x0=start_x, x1=end_x, y0=start_y, y1=end_y, col=tmpcolor, lwd=lwd_tmp, lty=lty, lend="butt")
					}

				current_time_in_subbranch = branch_event_TOP_time_bp
				
				list_of_each_state = c(list_of_each_state, statenum_1based)
				lengths_in_each_state = c(lengths_in_each_state, end_x - start_x)
				}
			
			# Loop through the events
			for (j in 1:nrow(events_table_for_branch))
				{
				# Time above the bottom of the branch
				#event_time = branch_piece_BOT_time_bp - #as.numeric(events_table_for_branch$event_time[j])
				#as.numeric(events_table_for_branch$event_time[j])
			
				# Calculate the time_bp of event in x
				if (old_painting_method == TRUE)
					{
					branch_event_BOT_time_bp = branch_piece_BOT_time_bp - current_time_in_subbranch
					branch_event_TOP_time_bp = as.numeric(events_table_for_branch$abs_event_time[j])
					} else {
					# Event bottom
					#branch_event_BOT_time_bp = branch_piece_BOT_time_bp - current_time_in_subbranch
					branch_event_BOT_time_bp = as.numeric(events_table_for_branch$abs_event_time[j])
					# Event top
					if (j < nrow(events_table_for_branch))
						{
						branch_event_TOP_time_bp = as.numeric(events_table_for_branch$abs_event_time[j+1])
						} else {
						branch_event_TOP_time_bp = as.numeric(branch_piece_TOP_time_bp)
						}
					}
			
				#cat("\n\n")
				#cat(events_table_for_branch$nodenum_at_top_of_branch[j], ": ", branch_event_BOT_time_bp, " -- ", branch_event_TOP_time_bp, sep="")
				#cat("\n\n")
			
				
			
				# Get the corresponding plot coordinates
				start_x = max_x - branch_event_BOT_time_bp
				end_x = max_x - branch_event_TOP_time_bp
			
				# Plot the segment
				# The old method used the current rangenum
				# The new method uses the new rangenum
				if (old_painting_method == TRUE)
					{
					statenum_1based = as.numeric(events_table_for_branch$current_rangenum_1based[j])
					} else {
					statenum_1based = as.numeric(events_table_for_branch$new_rangenum_1based[j])
					}
				
				#######################################################
				# PAINT THE NEXT EVENT ON THE BRANCH (new painting method)
				#######################################################
				# Paint the segment
				tmpcolor = colors_list_for_states[statenum_1based]
				lwd_tmp = rangesizes[statenum_1based]

				# For a shift to NULL, don't plot line
				if (include_null_range==TRUE && statenum_1based==1)
					{
					segments(x0=end_x, x1=end_x, y0=start_y, y1=end_y, col=tmpcolor, lwd=lwd_tmp, lty=lty, lend="butt")
					} else {
					segments(x0=start_x, x1=end_x, y0=start_y, y1=end_y, col=tmpcolor, lwd=lwd_tmp, lty=lty, lend="butt")
					}

				
				list_of_each_state = c(list_of_each_state, statenum_1based)
				lengths_in_each_state = c(lengths_in_each_state, end_x - start_x)


				if (plot_clado_1desc_points == TRUE)
					{
					# Default is NA (don't plot)
					pointchar = NA
					
					if ((events_table_for_branch$event_type[j] == "sympatry (y)") || (events_table_for_branch$event_type == "y") )
						{
						# Character to plot
						pointchar = "y\n|"
						}
					if ((events_table_for_branch$event_type[j] == "vicariance (v)") || (events_table_for_branch$event_type == "v") )
						{
						# Character to plot
						pointchar = "v\n|"
						}
					if ((events_table_for_branch$event_type[j] == "subset (s)") || (events_table_for_branch$event_type == "s") )
						{
						# Character to plot
						pointchar = "s\n|"
						}
					if ((events_table_for_branch$event_type[j] == "founder (j)") || (events_table_for_branch$event_type == "j") )
						{
						# Character to plot
						pointchar = "j\n|"
						}
					if (events_table_for_branch$event_type[j] == "d")
						{
						# Character to plot
						pointchar = "d\n|"
						}
					if (events_table_for_branch$event_type[j] == "e")
						{
						# Character to plot
						pointchar = "e\n|"
						}
					if (events_table_for_branch$event_type[j] == "a")
						{
						# Character to plot
						pointchar = "a\n|"
						}
					
					# Plot the point
					if (!is.na(pointchar))
						{
						# Plot points
						# points(x=start_x, y=start_y, pch="+", col="black", cex=1)
						
						# Plot some text
						text(x=start_x, y=start_y, labels=pointchar, adj=c(0.5,0.15), col=tmpcolor, cex=0.85)
						#text(x=start_x, y=start_y, labels=pointchar, adj=c(0.5,0.15), col="black", cex=0.85)
						} # END if (!is.na(pointchar))
					} # END if (plot_clado_1desc_points == TRUE)

				# Update the time of the base of the next event
				current_time_in_subbranch = branch_event_BOT_time_bp - branch_event_TOP_time_bp
				} # END for (j in 1:nrow(events_table_for_branch))
		

			# Calculate the time_bp of event in x
			branch_event_BOT_time_bp = branch_event_TOP_time_bp
			branch_event_TOP_time_bp = branch_piece_TOP_time_bp

			# This shouldn't be needed anymore
			if (old_painting_method == TRUE)
				{
				# Then, after the last change event, you still HAVE to plot the last chunk
				# of branch events
				j = nrow(events_table_for_branch)

				# Time above the bottom of the branch
				#event_time = as.numeric(events_table_for_branch$event_time[j])
				#event_time = branch_event_TOP_time_bp
		
				#cat("\n\n")
				#cat(events_table_for_branch$nodenum_at_top_of_branch[j], ": ", branch_event_BOT_time_bp, " -- ", branch_event_TOP_time_bp, sep="")
				#cat("\n\n")

		
		
				# Get the corresponding plot coordinates
				start_x = max_x - branch_event_BOT_time_bp
				end_x = max_x - branch_event_TOP_time_bp
		
				#######################################################
				# PAINT the remaining part of the branch (old painting method)
				#######################################################
				# Plot the segment
				statenum_1based = as.numeric(events_table_for_branch$new_rangenum_1based[j])
				tmpcolor = colors_list_for_states[statenum_1based]
				lwd_tmp = rangesizes[statenum_1based]
				
				# For a shift to NULL, don't plot line
				if (include_null_range==TRUE && statenum_1based==1)
					{
					segments(x0=end_x, x1=end_x, y0=start_y, y1=end_y, col=tmpcolor, lwd=lwd_tmp, lty=lty, lend="butt")
					} else {
					segments(x0=start_x, x1=end_x, y0=start_y, y1=end_y, col=tmpcolor, lwd=lwd_tmp, lty=lty, lend="butt")
					}

				list_of_each_state = c(list_of_each_state, statenum_1based)
				lengths_in_each_state = c(lengths_in_each_state, end_x - start_x)
				} # END if (old_painting_method == TRUE)
			} # END if anagenetic events
		} # END for (rownum in 1:nrow(master_table_cladogenetic_events))
		  # (ENDING loop through the branches and subbranches in
		  #  master_table_cladogenetic_events)
	
	times_in_each_state = NULL
	times_in_each_state$list_of_each_state = list_of_each_state
	times_in_each_state$lengths_in_each_state = lengths_in_each_state

	return(times_in_each_state)
	} # END function



#######################################################
# Plot stochastic maps
#######################################################

plot_BSM <- function(results_object, clado_events_table, stratified, analysis_titletxt="Stochastic map", addl_params=list(), plotwhat="text", label.offset=NULL, tipcex=0.8, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, plotlegend=FALSE, legend_ncol=NULL, legend_cex=1, cornercoords_loc="manual", tr=NULL, tipranges=NULL, if_ties="takefirst", pie_tip_statecex=0.7, juststats=FALSE, xlab="Millions of years ago", root.edge=TRUE, colors_list_for_states, skiptree=FALSE, show.tip.label=TRUE, tipcol="black", dej_params_row=NULL, plot_max_age=NULL, skiplabels=FALSE, plot_stratum_lines=TRUE, include_null_range=NULL, plot_null_range=FALSE)
	{
	scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
	
	# Convert the BSM into a modified res object
	resmod = stochastic_map_states_into_res(res=results_object, master_table_cladogenetic_events=clado_events_table, stratified=stratified)

	# Plot the tree and states at nodes/corners
	# (copying everything from the inputs; mostly these should be kept on defaults)
	# (skiptree=FALSE the first time, TRUE the second time)
	plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt=analysis_titletxt, addl_params=addl_params, plotwhat=plotwhat, label.offset=label.offset, tipcex=tipcex, statecex=statecex, splitcex=splitcex, titlecex=titlecex, plotsplits=plotsplits, plotlegend=plotlegend, legend_ncol=legend_ncol, legend_cex=legend_cex, cornercoords_loc=cornercoords_loc, tr=tr, tipranges=tipranges, if_ties=if_ties, pie_tip_statecex=pie_tip_statecex, juststats=juststats, xlab=xlab, root.edge=root.edge, colors_list_for_states=colors_list_for_states, skiptree=FALSE, show.tip.label=show.tip.label, tipcol=tipcol, dej_params_row=dej_params_row, plot_max_age=plot_max_age, skiplabels=skiplabels, plot_stratum_lines=plot_stratum_lines, include_null_range=include_null_range, plot_null_range=plot_null_range)
	
	# Paint on the branch states
	paint_stochastic_map_branches(res=resmod, master_table_cladogenetic_events=clado_events_table, colors_list_for_states=colors_list_for_states, lwd=5, lty=par("lty"), root.edge=TRUE, stratified=stratified)

	# Re-plot the tree to get the states on top
	# (skiptree=TRUE this time)
	plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt=analysis_titletxt, addl_params=addl_params, plotwhat=plotwhat, label.offset=label.offset, tipcex=tipcex, statecex=statecex, splitcex=splitcex, titlecex=titlecex, plotsplits=plotsplits, plotlegend=plotlegend, legend_ncol=legend_ncol, legend_cex=legend_cex, cornercoords_loc=cornercoords_loc, tr=tr, tipranges=tipranges, if_ties=if_ties, pie_tip_statecex=pie_tip_statecex, juststats=juststats, xlab=xlab, root.edge=root.edge, colors_list_for_states=colors_list_for_states, skiptree=TRUE, show.tip.label=show.tip.label, tipcol=tipcol, dej_params_row=dej_params_row, plot_max_age=plot_max_age, skiplabels=skiplabels, plot_stratum_lines=plot_stratum_lines, include_null_range=include_null_range, plot_null_range=plot_null_range)	
	
	return(NULL)
	} # END plot_BSM


# Various problems emerge from "ladderize" in some versions
# https://www.mail-archive.com/r-sig-phylo@r-project.org/msg04176.html
ladderize_and_reorder <- function(phy, right=TRUE)
	{
	ltr = ladderize(phy, right=right)
	# MAKE FREAKING SURE that this tree has the right node order etc.
	ltr = read.tree(file="", text=write.tree(phy=ltr, file="") )
	return(ltr)
	} # END ladderize_and_reorder <- function(tr, right=TRUE)


# Returns:
# indexes_to_convert_tr2nodes_to_tr1
# NA for non-matches

# E.g. 
# tr1 = new tree
# tr2 = original tree, used the BioGeoBEARS analysis
ordernodes <- function(tr1, tr2)
	{
	tr1_table = prt(tr1, printflag=FALSE, get_tipnames=TRUE)
	tr2_table = prt(tr2, printflag=FALSE, get_tipnames=TRUE)
	indexes_to_convert_tr2nodes_to_tr1 = match(x=tr1_table$tipnames, table=tr2_table$tipnames)
	
	if (any(is.na(indexes_to_convert_tr2nodes_to_tr1)) == TRUE)
		{
		txt = "WARNING in ordernodes() or BioGeoBEARS_reorder() -- not all nodes match between tr1 and tr2. Any non-matching nodes will get NAs."
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		
		warning(txt)
		} # END if (any(is.na(indexes_to_convert_tr2nodes_to_tr1)) == TRUE)
	
	return(indexes_to_convert_tr2nodes_to_tr1)
	} # END ordernodes <- function(tr1, tr2)


# Convert BioGeoBEARS object from one tree to another
# tr1 = new tree
# tr2 = original tree, used the BioGeoBEARS analysis
BioGeoBEARS_reorder <- function(res, tr1, tr2, trfn_for_BGB_inputs=NULL)
	{
	indexes_to_convert_tr2nodes_to_tr1 = ordernodes(tr1, tr2)
	res2 = res
	
	# Vector
	res2$computed_likelihoods_at_each_node = res$computed_likelihoods_at_each_node[indexes_to_convert_tr2nodes_to_tr1]
	
	# Matrices
	res2$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS = res$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[indexes_to_convert_tr2nodes_to_tr1,]
	
	res2$condlikes_of_each_state = res$condlikes_of_each_state[indexes_to_convert_tr2nodes_to_tr1,]
	
	res2$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS = res$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[indexes_to_convert_tr2nodes_to_tr1,]
	
	res2$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS = res$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS[indexes_to_convert_tr2nodes_to_tr1,]
	
	res2$relative_probs_of_each_state_at_branch_top_AT_node_UPPASS = res$relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[indexes_to_convert_tr2nodes_to_tr1,]
	
	res2$ML_marginal_prob_each_state_at_branch_bottom_below_node = res$ML_marginal_prob_each_state_at_branch_bottom_below_node[indexes_to_convert_tr2nodes_to_tr1,]
	
	res2$ML_marginal_prob_each_state_at_branch_top_AT_node = res$ML_marginal_prob_each_state_at_branch_top_AT_node[indexes_to_convert_tr2nodes_to_tr1,]
	
	# Input the tree filename into BioGeoBEARS inputs, if desired
	if (is.null(trfn_for_BGB_inputs) == FALSE)
		{
		res2$inputs$trfn = trfn_for_BGB_inputs
		}
	
	return(res2)
	} # END BioGeoBEARS_reorder <- function(res, tr1, tr2)



run_BioGeoBEARS_reorder <- function()
		{
	#######################################################
	# Reorder Robin Becks tree so that the BioGeoBEARS plots look better
	#######################################################
	# Convert BioGeoBEARS object from one tree to another
	# tr1 = new tree
	# tr2 = original tree, used the BioGeoBEARS analysis
	wd  = "/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/z_help/Robin_Beck/"
	setwd(wd)

	# Original tree
	trfn = "PBG_TE_no_Alo_IGR_topcons_MCC_no_zero.tre"
	tr2 = read.tree(trfn)
	plot(tr2)

	# Ladderize the tree
	# GOOD LADDERIZE
	tr1 = ladderize_and_reorder(phy=tr2, right=FALSE)
	plot(tr1)

	tr1 = rotate_tips(tr=tr1, tipname1="Kulbeckia", tipname2="Paranyctoides")
	plot(tr1)

	tr1 = rotate_tips(tr=tr1, tipname1="Paranyctoides", tipname2="Sheikhdzheilia")
	plot(tr1)




	new_trfn = "tree_ladderized.newick"
	write.tree(phy=tr1, file=new_trfn)



	resfns = c("Mammals_BAYAREALIKE+J_M0_unconstrained_v1.Rdata",
	"Mammals_BAYAREALIKE_M0_unconstrained_v1.Rdata",
	"Mammals_DIVALIKE+J_M0_unconstrained_v1.Rdata",
	"Mammals_DIVALIKE_M0_unconstrained_v1.Rdata",
	"Mammals_DEC+J_M0_unconstrained_v1.Rdata",
	"Mammals_DEC_M0_unconstrained_v1.Rdata")

	for (i in 1:length(resfns))
		{
		# Convert each file
	
		# Loads to res
		load(file=resfns[i])
		res2 = BioGeoBEARS_reorder(res, tr1, tr2, trfn_for_BGB_inputs=new_trfn)
	
		new_resfn = gsub(pattern="\\.Rdata", replacement="_ladderized.Rdata", x=resfns[i])
		res = res2
		# Loads to res
		save(res, file=new_resfn)
		} # END for (i in 1:length(resfns))

	} # END run_BioGeoBEARS_reorder










#######################################################
# These functions fix graphics issues that arise in APE 5.0
# 
# They were pointed out by Liz, here:
#
# https://groups.google.com/forum/#!topic/biogeobears/gQ4bZ4U3FU8
# 
# BioGeoBEARS ›
# node_height error
# 1 post by 1 author  
# 
# Liz	
# 
# 12:25 PM (2 hours ago)
# 
# Hi all,
# 
# I'm a new user of BioGeoBEARS and I'm just trying to get the sample 
# script running.  I'm hoping the following is a simple issue....  
# 
# When I run the test script, everything runs fine until this command:
# 
# res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
# 
# I get the following traceback error in R Studio:
# 
#  Error in .C("node_height", as.integer(Ntip), as.integer(Nnode), 
#               as.integer(edge[,  : 
#   Incorrect number of arguments (6), expecting 4 for 'node_height' 
# 6. .nodeHeight(Ntip, Nnode, z$edge, Nedge, yy) at plot_phylo3_nodecoords.R#156
# 5. plot_phylo3_nodecoords(tr, plot = FALSE, root.edge = root.edge) at <text>#1
# 4. eval(parse(text = cmdstr)) 
# 3. eval(parse(text = cmdstr)) at BioGeoBEARS_plots_v1.R#1871
# 2. node_coords(tr, tmplocation = cornercoords_loc, root.edge = root.edge) at BioGeoBEARS_plots_v1.R#457
# 1. plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params = list("j"), 
#  plotwhat = "text", label.offset = 0.45, tipcex = 0.7, statecex = 0.7, 
#  splitcex = 0.6, titlecex = 0.8, plotsplits = TRUE, 
#  cornercoords_loc=scriptdir, 
#  include_null_range = TRUE, tr = tr, tipranges = tipranges) 
# 
# Please note I am using ape version 5.0.
# 
# Thanks for your help!
# Liz
#######################################################


# Access par (graphics parameters) without opening a &%$@ plot!
# https://stackoverflow.com/questions/20363266/how-can-i-suppress-the-creation-of-a-plot-while-calling-a-function-in-r
par_invisible <- function(parname)
	{
	setup='
	pin1 = par_invisible(parname="pin")
	'

	ff <- tempfile()
	png(filename=ff)
	tmp = do.call(par, args=list(parname))
	res = tmp[1]
	dev.off()
	unlink(ff)
	return(res)
	}



plot_phylo3_nodecoords_APE5 <- function (x, type = "phylogram", use.edge.length = TRUE, node.pos = NULL, 
    show.tip.label = TRUE, show.node.label = FALSE, edge.color = "black", 
    edge.width = 1, edge.lty = 1, font = 3, 
    adj = NULL, srt = 0, no.margin = FALSE, root.edge = FALSE, 
    label.offset = 0, underscore = FALSE, x.lim = NULL, y.lim = NULL, 
    direction = "rightwards", lab4ut = NULL, tip.color = "black", 
    plot = TRUE, rotate.tree = 0, open.angle = 0, node.depth = 1, 
    align.tip.label = FALSE, cex_original=FALSE, ...) 
{
    
    # Original behavior in APE 5.0 plot.phylo
    if (cex_original == TRUE)
    	{
    	cex = par("cex")
    	} else {
    	cex = par_invisible(parname="cex")
    	}
    
    Ntip <- length(x$tip.label)
    if (Ntip < 2) {
        warning("found less than 2 tips in the tree")
        return(NULL)
    }
    .nodeHeight <- function(edge, Nedge, yy) .C(node_height, 
        as.integer(edge[, 1]), as.integer(edge[, 2]), as.integer(Nedge), 
        as.double(yy))[[4]]
    .nodeDepth <- function(Ntip, Nnode, edge, Nedge, node.depth) .C(node_depth, 
        as.integer(Ntip), as.integer(edge[, 1]), as.integer(edge[, 
            2]), as.integer(Nedge), double(Ntip + Nnode), as.integer(node.depth))[[5]]
    .nodeDepthEdgelength <- function(Ntip, Nnode, edge, Nedge, 
        edge.length) .C(node_depth_edgelength, as.integer(edge[, 
        1]), as.integer(edge[, 2]), as.integer(Nedge), as.double(edge.length), 
        double(Ntip + Nnode))[[5]]
    Nedge <- dim(x$edge)[1]
    Nnode <- x$Nnode
    if (any(x$edge < 1) || any(x$edge > Ntip + Nnode)) 
        stop("tree badly conformed; cannot plot. Check the edge matrix.")
    ROOT <- Ntip + 1
    type <- match.arg(type, c("phylogram", "cladogram", "fan", 
        "unrooted", "radial"))
    direction <- match.arg(direction, c("rightwards", "leftwards", 
        "upwards", "downwards"))
    if (is.null(x$edge.length)) {
        use.edge.length <- FALSE
    } else {
        if (use.edge.length && type != "radial") {
            tmp <- sum(is.na(x$edge.length))
            if (tmp) {
                warning(paste(tmp, "branch length(s) NA(s): branch lengths ignored in the plot"))
                use.edge.length <- FALSE
            }
        }
    }
    if (is.numeric(align.tip.label)) {
        align.tip.label.lty <- align.tip.label
        align.tip.label <- TRUE
    } else {
        if (align.tip.label) 
            align.tip.label.lty <- 3
    }
    if (align.tip.label) {
        if (type %in% c("unrooted", "radial") || !use.edge.length || 
            is.ultrametric(x)) 
            align.tip.label <- FALSE
    }
    if (type %in% c("unrooted", "radial") || !use.edge.length || 
        is.null(x$root.edge) || !x$root.edge) 
        root.edge <- FALSE
    phyloORclado <- type %in% c("phylogram", "cladogram")
    horizontal <- direction %in% c("rightwards", "leftwards")
    xe <- x$edge
    if (phyloORclado) {
        phyOrder <- attr(x, "order")
        if (is.null(phyOrder) || phyOrder != "cladewise") {
            x <- reorder(x)
            if (!identical(x$edge, xe)) {
                ereorder <- match(x$edge[, 2], xe[, 2])
                if (length(edge.color) > 1) {
                  edge.color <- rep(edge.color, length.out = Nedge)
                  edge.color <- edge.color[ereorder]
                }
                if (length(edge.width) > 1) {
                  edge.width <- rep(edge.width, length.out = Nedge)
                  edge.width <- edge.width[ereorder]
                }
                if (length(edge.lty) > 1) {
                  edge.lty <- rep(edge.lty, length.out = Nedge)
                  edge.lty <- edge.lty[ereorder]
                }
            }
        }
        yy <- numeric(Ntip + Nnode)
        TIPS <- x$edge[x$edge[, 2] <= Ntip, 2]
        yy[TIPS] <- 1:Ntip
    }
    z <- reorder(x, order = "postorder")
    if (phyloORclado) {
        if (is.null(node.pos)) 
            node.pos <- if (type == "cladogram" && !use.edge.length) 
                2
            else 1
        if (node.pos == 1) 
            yy <- .nodeHeight(z$edge, Nedge, yy)
        else {
            ans <- .C(node_height_clado, as.integer(Ntip), as.integer(z$edge[, 
                1]), as.integer(z$edge[, 2]), as.integer(Nedge), 
                double(Ntip + Nnode), as.double(yy))
            xx <- ans[[5]] - 1
            yy <- ans[[6]]
        }
        if (!use.edge.length) {
            if (node.pos != 2) 
                xx <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, 
                  node.depth) - 1
            xx <- max(xx) - xx
        } else {
            xx <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, Nedge, 
                z$edge.length)
        }
    } else {
        twopi <- 2 * pi
        rotate.tree <- twopi * rotate.tree/360
        if (type != "unrooted") {
            TIPS <- x$edge[which(x$edge[, 2] <= Ntip), 2]
            xx <- seq(0, twopi * (1 - 1/Ntip) - twopi * open.angle/360, 
                length.out = Ntip)
            theta <- double(Ntip)
            theta[TIPS] <- xx
            theta <- c(theta, numeric(Nnode))
        }
        switch(type, fan = {
            theta <- .nodeHeight(z$edge, Nedge, theta)
            if (use.edge.length) {
                r <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, 
                  Nedge, z$edge.length)
            } else {
                r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
                r <- 1/r
            }
            theta <- theta + rotate.tree
            if (root.edge) r <- r + x$root.edge
            xx <- r * cos(theta)
            yy <- r * sin(theta)
        }, unrooted = {
            nb.sp <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
            XY <- if (use.edge.length) unrooted.xy(Ntip, Nnode, 
                z$edge, z$edge.length, nb.sp, rotate.tree) else unrooted.xy(Ntip, 
                Nnode, z$edge, rep(1, Nedge), nb.sp, rotate.tree)
            xx <- XY$M[, 1] - min(XY$M[, 1])
            yy <- XY$M[, 2] - min(XY$M[, 2])
        }, radial = {
            r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
            r[r == 1] <- 0
            r <- 1 - r/Ntip
            theta <- .nodeHeight(z$edge, Nedge, theta) + rotate.tree
            xx <- r * cos(theta)
            yy <- r * sin(theta)
        })
    }
    
    
    
    if (phyloORclado) {
        if (!horizontal) {
            tmp <- yy
            yy <- xx
            xx <- tmp - min(tmp) + 1
        }
        if (root.edge) {
            if (direction == "rightwards") 
                xx <- xx + x$root.edge
            if (direction == "upwards") 
                yy <- yy + x$root.edge
        }
    }
    if (no.margin)
    	{
			# Original behavior in APE 5.0 plot.phylo
			if (cex_original == TRUE)
				{
				par(mai = rep(0, 4))
				}
    	}

    if (show.tip.label) 
        nchar.tip.label <- nchar(x$tip.label)
    max.yy <- max(yy)
    getLimit <- function(x, lab, sin, cex, cex_original=TRUE) {
        if (cex_original == TRUE)
        	{
	        s <- strwidth(lab, "inches", cex = cex)
	        } else {
	        ff <- tempfile()
					png(filename=ff)
					s = strwidth(lab, "inches", cex = cex)
					dev.off()
					unlink(ff)
	        }
        if (any(s > sin)) 
            return(1.5 * max(x))
        Limit <- 0
        while (any(x > Limit)) {
            i <- which.max(x)
            alp <- x[i]/(sin - s[i])
            Limit <- x[i] + alp * s[i]
            x <- x + alp * s
        }
        Limit
    }
    
    
    if (is.null(x.lim)) {
        if (phyloORclado) {
            if (horizontal) {
                xx.tips <- xx[1:Ntip]
                if (show.tip.label) {
									# Original behavior in APE 5.0 plot.phylo
									if (cex_original == TRUE)
										{
	                  pin1 <- par("pin")[1]
	                  } else {
										res = par_invisible(parname="pin")
										pin1 = res[1]
	                  }
                  tmp <- getLimit(xx.tips, x$tip.label, pin1, 
                    cex, cex_original=cex_original)
                  tmp <- tmp + label.offset
                }
                else tmp <- max(xx.tips)
                x.lim <- c(0, tmp)
            }
            else x.lim <- c(1, Ntip)
        } else switch(type, fan = {
            if (show.tip.label) {
                offset <- max(nchar.tip.label * 0.018 * max.yy * 
                  cex)
                x.lim <- range(xx) + c(-offset, offset)
            } else x.lim <- range(xx)
        }, unrooted = {
            if (show.tip.label) {
                offset <- max(nchar.tip.label * 0.018 * max.yy * 
                  cex)
                x.lim <- c(0 - offset, max(xx) + offset)
            } else x.lim <- c(0, max(xx))
        }, radial = {
            if (show.tip.label) {
                offset <- max(nchar.tip.label * 0.03 * cex)
                x.lim <- c(-1 - offset, 1 + offset)
            } else x.lim <- c(-1, 1)
        })
    } else if (length(x.lim) == 1) {
        x.lim <- c(0, x.lim)
        if (phyloORclado && !horizontal) 
            x.lim[1] <- 1
        if (type %in% c("fan", "unrooted") && show.tip.label) 
            x.lim[1] <- -max(nchar.tip.label * 0.018 * max.yy * 
                cex)
        if (type == "radial") 
            x.lim[1] <- if (show.tip.label) 
                -1 - max(nchar.tip.label * 0.03 * cex)
            else -1
    }
    if (phyloORclado && direction == "leftwards") 
        xx <- x.lim[2] - xx
    if (is.null(y.lim)) {
        if (phyloORclado) {
            if (horizontal) 
                y.lim <- c(1, Ntip)
            else {
								# Original behavior in APE 5.0 plot.phylo
								if (cex_original == TRUE)
									{
									pin2 <- par("pin")[2]
									} else {
									res = par_invisible(parname="pin")
									pin2 = res[2]
									}

                yy.tips <- yy[1:Ntip]
                if (show.tip.label) {
                  tmp <- getLimit(yy.tips, x$tip.label, pin2, 
                    cex, cex_original=cex_original)
                  tmp <- tmp + label.offset
                }
                else tmp <- max(yy.tips)
                y.lim <- c(0, tmp)
            }
        } else switch(type, fan = {
            if (show.tip.label) {
                offset <- max(nchar.tip.label * 0.018 * max.yy * 
                  cex)
                y.lim <- c(min(yy) - offset, max.yy + offset)
            } else y.lim <- c(min(yy), max.yy)
        }, unrooted = {
            if (show.tip.label) {
                offset <- max(nchar.tip.label * 0.018 * max.yy * 
                  cex)
                y.lim <- c(0 - offset, max.yy + offset)
            } else y.lim <- c(0, max.yy)
        }, radial = {
            if (show.tip.label) {
                offset <- max(nchar.tip.label * 0.03 * cex)
                y.lim <- c(-1 - offset, 1 + offset)
            } else y.lim <- c(-1, 1)
        })
    } else if (length(y.lim) == 1) {
        y.lim <- c(0, y.lim)
        if (phyloORclado && horizontal) 
            y.lim[1] <- 1
        if (type %in% c("fan", "unrooted") && show.tip.label) 
            y.lim[1] <- -max(nchar.tip.label * 0.018 * max.yy * 
                cex)
        if (type == "radial") 
            y.lim[1] <- if (show.tip.label) 
                -1 - max(nchar.tip.label * 0.018 * max.yy * cex)
            else -1
    }
    if (phyloORclado && direction == "downwards") 
        yy <- y.lim[2] - yy
    if (phyloORclado && root.edge) {
        if (direction == "leftwards") 
            x.lim[2] <- x.lim[2] + x$root.edge
        if (direction == "downwards") 
            y.lim[2] <- y.lim[2] + x$root.edge
    }
    asp <- if (type %in% c("fan", "radial", "unrooted")) 
        1
    else NA
    
    # 2017-11-02_NJM add for plot_phylo3_nodecoords_APE5
    if (plot)
    	{
	    plot.default(0, type = "n", xlim = x.lim, ylim = y.lim, xlab = "", 
        ylab = "", axes = FALSE, asp = asp, ...)
      }
    if (plot) {
        if (is.null(adj)) 
            adj <- if (phyloORclado && direction == "leftwards") 
                1
            else 0
        if (phyloORclado && show.tip.label) {
            MAXSTRING <- max(strwidth(x$tip.label, cex = cex))
            loy <- 0
            if (direction == "rightwards") {
                lox <- label.offset + MAXSTRING * 1.05 * adj
            }
            if (direction == "leftwards") {
                lox <- -label.offset - MAXSTRING * 1.05 * (1 - 
                  adj)
            }
            if (!horizontal) {
                # Original behavior in APE 5.0 plot.phylo
								if (cex_original == TRUE)
									{
									psr <- par("usr")
									} else {
									psr = par_invisible(parname="usr")
									}
                MAXSTRING <- MAXSTRING * 1.09 * (psr[4] - psr[3])/(psr[2] - 
                  psr[1])
                loy <- label.offset + MAXSTRING * 1.05 * adj
                lox <- 0
                srt <- 90 + srt
                if (direction == "downwards") {
                  loy <- -loy
                  srt <- 180 + srt
                }
            }
        }
        if (type == "phylogram") {
            phylogram.plot(x$edge, Ntip, Nnode, xx, yy, horizontal, 
                edge.color, edge.width, edge.lty)
        } else {
            if (type == "fan") {
                ereorder <- match(z$edge[, 2], x$edge[, 2])
                if (length(edge.color) > 1) {
                  edge.color <- rep(edge.color, length.out = Nedge)
                  edge.color <- edge.color[ereorder]
                }
                if (length(edge.width) > 1) {
                  edge.width <- rep(edge.width, length.out = Nedge)
                  edge.width <- edge.width[ereorder]
                }
                if (length(edge.lty) > 1) {
                  edge.lty <- rep(edge.lty, length.out = Nedge)
                  edge.lty <- edge.lty[ereorder]
                }
                circular.plot(z$edge, Ntip, Nnode, xx, yy, theta, 
                  r, edge.color, edge.width, edge.lty)
            }
            else cladogram.plot(x$edge, xx, yy, edge.color, edge.width, 
                edge.lty)
        }
        if (root.edge) {
            rootcol <- if (length(edge.color) == 1) 
                edge.color
            else "black"
            rootw <- if (length(edge.width) == 1) 
                edge.width
            else 1
            rootlty <- if (length(edge.lty) == 1) 
                edge.lty
            else 1
            if (type == "fan") {
                tmp <- polar2rect(x$root.edge, theta[ROOT])
                segments(0, 0, tmp$x, tmp$y, col = rootcol, lwd = rootw, 
                  lty = rootlty)
            }
            else {
                switch(direction, rightwards = segments(0, yy[ROOT], 
                  x$root.edge, yy[ROOT], col = rootcol, lwd = rootw, 
                  lty = rootlty), leftwards = segments(xx[ROOT], 
                  yy[ROOT], xx[ROOT] + x$root.edge, yy[ROOT], 
                  col = rootcol, lwd = rootw, lty = rootlty), 
                  upwards = segments(xx[ROOT], 0, xx[ROOT], x$root.edge, 
                    col = rootcol, lwd = rootw, lty = rootlty), 
                  downwards = segments(xx[ROOT], yy[ROOT], xx[ROOT], 
                    yy[ROOT] + x$root.edge, col = rootcol, lwd = rootw, 
                    lty = rootlty))
            }
        }
        if (show.tip.label) {
            if (is.expression(x$tip.label)) 
                underscore <- TRUE
            if (!underscore) 
                x$tip.label <- gsub("_", " ", x$tip.label)
            if (phyloORclado) {
                if (align.tip.label) {
                  xx.tmp <- switch(direction, rightwards = max(xx[1:Ntip]), 
                    leftwards = min(xx[1:Ntip]), upwards = xx[1:Ntip], 
                    downwards = xx[1:Ntip])
                  yy.tmp <- switch(direction, rightwards = yy[1:Ntip], 
                    leftwards = yy[1:Ntip], upwards = max(yy[1:Ntip]), 
                    downwards = min(yy[1:Ntip]))
                  segments(xx[1:Ntip], yy[1:Ntip], xx.tmp, yy.tmp, 
                    lty = align.tip.label.lty)
                }
                else {
                  xx.tmp <- xx[1:Ntip]
                  yy.tmp <- yy[1:Ntip]
                }
                text(xx.tmp + lox, yy.tmp + loy, x$tip.label, 
                  adj = adj, font = font, srt = srt, cex = cex, 
                  col = tip.color)
            }
            else {
                angle <- if (type == "unrooted") 
                  XY$axe
                else atan2(yy[1:Ntip], xx[1:Ntip])
                lab4ut <- if (is.null(lab4ut)) {
                  if (type == "unrooted") 
                    "horizontal"
                  else "axial"
                }
                else match.arg(lab4ut, c("horizontal", "axial"))
                xx.tips <- xx[1:Ntip]
                yy.tips <- yy[1:Ntip]
                if (label.offset) {
                  xx.tips <- xx.tips + label.offset * cos(angle)
                  yy.tips <- yy.tips + label.offset * sin(angle)
                }
                if (lab4ut == "horizontal") {
                  y.adj <- x.adj <- numeric(Ntip)
                  sel <- abs(angle) > 0.75 * pi
                  x.adj[sel] <- -strwidth(x$tip.label)[sel] * 
                    1.05
                  sel <- abs(angle) > pi/4 & abs(angle) < 0.75 * 
                    pi
                  x.adj[sel] <- -strwidth(x$tip.label)[sel] * 
                    (2 * abs(angle)[sel]/pi - 0.5)
                  sel <- angle > pi/4 & angle < 0.75 * pi
                  y.adj[sel] <- strheight(x$tip.label)[sel]/2
                  sel <- angle < -pi/4 & angle > -0.75 * pi
                  y.adj[sel] <- -strheight(x$tip.label)[sel] * 
                    0.75
                  text(xx.tips + x.adj * cex, yy.tips + y.adj * 
                    cex, x$tip.label, adj = c(adj, 0), font = font, 
                    srt = srt, cex = cex, col = tip.color)
                }
                else {
                  if (align.tip.label) {
                    POL <- rect2polar(xx.tips, yy.tips)
                    POL$r[] <- max(POL$r)
                    REC <- polar2rect(POL$r, POL$angle)
                    xx.tips <- REC$x
                    yy.tips <- REC$y
                    segments(xx[1:Ntip], yy[1:Ntip], xx.tips, 
                      yy.tips, lty = align.tip.label.lty)
                  }
                  if (type == "unrooted") {
                    adj <- abs(angle) > pi/2
                    angle <- angle * 180/pi
                    angle[adj] <- angle[adj] - 180
                    adj <- as.numeric(adj)
                  }
                  else {
                    s <- xx.tips < 0
                    angle <- angle * 180/pi
                    angle[s] <- angle[s] + 180
                    adj <- as.numeric(s)
                  }
                  font <- rep(font, length.out = Ntip)
                  tip.color <- rep(tip.color, length.out = Ntip)
                  cex <- rep(cex, length.out = Ntip)
                  for (i in 1:Ntip) text(xx.tips[i], yy.tips[i], 
                    x$tip.label[i], font = font[i], cex = cex[i], 
                    srt = angle[i], adj = adj[i], col = tip.color[i])
                }
            }
        }
        if (show.node.label) 
            text(xx[ROOT:length(xx)] + label.offset, yy[ROOT:length(yy)], 
                x$node.label, adj = adj, font = font, srt = srt, 
                cex = cex)
    }
    # 2017-11-02_NJM add for plot_phylo3_nodecoords_APE5
    # added: , edge = xe, xx = xx, yy = yy
    L <- list(type = type, use.edge.length = use.edge.length, 
        node.pos = node.pos, node.depth = node.depth, show.tip.label = show.tip.label, 
        show.node.label = show.node.label, font = font, cex = cex, 
        adj = adj, srt = srt, no.margin = no.margin, label.offset = label.offset, 
        x.lim = x.lim, y.lim = y.lim, direction = direction, 
        tip.color = tip.color, Ntip = Ntip, Nnode = Nnode, root.time = x$root.time, 
        align.tip.label = align.tip.label, edge = xe, xx = xx, yy = yy)
    assign("last_plot.phylo", c(L, list(edge = xe, xx = xx, yy = yy)), 
        envir = .PlotPhyloEnv)
    invisible(L)
}









