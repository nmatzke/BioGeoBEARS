

#######################################################
# Make nice plots
#######################################################


#######################################################
# plot_BioGeoBEARS_results
#######################################################
#' Plot the results of a BioGeoBEARS run
#' 
#' This function plots on a tree the highest-probability ancestral states (ranges), splits if desired (these are the ranges/states just after cladogenesis, and are 
#' plotted on the corners of a tree), and/or pie charts at nodes.  A legend tying the relationship between colors and states/ranges is also optionally plotted.
#'
#' 
#' The legend is plotted on a separate plot, as it is very difficult to predict whether or not there will be space on any given tree plot.  The utility of the legend 
#' is also debatable, as \code{plot_BioGeoBEARS_results} plots the colors and state/range names directly onto the plot.  Any legend will get unwieldy above perhaps 32
#' states, which is just 5 areas with no constraints (see \code{\link[cladoRcpp]{numstates_from_numareas}}, or type \code{numstates_from_numareas(numareas5, maxareas5, include_null_range=TRUE)}.
#' 
#' Note that this assumes
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
#' @export
#' @seealso \code{\link{get_leftright_nodes_matrix_from_results}}, \code{\link{corner_coords}}, \code{\link[ape]{plot.phylo}}, \code{\link[ape]{plot.phylo}}, \code{\link[ape]{tiplabels}}, \code{\link[graphics]{legend}}, \code{\link[base]{floor}}, \code{\link[base]{ceiling}}, \code{\link[base]{floor}}, \code{\link[cladoRcpp]{numstates_from_numareas}}, \code{\link[base]{system.file}}, \code{\link[base]{list.files}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @examples
#' test=1
#' 
plot_BioGeoBEARS_results <- function(results_object, analysis_titletxt=NULL, addl_params=list(), plotwhat="text", label.offset=NULL, tipcex=0.8, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, plotlegend=FALSE, legend_ncol=NULL, legend_cex=1, cornercoords_loc="manual", include_null_range=TRUE, tr=NULL, tipranges=NULL)
	{
	
	junk='
	scriptdir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/a_scripts/"
	plot_BioGeoBEARS_results(results_object, analysis_titletxt=NULL, addl_params=list(), plotwhat="text", label.offset=NULL, tipcex=0.8, statecex=0.8, splitcex=0.8, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=NULL, tipranges=NULL)
	' # endjunk
	
	#######################################################
	# Plot ancestral states - DEC
	#######################################################


	# Setup
	#results_object = resDEC_strat
	BioGeoBEARS_run_object = results_object$inputs
	
	# Read the tree from file, if needed
	if (is.null(tr))
		{
		tr = read.tree(BioGeoBEARS_run_object$trfn)
		}
	tr2 = reorder(tr, "pruningwise")
	
	# Basic tree info
	tips = 1:length(tr2$tip.label)
	nodes = (length(tr2$tip.label)+1):(length(tr2$tip.label)+tr$Nnode)


	
	# Read the tipranges from file, if needed.
	if (is.null(tipranges))
		{
		tipranges = getranges_from_LagrangePHYLIP(BioGeoBEARS_run_object$geogfn)
		}
	
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
		numstates = numstates_from_numareas(numareas=length(areas), maxareas=max_range_size, include_null_range=include_null_range)
		numstates
		states_list = areas_list_to_states_list_new(areas, maxareas=max_range_size, include_null_range = include_null_range)
		#states_list
		states_list_0based_index = rcpp_areas_list_to_states_list(areas, maxareas=max_range_size, include_null_range = include_null_range)
		#states_list_0based_index
		} else {
		states_list_0based_index = results_object$inputs$states_list
		states_list = states_list_indexes_to_areastxt(states_list=states_list_0based_index, areanames=areas, counting_base=0, concat=FALSE, sep="")
		}


	# calculate the ML marginal probs of states at the base of each branch
	# above each split (local, non-joint probs, global model)
	results_object = get_MLsplitprobs_from_results(results_object)
	names(results_object)

	LnL = round(as.numeric(results_object$optim_result$fvalues), digits=1)
	LnL

	params_table = results_object$outputs@params_table
	params_table

	get_perEvent_probs(params_table)




	# PLOT TITLE
	# What should be on the plot title
	
	params_free_TF = params_table$type == "free"
	params_free = (rownames(params_table))[params_free_TF]
	numparams = sum(params_free_TF)
	
	# Add additional user-specified parameters, if desired
	if (length(addl_params) > 0)
		{
		params_free = c(params_free, unlist(addl_params))
		params_free = unique(params_free)
		}
	
	# Write the string of FREE parameters
	paramstrs = rep("", length(params_free))
	param_names = NULL
	param_ests = NULL
	for (i in 1:length(params_free))
		{
		param_name = params_free[i]
		param_est = params_table[param_name,"est"]
		param_print = round(param_est, digits=3)
		paramstrs[i] = paste(param_name, "=", param_print, "; ", sep="")
		
		# Store for output
		param_names = c(param_names, param_name)
		param_ests = c(param_ests, param_est)
		}
	paramstrs = c(paramstrs, "LnL=", LnL)
	paramstr = paste0(paramstrs, collapse="")
	paramstr
	
	# Store for output
	param_names = c("LnL", "numparams", param_names)
	param_ests = c(LnL, numparams, param_ests)
	
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
	
	analysis_titletxt = paste(analysis_titletxt, "\n", "anstates: global optim, ", max_range_size, " areas max. ", paramstr, sep="")
	analysis_titletxt




	#######################################################
	# Get the marginal probs of the splits (global ML probs, not local)
	# (These are marginal, rather than joint probs; but not local optima)
	#######################################################
	leftright_nodes_matrix = get_leftright_nodes_matrix_from_results(tr2, results_object)

	# This gets you the prob. of each state at the left base above each node, and
	# at the right base above each node
	marprobs = results_object$ML_marginal_prob_each_state_at_branch_bottom_below_node
	left_ML_marginals_by_node = marprobs[leftright_nodes_matrix[, 2], ]
	right_ML_marginals_by_node = marprobs[leftright_nodes_matrix[, 1], ]
	right_ML_marginals_by_node




	#######################################################
	# Extract the outputs ancestral states at nodes, and plot!
	#######################################################
	relprobs_matrix = results_object$ML_marginal_prob_each_state_at_branch_top_AT_node
	relprobs_matrix_for_internal_states = relprobs_matrix[nodes,]	# subset to just internal nodes
	relprobs_matrix

	statenames = areas_list_to_states_list_new(areas, maxareas=max_range_size, include_null_range=include_null_range, split_ABC=FALSE)
	statenames


	MLprobs = get_ML_probs(relprobs_matrix)
	MLprobs
	MLstates = get_ML_states_from_relprobs(relprobs_matrix, statenames, returnwhat="states", if_ties="takefirst")


	# Set up colors
	colors_matrix = get_colors_for_numareas(length(areas))
	colors_list_for_states = mix_colors_for_states(colors_matrix, states_list_0based_index, exclude_null=(include_null_range==FALSE))
	colors_list_for_states

	possible_ranges_list_txt = areas_list_to_states_list_new(areas,  maxareas=max_range_size, split_ABC=FALSE, include_null_range=include_null_range)
	cols_byNode = rangestxt_to_colors(possible_ranges_list_txt, colors_list_for_states, MLstates)

	# Legend, if desired
	if (plotlegend == TRUE)
		{
		colors_legend(possible_ranges_list_txt, colors_list_for_states, legend_ncol=legend_ncol, legend_cex=legend_cex)
		}

	# Plot to screen
	if (is.null(label.offset))
		{
		label.offset = 0.007 * get_max_height_tree(tr)
		}
	
	xlims = c(0, 1.25*get_max_height_tree(tr))
	
	plot(tr2, x.lim=xlims, show.tip.label=FALSE, label.offset=label.offset, cex=tipcex, no.margin=FALSE)
	tiplabels_to_plot = sapply(X=tr2$tip.label, FUN=substr, start=1, stop=30)
	tiplabels(text=tiplabels_to_plot, tip=tips, cex=tipcex, adj=0, bg="white", frame="n", pos=4, offset=label.offset)	# pos=4 means labels go to the right of coords
	
	axisPhylo()
	mtext(text="Millions of years ago", side=1, line=2)
	
	# Add states / piecharts
	if (plotwhat == "text")
		{
		nodelabels(text=MLstates[nodes], node=nodes, bg=cols_byNode[nodes], cex=statecex)
		}
	if (plotwhat == "pie")
		{
		nodelabels(pie=relprobs_matrix_for_internal_states, node=nodes, piecol=colors_list_for_states, cex=statecex)
		}
	tiplabels(text=MLstates[tips], tip=tips, bg=cols_byNode[tips], cex=statecex)	# plot tiplabels either way
	title(analysis_titletxt)


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
		coords = corner_coords(tr, tmplocation=cornercoords_loc)
		coords

		# LEFT SPLITS
		relprobs_matrix = left_ML_marginals_by_node
		
		if (plotwhat == "text")
			{
			MLprobs = get_ML_probs(relprobs_matrix)
			MLprobs
			MLstates = get_ML_states_from_relprobs(relprobs_matrix, statenames, returnwhat="states", if_ties="takefirst")
			MLstates
			length(MLstates)
		
			# Set up colors
			possible_ranges_list_txt = areas_list_to_states_list_new(areas,  maxareas=max_range_size, split_ABC=FALSE, include_null_range=include_null_range)
			cols_byNode = rangestxt_to_colors(possible_ranges_list_txt, colors_list_for_states, MLstates)
			cornerlabels(text=MLstates, coords=coords$leftcorns, bg=cols_byNode, cex=splitcex)
			}
		
		if (plotwhat == "pie")
			{
			cornerpies(pievals=relprobs_matrix, coords$leftcorns, piecol=colors_list_for_states, cex=splitcex)
			}



		# RIGHT SPLITS
		relprobs_matrix = right_ML_marginals_by_node

		if (plotwhat == "text")
			{
			MLprobs = get_ML_probs(relprobs_matrix)
			MLprobs
			MLstates = get_ML_states_from_relprobs(relprobs_matrix, statenames, returnwhat="states", if_ties="takefirst")
			MLstates
			length(MLstates)

			# Set up colors
			possible_ranges_list_txt = areas_list_to_states_list_new(areas,  maxareas=max_range_size, split_ABC=FALSE, include_null_range=include_null_range)
			cols_byNode = rangestxt_to_colors(possible_ranges_list_txt, colors_list_for_states, MLstates)

			cornerlabels(text=MLstates, coords=coords$rightcorns, bg=cols_byNode, cex=splitcex)
			}

		if (plotwhat == "pie")
			{
			cornerpies(pievals=relprobs_matrix, coords$rightcorns, piecol=colors_list_for_states, cex=splitcex)			
			}
		}



	# Handy summary outputs
	param_ests = matrix(data=param_ests, nrow=1)
	param_ests = adf2(param_ests)
	names(param_ests) = param_names
	
	param_ests = dfnums_to_numeric(param_ests)
	
	return(param_ests)
	}






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
#' @note Go BEARS!
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
		area_colors = c("black", "white")
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
#' @param exclude_null If TRUE, null ranges are excluded (however coded). Default TRUE.
#' @return \code{colors_list_for_states} The colors for the ML states
#' @export
#' @seealso \code{\link[stats]{optim}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' 
mix_colors_for_states <- function(colors_matrix, states_list_0based_index, exclude_null=TRUE)
	{
	if (exclude_null == TRUE)
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
#' @note Go BEARS!
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
#' @return Nothing
#' @export
#' @seealso \code{\link[graphics]{legend}}, \code{\link[base]{floor}}, \code{\link[base]{ceiling}}, \code{\link[base]{floor}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' 
colors_legend <- function(possible_ranges_list_txt, colors_list_for_states, legend_ncol=NULL, legend_cex=1)
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
	plot(1:10, 1:10, pch=".", col="white", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
	#lines(1:4,4:1, col="blue") 
	#legend("top", leg=c("a","b"),col=c("black","blue"), fill=TRUE) 
	legend("top", legend=possible_ranges_list_txt, fill=colors_list_for_states, ncol=legend_ncol, title="Legend", cex=legend_cex)#, fill=TRUE) 

	
	}










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
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
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
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
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
#' @return \code{statesColors_table} A table giving the colors for each state.
#' @export
#' @seealso \code{\link{get_lagrange_nodenums}}, \code{\link{LGpy_splits_fn_to_table}}, \code{\link{LGcpp_splits_fn_to_table}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @examples
#' test=1
#' 
get_statesColors_table <- function(areanames=c("K","O","M","H"))
	{
	# Make the color matrix for the individual areas
	colors_matrix = get_colors_for_numareas(numareas=length(areanames), use_rainbow=FALSE)
	colors_matrix
	
	# Get the states
	states_list_0based_index = rcpp_areas_list_to_states_list(areas=areanames, include_null_range=FALSE)
	
	colors_list_for_states = mix_colors_for_states(colors_matrix, states_list_0based_index, exclude_null=TRUE)
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
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @examples
#' test=1
#' 
map_LG_MLsplits_to_tree_corners <- function(MLsplits, tr, tipranges, removechar=NULL, type="C++", statesColors_table="default", bgcol="green3", areanames="default", newplot=TRUE, ...)
	{
	defaults='
	MLsplits=MLsplits_LGpy
	type="python"
	'
	# Get corner coordinates
	corners_list = corner_coords(tr)
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
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
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
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
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
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
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
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
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
#' @note Go BEARS!
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
#' @return \code{corners_list} 
#' @export
#' @seealso \code{\link[ape]{phylo}}, \code{\link{get_nodenums}}
#' @note Go BEARS!
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
corner_coords <- function(tr, coords_fun="plot_phylo3_nodecoords", tmplocation="manual")
	{
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
	cmdstr = paste("trcoords = ", coords_fun, "(tr, plot=FALSE)", sep="")
	eval(parse(text=cmdstr))
	
	# X and Y coords for nodes, 1-37
	nodecoords = cbind(trcoords$xx, trcoords$yy)
	
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
#' @note Go BEARS!
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
	TF = (rowSums(dfnums_to_numeric(tipranges@df))) > max_range_size
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
		states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
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
	
	# Apply element-wise
	dispersal_multipliers_matrix = dispersal_multipliers_matrix * manual_dispersal_multipliers_matrix

	#######################################################
	# multiply parameter d by dispersal_multipliers_matrix
	#######################################################
	dmat = dispersal_multipliers_matrix * matrix(d, nrow=length(areas), ncol=length(areas))
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
	val_txt = c("d", "e", "a", "x", "u", "b")
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
	val_txt = c("d", "e", "a", "x", "u")
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
#' @note Go BEARS!
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
		states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
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
	
	# Apply element-wise
	dispersal_multipliers_matrix = dispersal_multipliers_matrix * manual_dispersal_multipliers_matrix


	
	# Set up the instantaneous rate matrix (Q matrix)
	# someday we'll have to put "a" (anagenic range-switching) in here...
	#Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, amat, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)

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
	# e.g. the j events, NOT the dmat above which is d*dispersal_multipliers_matrix
	spPmat_inputs$dmat = dispersal_multipliers_matrix
	dmat = dispersal_multipliers_matrix

	states_indices = states_list
	states_indices[1] = NULL	# shorten the states_indices by 1 (cutting the null range state from the speciation matrix)
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
		tmpvals = get_perEvent_probs(params_table, plotwhat=plotwhat)
		
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
	
	}



