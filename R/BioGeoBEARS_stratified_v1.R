require("ape")
require("rexpokit")
require("cladoRcpp")


#######################################################
# Do a stratified or other constrained analysis
#######################################################

# Based on:
# /_examples/changing_geog_v1.R




#######################################################
# section_the_tree
#######################################################
#' Section a tree for stratified analysis
#'
#' A utility function for stratified analysis.  Sections the tree into a series of strata. 
#' Each stratum may have one or more subtrees (APE phylo3 objects, *WITH* root edges) and/or
#' branch segments (which are just represented as numeric values, indicating the length of the sub-branch,
#' i.e. the time-width of the stratum, if the branch crosses the whole stratum.
#' 
#' @param inputs The list of inputs for stratified analysis
#' @param make_master_table If desired, make an \code{inputs$master_table} containing the
#' correspondance between the original tree and the sectioned pieces.
#' @param plot_pieces If \code{TRUE}, plot the tree chunks (but not isolated branch segments) as they are created.
#' @param cut_fossils If \code{TRUE} (default), the program is stopped if there are fossils, i.e. tips older than 0.6 my (default).  Users should
#' use code{\link[ape]{drop.tip}} or an external program to clip fossils out of the tree. PLEASE NOTE that several times I have experienced miserable long nights
#' due, apparently, to \code{\link[ape]{drop.tip}} producing weird tree structures, resulting in weird Newick files, without me realizing it.  The solution is usually to 
#' open the Newick file in something like \code{FigTree}, resort the branches, and save to a new Newick file.
#' Fossils have now been implemented in stratified analysis; this was complicated, as it involves inserting new branches in chopped trees.
#' @param fossils_older_than Tips that are older than \code{fossils_older_than} will be marked as \code{TRUE} in a column called \code{fossil}.
#' This is not currently set to 0, because Newick files can have slight precision issues etc. that mean not all tips quite come to zero.  You 
#' can attempt to fix this with \code{\link{extend_tips_to_ultrametricize}} (but make sure you do not inappropriately average in fossils!!).
#' @return \code{inputs} with \code{inputs$tree_sections_list} added.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}, \code{\link[ape]{drop.tip}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
section_the_tree <- function(inputs, make_master_table=FALSE, plot_pieces=TRUE, cut_fossils=TRUE, fossils_older_than=0.6)
	{
	runjunk='
	make_master_table=TRUE; plot_pieces=FALSE; cut_fossils=TRUE; fossils_older_than=0.6;
	'
	
	
	# Fixing nodes, for marginal local optimum ancestral state reconstruction, is COMPLICATED when you are
	# chopping up an APE tree.  Somehow we would have to keep track of which node.  So, save this for later.
	
	
	orig_timeperiods = inputs$timeperiods
	timeperiods = orig_timeperiods
	original_tree = read.tree(inputs$trfn)
	phy_as_it_is_chopped_down = original_tree
	
	# Make the tree table for the original tree
	orig_tr_table = prt(original_tree, printflag=FALSE, get_tipnames=TRUE)
	orig_tr_table

	# Identify fossils
	tipnums = 1:length(original_tree$tip.label)
	fossils_TF = orig_tr_table$time_bp[tipnums] >= fossils_older_than
	numfossils = sum(fossils_TF)
	fossil_names = original_tree$tip.label[fossils_TF]

	if (numfossils > 0)
		{
		if (cut_fossils == TRUE)
			{
			# Stop the analysis so that the user may cut the fossils.
			stoptxt = cat("\n\nFATAL ERROR in section_the_tree(): Your tree has ", numfossils, " fossil tips older than ", fossils_older_than, " my!\n",
							"But you have not turned on fossils by setting cut_fossils=FALSE in section_the_tree()\n", 
							"Fossil tipnames listed below:\n", sep="")
			cat(stoptxt)
			print(fossil_names)
			
			# Warn about drop.tip
			cat("\n\nAlso: PLEASE NOTE that several times I have experienced miserable long nights due, apparently, to drop.tip producing weird tree structures, resulting in weird Newick files, without me realizing it.  The solution is usually to open the Newick file in something like FigTree, resort the branches, and save to a new Newick file.\n\n")
			
			stop(stoptxt)
		
		
			junk='
			tr_nofossils = drop.tip(original_tree, fossil_names)
			write.tree(tr_nofossils, file="venerid_tree_for_biogeog_v1.newick")
			'
			} else {
			# The simplest approach to INCLUDING fossils is to artificially extend the branchlengths
			warntxt = cat("\n\nWARNING: Your tree has ", numfossils, " fossil tips older than ", fossils_older_than, " my!\n",
							"Make sure that 'fossils_older_than' is set to capture all of the fossil tips in your tree!\n", 
							"(default: fossils_older_than=0.6)\n", 
							"Fossil tipnames listed below:\n", sep="")
			cat(warntxt)
			cat(paste(fossil_names, collapse="\n", sep=""))
			cat("\n\n")
			
			# This will extend ALL tips up to time_bp=0 my.  Keep track of true tip age through orig_tr_table$fossils and orig_tr_table$time_bp
			phy_as_it_is_chopped_down = extend_tips_to_ultrametricize(obj=phy_as_it_is_chopped_down, age_of_root=0, tips_end_at_this_date=NA)
			}
		}
	
	# Make a master table of how the pieces correspond to the original tree!
	if (make_master_table == TRUE)
		{
		master_table = NULL
		}
	
	
	if (plot_pieces == TRUE)
		{
		plot(phy_as_it_is_chopped_down)
		#abline(v=timeperiods)
		axisPhylo()
		}
	# CHECK THIS FUNCTION
	#phy_as_it_is_chopped_down$edge.length = phy_as_it_is_chopped_down$edge.length + 0.0001
	
	
	tree_sections_list = NULL
	tnum = 0
	
	if (length(timeperiods) <= 1)
		{
		chainsaw_result = list()
		chainsaw_result$tree_to_chainsaw = phy_as_it_is_chopped_down
		chainsaw_result$return_pieces_list[[1]] = phy_as_it_is_chopped_down
		
		# Merge THEN split THEN sort!!
		# Make sure to sort the names before merging
		tmp_sorted_names_merge = paste(phy_as_it_is_chopped_down$tip.label, collapse=",", sep="")
		tmp_sorted_names_split = strsplit(x=tmp_sorted_names_merge, split=",")[[1]]
		chainsaw_result$return_pieces_basenames[[1]] = paste(sort(tmp_sorted_names_split), collapse=",", sep="")
		attr(chainsaw_result, "class") = "chainsaw_result"
		tree_sections_list[[1]] = chainsaw_result
		} else {
		
		# Instead of using these column names, which might change:
		# c(1,4:8,10)
		# ...use col headings (and add fossils)
		table_colnames = c("node", "node.type", "parent_br", "edge.length", "ancestor", "daughter_nds", "time_bp", "fossils", "label")
		SUBtable_colnames = paste("SUB", table_colnames, sep="")
		
		# Put the tips into the master condlikes table (important if we have ambiguous tips)
		if (make_master_table == TRUE)
			{
			# First, put in the original tree tips
			orig_tips_table = orig_tr_table[1:length(original_tree$tip.label), ]
			subtree_table = orig_tips_table
			names(subtree_table) = paste("SUB", names(subtree_table), sep="")
			
			# Get the relative timepoint
			stratum = 0
			reltimept = 0
			time_bot = 0
			time_top = 0
			piecenum = 0
			piececlass = "orig_tip"
			subtree_table = cbind(stratum, time_top, time_bot, reltimept, piecenum, piececlass, subtree_table[,SUBtable_colnames])
			subtree_table$SUBnode.type = "orig_tip"		
			
			tmp_join_table = cbind(orig_tips_table[, table_colnames], subtree_table)
			tmp_join_table	
			master_table = rbind(master_table, tmp_join_table)
			}

		for (i in 1:(length(timeperiods)))
		#for (i in 1:3))
			{
			# Label the stratum
			stratum = i
			
			
			cat("\n", i, ": ", timeperiods[i], "\n", "\n", sep="")
			# Chainsaw the top off the tree
			if (i == 1)
				{
				timepoint = timeperiods[i] - 0
				} else {
				timepoint = timeperiods[i]# - timeperiods[i-1]
				}
			# Update timepoints so you are subtracting the right amount!!!!!!!!
			timeperiods = timeperiods - timepoint
			timeperiods
		  
			# Check if you are in the last timeperiod
			if (i < length(timeperiods))
				{
				# Otherwise, CHAINSAW the sucker!
				chainsaw_result = chainsaw2(phy_as_it_is_chopped_down, timepoint=timepoint, return_pieces=TRUE)
				#print(chainsaw_result)
				} else {
				# If it's the last piece, just use the remaining leftover tree chunk
				chainsaw_result = list()
				chainsaw_result$tree_to_chainsaw = phy_as_it_is_chopped_down
				chainsaw_result$return_pieces_list[[1]] = phy_as_it_is_chopped_down
				
				# Merge THEN split THEN sort!!
				# This may not be necessary; but what the heck.
				# Make sure to sort the names before merging
				tmp_sorted_names_merge = paste(phy_as_it_is_chopped_down$tip.label, collapse=",", sep="")
				tmp_sorted_names_split = strsplit(x=tmp_sorted_names_merge, split=",")[[1]]
				
				chainsaw_result$return_pieces_basenames[[1]] = paste(sort(tmp_sorted_names_split), collapse=",", sep="")
				attr(chainsaw_result, "class") = "chainsaw_result"
				}
			
			# Store the chainsaw result
			tree_sections_list[[(tnum=tnum+1)]] = chainsaw_result

			# Make a master table of how the pieces correspond to the original tree!
			if (make_master_table == TRUE)
				{
				# Update the corresponding table
				tipnames_above_cutpoints = unlist(chainsaw_result$return_pieces_basenames)
				tipnames_above_cutpoints
	
				# Find the position of this subbranch (its top node) in the overall tree
				pos_of_1st_in_2nd = match(tipnames_above_cutpoints, orig_tr_table$tipnames)
				pos_of_1st_in_2nd
	
				classes_of_pieces = sapply(X=chainsaw_result$return_pieces_list, FUN=class)
				classes_of_pieces[classes_of_pieces == "numeric"] = "subbranch"
				classes_of_pieces[classes_of_pieces == "phylo"] = "subtree"
				
				
				# Get the tree structure as the tree is chopped down
				print(i)
				phy_chopped_down_table = prt(phy_as_it_is_chopped_down, printflag=FALSE, get_tipnames=TRUE)
				
				# re-sort the tipnames
				for (rownum in 1:nrow(phy_chopped_down_table))
					{
					temp_tipnames = phy_chopped_down_table$tipnames[rownum]
					words = strsplit(temp_tipnames, split=",")[[1]]
					words = sort(words)
					phy_chopped_down_table$tipnames[rownum] = paste(words, collapse=",", sep="")
					}
				
				
				# Get the relative timepoint
				reltimept = timepoint
				time_bot = orig_timeperiods[i]
				time_top = time_bot - reltimept
				
				# Accumulate the rows of the table
				# Go through the pieces
				for (p in 1:length(classes_of_pieces))
					{
					# For subtrees, get all the corresponding node info
					if (classes_of_pieces[p] == "subtree")
						{
						# Get the nodenums in the subtree that's been removed
						tmp_subtree = chainsaw_result$return_pieces_list[[p]]
						subtree_table = prt(tmp_subtree, printflag=FALSE, get_tipnames=TRUE)
						
						# re-sort the tipnames
						for (rownum in 1:nrow(subtree_table))
							{
							temp_tipnames = subtree_table$tipnames[rownum]
							words = strsplit(temp_tipnames, split=",")[[1]]
							words = sort(words)
							subtree_table$tipnames[rownum] = paste(words, collapse=",", sep="")
							}
						
						names(subtree_table) = paste("SUB", names(subtree_table), sep="")
						subtree_table
						
						# Identify the corresponding nodes						
						tree_piece_nodenums = subtree_table$SUBnode
						tiplabels_for_each_node_in_tree_piece = subtree_table$SUBtipnames
						
						pos_of_1st_in_2nd = match(tiplabels_for_each_node_in_tree_piece, orig_tr_table$tipnames)
						pos_of_1st_in_2nd
						
						
						# Add the pieces identifiers
						piecenum = p
						piececlass = classes_of_pieces[p]
						subtree_table = cbind(stratum, time_top, time_bot, reltimept, piecenum, piececlass, subtree_table[,SUBtable_colnames])
						subtree_table

						tmp_join_table = cbind(orig_tr_table[pos_of_1st_in_2nd, table_colnames], subtree_table)
						tmp_join_table						

						# NA check
						if (is.na(tmp_join_table[1,1]) == TRUE)
							{
							stoptxt = "\n\nFATAL ERROR #1 produced in section_the_tree(): NAs in tmp_join_table.\n\n"
							cat(stoptxt)
							
							print("i")
							print(i)
							print("p")
							print(p)
							
							print(tmp_join_table)
							stop(stoptxt)
							}

						} else {
						# For sub-branches, just add 1 row
						# Get the nodenums in the subtree that's been removed
						tmp_subbranch = chainsaw_result$return_pieces_list[[p]]
# 						subtree_table = prt(phy_as_it_is_chopped_down, printflag=FALSE, get_tipnames=TRUE)
# 						names(subtree_table) = paste("SUB", names(subtree_table), sep="")
# 						subtree_table
						
						# Identify the corresponding nodes						
						tree_piece_nodenums = 1
						tmp_basenames = chainsaw_result$return_pieces_basenames[[p]]
						
						# This may not be necessary; but what the heck.
						tmp_basenames2 = paste(tmp_basenames, collapse=",", sep="")
						tmp_basenames3 = strsplit(x=tmp_basenames2, split=",")[[1]]
						tiplabels_for_each_node_in_tree_piece = paste(sort(tmp_basenames3), collapse=",", sep="")
						
						pos_of_1st_in_2nd = match(tiplabels_for_each_node_in_tree_piece, phy_chopped_down_table$tipnames)
						pos_of_1st_in_2nd
						
						# Use the chopped-down-tree to reference isolated branches (may not matter)
						subtree_table = phy_chopped_down_table
						names(subtree_table) = paste("SUB", names(subtree_table), sep="")
						subtree_table
						
						# Add the pieces identifiers
						piecenum = p
						piececlass = classes_of_pieces[p]
						subtree_table = cbind(stratum, time_top, time_bot, reltimept, piecenum, piececlass, subtree_table[pos_of_1st_in_2nd, SUBtable_colnames])

						
						# Find the reference to the master tree
						pos_of_1st_in_2nd = match(tiplabels_for_each_node_in_tree_piece, orig_tr_table$tipnames)
						pos_of_1st_in_2nd
						
						tmp_join_table = cbind(orig_tr_table[pos_of_1st_in_2nd, table_colnames], subtree_table)
						tmp_join_table	
						
						# NA check
						if (is.na(tmp_join_table[1,1]) == TRUE)
							{
							stoptxt = "\n\nFATAL ERROR #2 produced in section_the_tree(): NAs in tmp_join_table.\n\n"
							cat(stoptxt)
							
							print("i")
							print(i)
							print("p")
							print(p)
							print(tmp_join_table)
							
							
							print(tiplabels_for_each_node_in_tree_piece)
							
							#print(orig_tr_table$tipnames)
							
							print(pos_of_1st_in_2nd)
							
							stop(stoptxt)
							}
						}
					
					master_table = rbind(master_table, tmp_join_table)
					}
				}
				


	
			# Convey the tree to the next round of chopping
			phy_as_it_is_chopped_down = chainsaw_result$tree_to_chainsaw
			
			if (plot_pieces == TRUE)
				{
				plot(phy_as_it_is_chopped_down)
				#axisPhylo2(side = 1, roundlabels=TRUE, minage=timeperiods[i] 
				axisPhylo()
				}
			
			}
		}

		
	# Append to inputs and return
	inputs$tree_sections_list = tree_sections_list

	# Also append the master table
	inputs$master_table = master_table

	return(inputs)
	}







#######################################################
# chainsaw2
#######################################################
#' Saw a tree off at a particular time before present
#' 
#' This function chops a tree like a hedge-trimmer, cutting straight across at a particular timepoint. 
#' The pieces are returned, as is the leftover tree, with branches shortened appropriately.  Pieces
#' that are mini-trees are returned as ape objects, whereas single branches are just lengths.
#'
#' This function is used during stratification, but could have other uses as well.
#' 
#' @param tr An ape phylo object.
#' @param timepoint The time at which the tree should be "chopped".
#' @param return_pieces Default TRUE, which means pieces should be returned
#' @return \code{chainsaw_result} (a list object with the pieces) or \code{tree_to_chainsaw}, just the leftover tree
#' @export
#' @seealso \code{\link{section_the_tree}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
chainsaw2 <- function(tr, timepoint=10, return_pieces=TRUE)
	{
	# Take a tree and saw it off evenly across a certain timepoint.
	# This removes any tips above the timepoint, and replaces them 
	# with a single tip representing the lineage crossing
	# the timepoint (with a new tip name).

	# Get the tree in a table
	tr_table = prt(tr, printflag=FALSE, get_tipnames=FALSE)
	tr_table

	# Find the tips that are less than 10 my old and drop them
	TF_exists_more_recently_than_10mya = tr_table$time_bp < timepoint
	
	# Get the corresponding labels
	labels_for_tips_existing_more_recently_than_10mya = tr_table$label[ TF_exists_more_recently_than_10mya == TRUE ]
	
	###########################################
	# Draft chainsaw function
	###########################################
	# loop through the branches that cross 10 mya
	
	# get a list of the edge start/stops in the phylogeny's edges
	edge_times_bp = get_edge_times_before_present(tr)
	
	# which of these branches cross 10 mya (or whatever timepoint)?
	edges_start_earlier_than_10mya = edge_times_bp[, 1] > timepoint
	edges_end_later_than_10mya = edge_times_bp[, 2] <= timepoint
	edges_to_chainsaw = edges_start_earlier_than_10mya + edges_end_later_than_10mya == 2
	
	# then, for each of these edges, figure out how many tips exist descending from it
	# these are the nodes ABOVE the cutoff line
	nodes_to_chainsaw = tr$edge[, 2][edges_to_chainsaw]
	

	# Take only internal nodes (? why ?)
	numtips = length(tr$tip.label)
	#nodes_to_chainsaw = nodes_to_chainsaw[nodes_to_chainsaw > numtips]
	
	# create a copy of the tree to chainsaw
	tree_to_chainsaw = tr
	
	if (return_pieces == TRUE)
		{
		return_pieces_list = as.list(rep(NA, length(nodes_to_chainsaw)))
		return_pieces_basenames = as.list(rep(NA, length(nodes_to_chainsaw)))
		chopTable = NULL
		
		}
	
	chainsaw_table = NULL
	
	for (i in 1:length(nodes_to_chainsaw))
		{
		# If this is a tip node on the current tree, shorten the branch rather than cut it off
		if (nodes_to_chainsaw[i] <= numtips)
			{
			# Here, chainsaw is cutting an internal node, so extract the sectioned branch before you cut it down
			# (the cutting happens after the forloop)
			# (This is easy, it is just the length of the timeslab;
			#  which you should UPDATE as you move down the tree
			if (return_pieces == TRUE)
				{
				# Record the length of the branch section, and the name of that tip
				# (which is also the name of that base)
				return_pieces_list[[i]] = timepoint
				tmp_tipname = tr$tip.label[nodes_to_chainsaw[i]]
				return_pieces_basenames[[i]] = tmp_tipname
				}
			# You don't have to do anything else, the chopping of single branches is 
			# covered after the forloop
			#cat("\ni=", i, "	ntips=", length(tree_to_chainsaw$tip.label), sep="")

			} else {
			# Here, it's an internal node, so extract the subtree before you drop it
			tmp_subtree = extract.clade(tr, nodes_to_chainsaw[i])
			#plot(tmp_subtree, root.edge=TRUE)
			# Also, record the branchlength below this node
			branchlength_below_subtree_LCA_node = timepoint - get_max_height_tree(tmp_subtree)
			# Add this to the bottom of the subtree
			tmp_subtree$root.edge = branchlength_below_subtree_LCA_node
			#plot(tmp_subtree, root.edge=TRUE)
			
			# Record the piece, if desired
			if (return_pieces == TRUE)
				{
				# Record the length of the branch section, and the name of that tip
				# (which is also the name of that base)
				return_pieces_list[[i]] = tmp_subtree

				# Merge THEN split THEN sort!!
				tmp_labels_merge = paste(tmp_subtree$tip.label, collapse=",", sep="")
				tmp_labels_split = strsplit(tmp_labels_merge, split=",")[[1]]
				new_labels = sort(tmp_labels_split)
				basename_after_cutting = paste(new_labels, collapse=",", sep="")
				return_pieces_basenames[[i]] = basename_after_cutting				
				}

			#print(tmp_subtree$tip.label)

			tmp_number_of_tips = length(tmp_subtree$tip.label)
			#print(tmp_number_of_tips)
			
			# number of tips to drop = (numtips -1)
			numtips_to_drop = tmp_number_of_tips - 1 
			
			# tips_to_drop
			tmp_labels = tmp_subtree$tip.label
			
			labels_to_drop = tmp_labels[1:numtips_to_drop]
			ordered_labels_to_make_into_new_name = sort(tmp_labels)
			name_new_tip = paste(ordered_labels_to_make_into_new_name, collapse=",", sep="")
			
			# new label
			label_kept_num = length(tmp_labels)
			label_kept = tmp_labels[label_kept_num]
			#new_label = paste("CA_", label_kept, "+", numtips_to_drop, "_tips", sep="")
			new_label = name_new_tip
			tree_to_chainsaw$tip.label[tree_to_chainsaw$tip.label == label_kept] = new_label

			# chop off e.g. 2 of the 3 tips
			tree_to_chainsaw = drop.tip(tree_to_chainsaw, labels_to_drop)
			#cat("\ni=", i, "	ntips=", length(tree_to_chainsaw$tip.label), sep="")
			} # end else
		} # end for loop
	#plot(tree_to_chainsaw)
	#axisPhylo()
	
	tree_to_chainsaw_table = prt(tree_to_chainsaw, printflag=FALSE)
	
	tree_to_chainsaw_table_tips_TF_time_bp_LT_10my = tree_to_chainsaw_table$time_bp < timepoint
	
	
	tmp_edge_lengths =  tree_to_chainsaw_table$edge.length[tree_to_chainsaw_table_tips_TF_time_bp_LT_10my]
	
	times_bp_for_edges_to_chainsaw = tree_to_chainsaw_table$time_bp[tree_to_chainsaw_table_tips_TF_time_bp_LT_10my]
	
	adjustment = times_bp_for_edges_to_chainsaw - timepoint
	
	revised_tmp_edge_lengths = tmp_edge_lengths + adjustment
	
	tree_to_chainsaw_table$edge.length[tree_to_chainsaw_table_tips_TF_time_bp_LT_10my] = revised_tmp_edge_lengths
	
	# revised
	ordered_nodenames = get_nodenums(tree_to_chainsaw)
	parent_branches = get_indices_where_list1_occurs_in_list2(ordered_nodenames, tree_to_chainsaw$edge[,2])
	
	NA_false = is.not.na(tree_to_chainsaw_table$edge.length)
	
	tree_to_chainsaw$edge.length[parent_branches[NA_false]] = tree_to_chainsaw_table$edge.length[NA_false]

	if (return_pieces == TRUE)
		{
		chainsaw_result = NULL
		chainsaw_result$tree_to_chainsaw = tree_to_chainsaw
		chainsaw_result$return_pieces_list = return_pieces_list
		chainsaw_result$return_pieces_basenames = return_pieces_basenames
		class(chainsaw_result) = "chainsaw_result"
		return(chainsaw_result)
		} else {
		return(tree_to_chainsaw)
		}
	}




#######################################################
# calc_loglike_sp_stratified:
#######################################################
#' Calculate log-likelihood with a transition matrix and speciation events, and with stratification
#'
#' This function is the stratified version of \code{\link{calc_loglike_sp}}.
#' 
#' @param tip_condlikes_of_data_on_each_state A numeric matrix with rows representing tips, and columns representing states/geographic ranges.  The cells
#' give the likelihood of the observation data under the assumption that the tip has that state; typically this means that the known geographic range gets a 
#' '1' and all other states get a 0.
#' @param phy A phylogeny object.  The function converts it to pruningwise order.
#' @param Qmat A Q transition matrix representing the along-branch model for the evolution of geographic range, using parameters \emph{d} (dispersal/range expansion), 
#' \emph{e} (extinction/range contraction/local extirpation), and perhaps others (e.g. distance).  This matrix can be input in either dense or sparse (COO) format, 
#' as specified by \code{input_is_COO}.
#' @param spPmat Default is \code{NULL}; users should usually use \code{spPmat_inputs}.  \code{spPmat} is A numeric matrix representing the probability of each
#' ancestor range-->(Left range, Right range) transition at cladogenesis events.  There are 
#' different ways to represent this matrix.  In the simplest representation, this is just a rectangular matrix with numstates rows (representing the ancestral
#' states) and numstates^2 columns (representing all possible descendant pairs).  Use of this type of matrix is specified by \code{cppSpMethod=1}. It is calculated
#' from a textual speciation matrix (typically \code{spmat} in the code) via \code{\link{symbolic_to_relprob_matrix_sp}}. However, this matrix gets huge and
#' slow for large numbers of states/ranges.  \code{cppSpMethod=2} and \code{cppSpMethod=3} implement successively more efficient and faster 
#' representation and processing of this matrix in COO-like formats.  See \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOprobs}} for the \code{cppSpMethod=2} 
#' method, and \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOweights_faster}} for the \code{cppSpMethod=3} method (the fastest).
#' @param min_branchlength Nodes with branches below this branchlength will not be treated as cladogenesis events; instead, they will be treated as 
#' if an OTU had been sampled from an anagenetic lineage, i.e. as if you had a direct ancestor.  This is useful for putting fossils into the biogeography analysis,
#' when you have fossil species that range through time. (Note: the proper way to obtain such trees, given that most phylogenetic methods force all OTUs to be tips 
#' rather than direct ancestors, is another question subject to active research.  However, one method might be to just set a branch-length cutoff, and treat any
#' branches sufficiently small as direct ancestors.)
#' @param return_what What should be returned to the user? Options are "loglike" (the log-likelihood of the data under the tree, model, and model parameters), 
#' "nodelikes" (the scaled conditional likelihoods at the nodes), "rootprobs" (the relative probability of the geographic ranges/states at the root), or "all"
#' (all of the above in a list).  Typically the user will only want to return "loglike" while doing ML optimization, but then return "all" once the ML parameter
#' values have been found.
#' @param probs_of_states_at_root The prior probability of the states/geographic ranges at the root.  The default, \code{NULL}, effectively means an equal probability
#' for each state (this is also what \code{LAGRANGE} assumes; and running with NULL will reproduce exactly the \code{LAGRANGE} parameter inferences and
#' log-likelihood).
#' @param rootedge  Should the root edge be included in the calculation (i.e., calculate to the bottom of the root), if a root edge is present?  Default \code{FALSE}.
#' @param sparse Should sparse matrix exponentiation be performed?  This should be faster for very large matrices (> 100-200 states), however, the calculations 
#' appear to be less accurate.  The function will transform a dense matrix to COO format (see \code{\link[rexpokit]{mat2coo}}) if necessary according to 
#' the \code{input_is_COO} parameter.
#' @param printlevel If >= 1, various amounts of intermediate output will be printed to screen.  Note: Intermediate outputs from C++ and FORTRAN functions have been
#' commented out, to meet CRAN guidelines.
#' @param use_cpp Should the C++ routines from \code{\link[cladoRcpp]{cladoRcpp}} be used to speed up calculations?  Default \code{TRUE}.
#' @param input_is_COO Is the input Q matrix a sparse, COO-formatted matrix (\code{TRUE}) or a standard dense matrix (\code{FALSE}). Default \code{FALSE}.
#' @param spPmat_inputs A list of parameters so that \code{spPmat} (the speciation transition probability matrix) can be calculated on-the-fly, according
#' to the method in \code{cppSpMethod}.  See example.
#' @param cppSpMethod Three C++ methods from cladoRcpp for calculating and using the cladogenesis probability matrix.  1 is slowest but easiest to understand; 3 is fastest.
#' If \code{spPmat_inputs} is given, the program will generate the appropriate spPmat on-the-fly, and the user does not have to input the full \code{spPmat} manually.
#' @param cluster_already_open If the user wants to distribute the matrix exponentiation calculations from all the branches across a number of processors/nodes on 
#' a cluster, specify the cluster here.  E.g. \code{cluster_already_open = makeCluster(rep("localhost",num_cores_to_use), type = "SOCK")}.  Note: this will work on 
#' most platforms, including Macs running R from command line, but will NOT work on Macs running the R GUI \code{R.app}, because parallel processing functions like
#' \code{MakeCluster} from e.g. \code{library(parallel)} for some reason crash R.app.  The program runs a check for R.app and will just run on 1 node if found. 
#' @param calc_ancprobs Should ancestral state estimation be performed (adds an uppass at the end).
#' @param null_range_allowed Does the state space include the null range?
#' Default is \code{NULL} which means running on a single processor.
#' @param fixnode If the state at a particular node is going to be fixed (e.g. for ML marginal ancestral states), give the node number. 
#' (Trial implementation for stratified analysis.)
#' @param fixlikes The state likelihoods to be used at the fixed node.  I.e. 1 for the fixed state, and 0 for the others.
#' (Trial implementation for stratified analysis.)
#' @param inputs A list of inputs containing the dispersal matrix for each time period, etc.
#' @param allareas A list of all the areas in the total analysis
#' @param all_states_list A list of all the stats in the total analysis (0-based coding - ?)
#' @param return_condlikes_table If \code{TRUE}, return the table of ALL conditional likelihood results, including at branch subsections
#' (only some should be used in calculating the final log-likelihood of the geography range data on the tree!)
#' @param calc_TTL_loglike_from_condlikes_table If TRUE, force making of the condlikes table, and use it to calculate the log-likelihood
#' (default=TRUE; matches LAGRANGE).
#' @return grand_total_likelihood The total log-likelihood of the data on the tree (default). Or, if 
#' \code{return_condlikes_table==TRUE}, the function returns \code{calc_loglike_sp_stratified_results}, with 
#' \code{calc_loglike_sp_stratified_results$condlikes_table} and \code{calc_loglike_sp_stratified_results$grand_total_likelihood}
#' as list items.  This can be useful for debugging stratified analyses, which have a lot of extra book-keeping that is easy to mess up.
#' @export
#' @seealso \code{\link{calc_loglike_sp}}, \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp}}, \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOprobs}}, 
#' \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOweights_faster}}, \code{\link[rexpokit]{mat2coo}}, 
#' \code{\link{rcpp_calc_anclikes_sp_COOweights_faster}}
#' @bibliography /Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#'	 @cite Landis_Matzke_etal_2013_BayArea
#' @note Go BEARS!
#' @note (COO = Coordinate list format for a matrix, see \url{http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29}
#' @author Nicholas Matzke \email{matzke@@berkeley.edu}
#' @examples
#' testval=1
#'
calc_loglike_sp_stratified <- function(tip_condlikes_of_data_on_each_state, phy, Qmat=NULL, spPmat=NULL, min_branchlength=1e-21, return_what="loglike", probs_of_states_at_root=NULL, rootedge=TRUE, sparse=FALSE, printlevel=0, use_cpp=TRUE, input_is_COO=FALSE, spPmat_inputs=NULL, cppSpMethod=3, cluster_already_open=NULL, calc_ancprobs=FALSE, null_range_allowed=TRUE, fixnode=NULL, fixlikes=NULL, inputs=inputs, allareas=allareas, all_states_list=all_states_list, return_condlikes_table=FALSE, calc_TTL_loglike_from_condlikes_table=TRUE)
	{
	defaults='
	Qmat=NULL; spPmat=NULL; min_branchlength=1e-21; return_what="loglike"; probs_of_states_at_root=NULL; rootedge=FALSE; sparse=FALSE; printlevel=1; use_cpp=TRUE; input_is_COO=FALSE; spPmat_inputs=NULL; cppSpMethod=3; cluster_already_open=NULL; calc_ancprobs=FALSE; null_range_allowed=TRUE; fixnode=NULL; fixlikes=NULL; inputs=inputs; allareas=allareas; all_states_list=all_states_list; return_condlikes_table=FALSE
	'
	
	defaults='
	maxareas = 4
	phy = read.tree(inputs$trfn)
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=np(inputs$geogfn))
	tip_condlikes_of_data_on_each_state = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges, phy, maxareas=maxareas)
	
	allareas = getareas_from_tipranges_object(tipranges)
	all_states_list = rcpp_areas_list_to_states_list(areas=allareas, include_null_range=TRUE, maxareas=maxareas)
	
	tmpres = calc_loglike_sp_stratified(tip_condlikes_of_data_on_each_state, phy, Qmat=NULL, spPmat=NULL, min_branchlength=1e-21, return_what="all", probs_of_states_at_root=NULL, rootedge=TRUE, sparse=FALSE, printlevel=0, use_cpp=TRUE, input_is_COO=FALSE, spPmat_inputs=NULL, cppSpMethod=3, cluster_already_open=NULL, calc_ancprobs=FALSE, null_range_allowed=TRUE, fixnode=NULL, fixlikes=NULL, inputs=inputs, allareas=allareas, all_states_list=all_states_list, return_condlikes_table=FALSE)
	tmpres
	
	min_branchlength=1e-21
	null_range_allowed=TRUE
	printlevel=0
	cppSpMethod=3
	return_condlikes_table=TRUE
	calc_TTL_loglike_from_condlikes_table=TRUE
	calc_ancprobs=TRUE
'
	defaults='
	tmpres = calc_loglike_sp_stratified(tip_condlikes_of_data_on_each_state, phy, Qmat=NULL, spPmat=NULL, min_branchlength=1e-21, return_what="all", probs_of_states_at_root=NULL, rootedge=TRUE, sparse=FALSE, printlevel=0, use_cpp=TRUE, input_is_COO=FALSE, spPmat_inputs=NULL, cppSpMethod=3, cluster_already_open=NULL, calc_ancprobs=TRUE, null_range_allowed=TRUE, fixnode=NULL, fixlikes=NULL, inputs=inputs, allareas=allareas, all_states_list=all_states_list, return_condlikes_table=TRUE, calc_TTL_loglike_from_condlikes_table=TRUE)
'



# 	defaults='
# 	# STANDARD DEBUGGING HERE
# 	tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state; phy=phy; Qmat=NULL; spPmat=NULL; min_branchlength=1e-21; return_what="loglike"; probs_of_states_at_root=NULL; rootedge=TRUE; sparse=FALSE; printlevel=0; use_cpp=TRUE; input_is_COO=FALSE; spPmat_inputs=NULL; cppSpMethod=3; cluster_already_open=NULL; calc_ancprobs=FALSE; null_range_allowed=TRUE; fixnode=fixnode; fixlikes=fixlikes; inputs=BioGeoBEARS_run_object; allareas=areas; all_states_list=states_list; return_condlikes_table=TRUE; calc_TTL_loglike_from_condlikes_table=TRUE;
# 	' # end junk


	if ((return_condlikes_table == TRUE) || (calc_TTL_loglike_from_condlikes_table == TRUE))
		{
		names_in_inputs = names(inputs)			# can't use exists() on list items; reasons explained here:
												# http://stackoverflow.com/questions/7719741/how-to-test-if-list-element-exists
		if ( ("master_table" %in% names_in_inputs) == TRUE)
			{
			condlikes_table = matrix(data=0, nrow=nrow(inputs$master_table), ncol=length(all_states_list))
			
			# Put in the conditional likelihoods at the tips
			tmprownums = nrow(tip_condlikes_of_data_on_each_state)
			condlikes_table[1:tmprownums, ] = tip_condlikes_of_data_on_each_state
			
			if (calc_ancprobs == TRUE)
				{
				relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_table = matrix(data=0, nrow=nrow(inputs$master_table), ncol=length(all_states_list))
				}
			
			} else {
			cat("\n\nWARNING: in 'calc_loglike_sp_stratified()', you set 'return_condlikes_table=TRUE'\n
			and/or calc_TTL_loglike_from_condlikes_table=TRUE, but this requires that\n
			'inputs$master_table' be available from the 'section_the_tree()' function. Try
			\nrunning 'inputs=section_the_tree(inputs, make_master_table=TRUE).\n", sep="")
			
			cat("\nAs a result, we are setting return_condlikes_table=FALSE\n\n", sep="")
			return_condlikes_table=FALSE
			}
		}

	# Get the timeperiods; if 1 time period, run once; if multiple, run 
	if (is.null(inputs$timeperiods) || length(inputs$timeperiods) == 1)
		{
		num_iterations = 1
		} else {
		# Multiple timeperiods
		timeperiods = inputs$timeperiods
		num_iterations = length(timeperiods)
		}
	
	# All areas in the total analysis
	allareas=allareas
	allareas_list = seq(0, length(allareas)-1, 1)		# 0-base indexes

	# This is used on the uppass -- it might change, if we re-write this to have changing # of areas within the stratum
	areas = allareas_list

	# All states in the total analysis (after e.g. limitation on total # of areas)
	all_states_list=all_states_list	
	
	# Other variables
	# sparse should probably be false for ancestral states/downpass/uppass considerations
	BioGeoBEARS_model_object = inputs$BioGeoBEARS_model_object
	force_sparse = sparse
	
	#######################################################
	# Set up the starting probabilities etc.
	#######################################################
	# Starting tip_relative_probs_of_each_state
	#current_condlikes_row = 0
	tip_relative_probs_of_each_state = tip_condlikes_of_data_on_each_state / rowSums(tip_condlikes_of_data_on_each_state)
	tip_relative_probs_of_each_state
	current_tip_relative_probs_of_each_state = tip_relative_probs_of_each_state
	current_condlikes_row = nrow(current_tip_relative_probs_of_each_state)
	
	
	# matrix to hold all of the relative probabilities; Making this purposely too big
	numnodes = phy$Nnode + length(phy$tip.label)
	all_relative_probs_of_each_state = matrix(0, ncol=length(all_states_list), nrow=(numnodes*length(timeperiods)))
	all_condlikes_of_each_state = matrix(0, ncol=length(all_states_list), nrow=(numnodes*length(timeperiods)))
	
	all_relative_probs_of_each_state[1:current_condlikes_row, ] = current_tip_relative_probs_of_each_state
	all_condlikes_of_each_state[1:current_condlikes_row, ] = current_tip_relative_probs_of_each_state
	
	
	#######################################################
	# Take the original tree and scale the branchlengths by b (branch-length exponent)
	# b=0, all branches=1; b=1, all branches normal
	#######################################################
	previous_timepoint = 0
	original_phy = phy
	
	# This is pointless in a stratified analysis
	# b_branch_length_exponent = inputs$BioGeoBEARS_model_object@params_table["b", "est"]
	# Branch-length exponent (must be applied *after* tree has been sectioned!)
	#original_phy$edge.length = original_phy$edge.length ^ b_branch_length_exponent
	phy_as_it_is_chopped_down = original_phy

	#tiplikes_to_delete = list()

	
	for (i in 1:num_iterations)
		{
		#i=1
		#cat("\ni=",i, sep="")

		# Set the dispersal and extinction rate
		d = BioGeoBEARS_model_object@params_table["d","est"]
		e = BioGeoBEARS_model_object@params_table["e","est"]
		a = BioGeoBEARS_model_object@params_table["a","est"]
		
		
		#######################################################
		# CONVERT NULL RANGE FROM "_" TO NA -- CRUCIAL, CAUSES CRASH OTHERWISE!!
		#######################################################
# 		if (null_range_allowed == TRUE)
# 			{
# 			TF = all_states_list == "_"
# 			all_states_list[TF] = NA
# 			} else {
# 			TF = all_states_list == "_"
# 			all_states_list[TF] = NA			
# 			}

		#######################################################
		# Cut down the number of areas, by what is allowed
		# (it would be more efficient to do this once during setup, but probably no biggie)
		#######################################################
		# states_to_use_TF: states to use in Qmat, speciation models, etc.
		# states_allowed_TF: use this to zero out impossible ancestral states according to
		#                    areas_allowed matrix
		# 
		
		if ( (is.null(inputs$list_of_areas_allowed_mats) == FALSE))
			{
			areas_allowed_mat = inputs$list_of_areas_allowed_mats[[i]]
	
			states_allowed_TF = sapply(X=all_states_list, FUN=check_if_state_is_allowed, areas_allowed_mat)
			#states_to_use_TF = all_states_list %in% tmp_states_list
			
			if (null_range_allowed == TRUE)
				{
				states_allowed_TF[1] = TRUE
				}
			# NO; use all areas for this
			# states_to_use_TF = states_allowed_TF
			
			} else {
			# Make no change
			pass = 1
			#states_list = states_list
			states_allowed_TF = rep(TRUE, length(all_states_list))
			}
		# Use this for regular calculations (Qmat, speciation models, etc.)
		states_to_use_TF = rep(TRUE, length(all_states_list))
		
			
		#####################################################
		# Make the dedf matrix for this time period
		#####################################################
		# If there is a distance matrix, use the first one 
		# (non-stratified analysis, here)

		# If there is a distance matrix, use the first one 
		# (non-stratified analysis, here)
		if ( (is.null(inputs$list_of_distances_mats) == FALSE))
			{
			distances_mat = inputs$list_of_distances_mats[[i]]
			} else {
			# Default is all areas effectively equidistant
			distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))
			}
	
		# Get the exponent on distance, apply to distances matrix
		x = BioGeoBEARS_model_object@params_table["x","est"]
		dispersal_multipliers_matrix = distances_mat ^ x
	
		# Apply manual dispersal multipliers, if any
		# If there is a manual dispersal multipliers matrix, use the first one 
		# (non-stratified analysis, here)
		if ( (is.null(inputs$list_of_dispersal_multipliers_mats) == FALSE))
			{
			manual_dispersal_multipliers_matrix = as.matrix(inputs$list_of_dispersal_multipliers_mats[[i]])
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
		if ( (is.null(inputs$list_of_area_of_areas) == FALSE))
			{
			area_of_areas = inputs$list_of_area_of_areas[[i]]
			} else {
			# Default is all areas effectively equidistant
			area_of_areas = rep(1, length(areas))
			}
			
		# Get the exponent on extinction, apply to extinction modifiers	
		u = BioGeoBEARS_model_object@params_table["u","est"]
		extinction_modifier_list = area_of_areas ^ (1 * u)
		
		# Apply to extinction rate
		elist = extinction_modifier_list * rep(e, length(areas))
		
		
		
		
		
		# Calculate the Q matrix
		# someday we'll have to put "a" (anagenic range-switching) in here...
# 		if (is.null(Qmat))
# 			{
			Qmat_tmp = rcpp_states_list_to_DEmat(areas_list=allareas_list, states_list=all_states_list[states_to_use_TF], 
			dmat=dmat, elist=elist, amat=amat, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)
# 			} else {
# 			# If Qmat is pre-specified
# 			Qmat_tmp = Qmat
# 			}
		
		
		# Now. IF you have a subtree structure, you need to run this with a cladogenesis matrix, 
		# through calc_loglike_sp(), like normal.
		
		# If there's just one tree, store it in the object
		if (is.null(inputs$timeperiods) || length(inputs$timeperiods) == 1)
			{
			tr = read.tree(inputs$trfn)
			
			tree_to_chainsaw = NULL
			tree_to_chainsaw[[1]] = tr

			return_pieces_list = NULL
			return_pieces_list[[1]] = tr
			
			return_pieces_basenames = NULL
			
			# Merge THEN split THEN sort!!
			tmp_labels_merge = paste(tr$tip.label, collapse=",", sep="")
			tmp_labels_split = strsplit(tmp_labels_merge, split=",")[[1]]
			return_pieces_basenames[[1]] = paste(sort(tmp_labels_split), collapse=",", sep="")
			
			chainsaw_object = list()
			chainsaw_object$tree_to_chainsaw = tree_to_chainsaw
			chainsaw_object$return_pieces_list = return_pieces_list
			chainsaw_object$return_pieces_basenames = return_pieces_basenames
			attr(chainsaw_object, "class") = "chainsaw_result"
			
			inputs$tree_sections_list[[1]] = chainsaw_object
			}
		
		
		# OK, if you have a tree here, do that
		# if not, exp the branch
		
		#######################################################
		# Cladogenic model 
		#######################################################
		j = BioGeoBEARS_model_object@params_table["j","est"]
		ysv = BioGeoBEARS_model_object@params_table["ysv","est"]
		ys = BioGeoBEARS_model_object@params_table["ys","est"]
		v = BioGeoBEARS_model_object@params_table["v","est"]
		y = BioGeoBEARS_model_object@params_table["y","est"]
		s = BioGeoBEARS_model_object@params_table["s","est"]
		sum_SPweights = y + s + j + v
	
		maxent_constraint_01 = BioGeoBEARS_model_object@params_table["mx01","est"]
		
		# Text version of speciation matrix	
		maxent_constraint_01v = BioGeoBEARS_model_object@params_table["mx01v","est"]
		#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
			
		# Set the parameter controlling the size distribution of 
		# the smaller descendant species
		maxent01s_param = BioGeoBEARS_model_object@params_table["mx01s","est"]
		maxent01v_param = BioGeoBEARS_model_object@params_table["mx01v","est"]
		maxent01j_param = BioGeoBEARS_model_object@params_table["mx01j","est"]
		maxent01y_param = BioGeoBEARS_model_object@params_table["mx01y","est"]
	
	
		# Cladogenesis model inputs
		spPmat_inputs = NULL

		# Note that this gets the dispersal multipliers matrix, which is applied to 
		# e.g. the j events, NOT the dmat above which is d*dispersal_multipliers_matrix
		spPmat_inputs$dmat = dispersal_multipliers_matrix

		states_indices = all_states_list[states_to_use_TF]
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
	
			
		#######################################################
		# Go through the tree pieces
		#######################################################
		chainsaw_result = inputs$tree_sections_list[[i]]
		
		# You will need the new tip likelihoods of the new tree:
		current_tip_relative_probs_of_each_state
		new_tip_likelihoods = matrix(0, nrow=length(chainsaw_result$return_pieces_list), ncol=length(all_states_list))

		for (jj in 1:length(chainsaw_result$return_pieces_list))
			{
			#cat("\njj=",jj, sep="")
			treepiece = chainsaw_result$return_pieces_list[[jj]]

		
			############################################
			# It's just a branch section
			############################################
			if (is.numeric(treepiece))
				{
				do_exponentiation = TRUE	# default
	
				# Check for fossil
				# If you are storing ALL of the conditional likelihoods that were calculated
				if ((return_condlikes_table == TRUE) || (calc_TTL_loglike_from_condlikes_table == TRUE))
					{
					# Find the row in the master_table
					TF1 = inputs$master_table$stratum == i
					TF2 = inputs$master_table$piecenum == jj
					TF3 = inputs$master_table$piececlass == "subbranch"
					TF = (TF1 + TF2 + TF3) == 3
					
					# Find the row
					rownum = (1:nrow(condlikes_table))[TF]
					tmp_master_table_row = inputs$master_table[rownum, ]
					
					# Error check
					if (nrow(tmp_master_table_row) != 1)
						{
						stoptxt = paste("\n\nFATAL ERROR in stratified loglike calculation at i=", i, "; jj=", jj, "; ", 'inputs$master_table$piececlass == "subbranch"', 
						          "\nnrow(tmp_master_table_row) should =1 but instead =", nrow(tmp_master_table_row), "\n", sep="")
						stop(stoptxt)
						}
					
					# Now check if it's a fossil that appears in this time bin
					master_tip_time_bp = tmp_master_table_row$time_bp
					time_top = tmp_master_table_row$time_top
					time_bot = tmp_master_table_row$time_bot
					is_fossil = tmp_master_table_row$fossils
					
					# If this is TRUE, there's a match and the fossil tip appears in this time period
					if ( (master_tip_time_bp >= time_top) && (master_tip_time_bp < time_bot) && (is_fossil == TRUE))
						{
						# Shorten the branchlength by master_tip_time_bp-time_top
						amount_to_shorten_by = master_tip_time_bp-time_top
						treepiece = treepiece - amount_to_shorten_by
						do_exponentiation = TRUE
						}

					# If this is TRUE, this fossil hasn't occurred yet, and you are looking at the "phantom limb".
					# In this case, DON'T do matrix exponentiation, just copy the likelihoods down!!
					if ( master_tip_time_bp < time_top )
						{
						do_exponentiation = FALSE
						}
						
					# If FALSE, you're below all this and hopefully don't care
					}



				tipname = chainsaw_result$return_pieces_basenames[[jj]]
				tip_TF = phy_as_it_is_chopped_down$tip.label == tipname
				relative_probs_of_each_state_at_the_tip_of_this_branch = current_tip_relative_probs_of_each_state[tip_TF, states_to_use_TF]
	
	
				if (do_exponentiation == TRUE)
					{
					# t = treepiece
					independent_likelihoods_at_branch_section_bottom = expokit_dgpadm_Qmat2(times=treepiece, Qmat=Qmat_tmp,  transpose_needed=TRUE)
					#independent_likelihoods_at_branch_section_bottom = expokit_dgpadm_Qmat(Qmat=Qmat_tmp,  t=treepiece, transpose_needed=FALSE)
				
				
					conditional_likelihoods_at_branch_section_bottom = matrix(independent_likelihoods_at_branch_section_bottom %*% relative_probs_of_each_state_at_the_tip_of_this_branch, nrow=1)
				
					# Zero out impossible states according to areas_allowed
					conditional_likelihoods_at_branch_section_bottom[states_allowed_TF==FALSE] = 0
					} else {
					# Copying the tip likelihoods down
					conditional_likelihoods_at_branch_section_bottom = matrix(relative_probs_of_each_state_at_the_tip_of_this_branch, nrow=1)
					}
				
				
				# Test forward exponentiation instead...NO
				# independent_likelihoods_at_branch_section_bottom = expokit_dgpadm_Qmat2(times=treepiece, Qmat=Qmat_tmp,  transpose_needed=FALSE)
# 				#independent_likelihoods_at_branch_section_bottom = expokit_dgpadm_Qmat(Qmat=Qmat_tmp,  t=treepiece, transpose_needed=FALSE)
# 				
# 				
# 				conditional_likelihoods_at_branch_section_bottom = matrix(independent_likelihoods_at_branch_section_bottom %*% relative_probs_of_each_state_at_the_tip_of_this_branch, nrow=1)
# 				conditional_likelihoods_at_branch_section_bottom[1] = 0
# 				
				
				
				
				# Also, store the conditional likelihoods for all nodes in this subtree
				chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]] = conditional_likelihoods_at_branch_section_bottom
	
				# Relative probabilities -- jjust the new tip
				chainsaw_result$relative_probs_of_each_state_at_bottom_of_root_branch[[jj]] = conditional_likelihoods_at_branch_section_bottom / sum(conditional_likelihoods_at_branch_section_bottom)
	
				# Relative probabilities -- all nodes plus branch bottom (jjust branch bottom, here)
				chainsaw_result$relative_probabilities_for_nodes_plus_bottom_in_this_section[[jj]] = conditional_likelihoods_at_branch_section_bottom / sum(conditional_likelihoods_at_branch_section_bottom)
				
				
				# If you are storing ALL of the conditional likelihoods that were calculated
				if ((return_condlikes_table == TRUE) || (calc_TTL_loglike_from_condlikes_table == TRUE))
					{
					# Find the row in the big conditional likelihoods table
					TF1 = inputs$master_table$stratum == i
					TF2 = inputs$master_table$piecenum == jj
					TF3 = inputs$master_table$piececlass == "subbranch"
					TF = (TF1 + TF2 + TF3) == 3
				
					rownum = (1:nrow(condlikes_table))[TF]
					condlikes_table[rownum, ] = conditional_likelihoods_at_branch_section_bottom
					}
				
				
				} else {
				############################################
				# Otherwise, treepiece is a subtree
				############################################
				tmp_subtree = treepiece


				# Check for fossils
				# If you are storing ALL of the conditional likelihoods that were calculated
				if ((return_condlikes_table == TRUE) || (calc_TTL_loglike_from_condlikes_table == TRUE))
					{
					tmp_subtree_tipnums = 1:length(tmp_subtree$tip.label)
					for (iter in 1:length(tmp_subtree_tipnums))
						{
						# Find the row in the master table corresponding to this subtree_tip
						subtree_tip = tmp_subtree_tipnums[iter]
						
						TF1 = inputs$master_table$stratum == i
						TF2 = inputs$master_table$piecenum == jj
						TF3 = inputs$master_table$piececlass == "subtree"
						TF4 = inputs$master_table$SUBnode == subtree_tip
						TF = (TF1 + TF2 + TF3 + TF4) == 4

						# Find the row
						rownum = (1:nrow(inputs$master_table))[TF]
						tmp_master_table_row = inputs$master_table[rownum, ]

						# Error check
						if (nrow(tmp_master_table_row) != 1)
							{
							stoptxt = paste("\n\nFATAL ERROR in stratified loglike calculation at i=", i, "; jj=", jj, "; ", 
									  'inputs$master_table$piececlass == "subtree"', "; subtree_tip=", subtree_tip, 
									  "\nnrow(tmp_master_table_row) should =1 but instead =", nrow(tmp_master_table_row), "\n", sep="")
							stop(stoptxt)
							}

						# Now check if it's a fossil that appears in this time bin
						master_tip_time_bp = tmp_master_table_row$time_bp
						time_top = tmp_master_table_row$time_top
						time_bot = tmp_master_table_row$time_bot
						is_fossil = tmp_master_table_row$fossils

						# If this is TRUE, there's a match and the fossil tip appears in this time period
						if ( (master_tip_time_bp >= time_top) && (master_tip_time_bp < time_bot) && (is_fossil == TRUE))
							{
							# Shorten the branchlength by master_tip_time_bp-time_top
							amount_to_shorten_by = master_tip_time_bp-time_top
							
							# Find the branch of the subtree!
							tmp2_edgeTF = tmp_subtree$edge[,2] == subtree_tip
							tmp2_edgenum = (1:nrow(tmp_subtree$edge))[tmp2_edgeTF]
							
							# Edit the length of the branch on this subtree tip
							tmp_subtree$edge.length[tmp2_edgenum] = tmp_subtree$edge.length[tmp2_edgenum] - amount_to_shorten_by
							# do_exponentiation = TRUE	# not needed here
							}
						} # end forloop through subtree tips
					} # End fossils check
				# That should be it, everything else works as normal
				# Except, do up-pass also


				
				# Get the names of the tips in this subtree
				tipnames = tmp_subtree$tip.label
				
				# Use the tipnames to get the conditional likelihoods at these tips
				tips_for_subtree_TF = phy_as_it_is_chopped_down$tip.label %in% tipnames
				subtree_tip_relative_probs_of_each_state = current_tip_relative_probs_of_each_state[tips_for_subtree_TF,states_to_use_TF]
				
				# Check if this subtree contains a fixed internal node on the master tree
				tmp_fixnode = NULL		# Default
				tmp_fixlikes = NULL		# Default
				if (!is.null(fixnode))
					{
					# e.g.
					# fixnode=20
					TF1 = inputs$master_table$node == fixnode
					TF2 = inputs$master_table$SUBnode.type == "root"
					TF3 = inputs$master_table$SUBnode.type == "internal"
					TF = ((TF1 + TF2 + TF3) == 2)
					tmprow = inputs$master_table[TF,]
					
					# Check if we're in the right stratum / piece / piececlass
					TF1 = tmprow$stratum == i
					TF2 = tmprow$piecenum == jj
					TF3 = tmprow$piececlass == "subtree"
					
					TF = ((TF1 + TF2 + TF3) == 3)
					if (TF == TRUE)
						{
						#txt = paste("Master tree node ", fixnode, " matched to i=", i, "; jj=", jj, "; piececlass=", piececlass, "; subtree subnode=", tmprow$SUBnode, sep="")
						#print(txt)
						#print(fixlikes)
						
						# Determine the number of the subnode in the subtree
						tmp_fixnode = tmprow$SUBnode
						tmp_fixlikes = fixlikes						
						}
					}
				
				
				# Calculate the likelihoods for this subtree
				calc_loglike_sp_results = calc_loglike_sp(
					tip_condlikes_of_data_on_each_state=subtree_tip_relative_probs_of_each_state, 
					phy=tmp_subtree, 
					Qmat=Qmat_tmp, 
					spPmat=NULL,
					return_what="all",
					probs_of_states_at_root=NULL,
					rootedge=TRUE,
					sparse=FALSE, 
					printlevel=printlevel,
					use_cpp=TRUE,
					input_is_COO=FALSE,
					spPmat_inputs=spPmat_inputs,
					cppSpMethod=cppSpMethod,
					cluster_already_open=cluster_already_open,
					calc_ancprobs=calc_ancprobs,	 # If TRUE, get e.g. relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS
					null_range_allowed=null_range_allowed,
					fixnode=tmp_fixnode,
					fixlikes=tmp_fixlikes,
					stratified=TRUE,		# This makes calc_loglike_sp skip UPPASS probs, which are irrelevant inside stratified analyses
					states_allowed_TF=states_allowed_TF
					)
	
				#chainsaw_result$conditional_likelihoods_at_branch_section_bottom[[jj]] = 
				
				# Also, store the conditional likelihoods for all nodes in this subtree
				# MINUS THE GODDAMN TIPS OF THE SUBTREE, THESE ARE ALREADY IN THERE
				tmp_tipnums = 1:length(tipnames)
				
				#tmp_tr_table = prt(tmp_subtree, printflag=FALSE, get_tipnames=FALSE)
				
				
				# If you are storing ALL of the conditional likelihoods that were calculated
				if ((return_condlikes_table == TRUE) || (calc_TTL_loglike_from_condlikes_table == TRUE))
					{
					for (rownum in 1:nrow(calc_loglike_sp_results$condlikes_of_each_state))
						{
						tmp_condlikes = calc_loglike_sp_results$condlikes_of_each_state[rownum,]
						
						subtree_node = rownum
						
						TF1 = inputs$master_table$stratum == i
						TF2 = inputs$master_table$piecenum == jj
						TF3 = inputs$master_table$piececlass == "subtree"
						TF4 = inputs$master_table$SUBnode == subtree_node
						TF = (TF1 + TF2 + TF3 + TF4) == 4
		
						condlikes_table_rownum = (1:nrow(condlikes_table))[TF]
						condlikes_table[condlikes_table_rownum, ] = tmp_condlikes
						
						if (calc_ancprobs == TRUE)
							{
							# Skip this, for the bottom of the root branch
							if (rownum <= nrow(calc_loglike_sp_results$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS))
								{
								# We also need the tip nodes in the subtree
								
								# Get relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS
								tmp = calc_loglike_sp_results$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[rownum,]
								relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_table[condlikes_table_rownum,] = tmp
								}
							if (rownum > nrow(calc_loglike_sp_results$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS))
								{
								# Get relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS
								tmp = NA
								relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_table[condlikes_table_rownum,] = tmp
								}

							}
						}
					}
				
				
				chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]] = matrix(data=calc_loglike_sp_results$condlikes_of_each_state[-tmp_tipnums, ], ncol=ncol(calc_loglike_sp_results$condlikes_of_each_state))
				
				# Matrix of tip likelihoods to delete so you don't repeat using them in the total
				# loglike
				#tiplikes_to_delete[[jj]] = calc_loglike_sp_results$condlikes_of_each_state[tmp_tipnums, ]
				
				# Relative probabilities -- all nodes plus branch bottom (jjust branch bottom, here)
				chainsaw_result$relative_probabilities_for_nodes_plus_bottom_in_this_section[[jj]] = calc_loglike_sp_results$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[-tmp_tipnums, ]
	
				# Relative probabilities -- just the new tip
				chainsaw_result$relative_probs_of_each_state_at_bottom_of_root_branch[[jj]] = calc_loglike_sp_results$relative_probs_of_each_state_at_bottom_of_root_branch


				# ONLY for the nodes in original tree, store the condlikes
				# Add these to the overall list of conditional likelihoods
				numrows_to_add = nrow(chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]])
				# also remove rootedge prob (fixbug)
				rownum_for_bottom_of_root = nrow(chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]])

				
				startrow = current_condlikes_row + 1
				endrow = current_condlikes_row + numrows_to_add
				all_relative_probs_of_each_state[startrow:endrow, states_to_use_TF] = chainsaw_result$relative_probabilities_for_nodes_plus_bottom_in_this_section[[jj]]
				
				all_condlikes_of_each_state[startrow:endrow, states_to_use_TF] = chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]]
				# fixbug, except for root
				#if (i != num_iterations)
				#	{
					#all_condlikes_of_each_state[endrow, states_to_use_TF] = matrix(data=0, nrow=1, ncol=sum(states_to_use_TF))
				#	}

				current_condlikes_row = current_condlikes_row + numrows_to_add
				
				} # End if/then on branch vs. subtree
	
			# Also, store the relative probabilities for the new tip
			#new_tip_likelihoods[jj, states_to_use_TF] = chainsaw_result$relative_probs_of_each_state_at_bottom_of_root_branch[[jj]]
			new_tip_likelihoods[jj, states_to_use_TF] = chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]][nrow(chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]]), ]
	
			
# Add these to the overall list of conditional likelihoods
# 			numrows_to_add = nrow(chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]])
# 			
# 			startrow = current_condlikes_row + 1
# 			endrow = current_condlikes_row + numrows_to_add
# 			all_relative_probs_of_each_state[startrow:endrow, states_to_use_TF] = chainsaw_result$relative_probabilities_for_nodes_plus_bottom_in_this_section[[jj]]
# 			
# 			all_condlikes_of_each_state[startrow:endrow, states_to_use_TF] = chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]]
	
			#print(chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]])
			#print(log(sum(chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]])))
			
			#rowSums(all_condlikes_of_each_state) != 0
			#tmp_all_condlikes_of_each_state = all_condlikes_of_each_state[rowSums(all_condlikes_of_each_state) != 0,]
			#currLnL = sum(log(rowSums(tmp_all_condlikes_of_each_state)))
			#cat("i=", i, "; jj=", jj, "; currLnL=", currLnL, "\n")
	
# 			current_condlikes_row = current_condlikes_row + numrows_to_add
	
	
			} # End loop through jj tree pieces WITHIN a stratum



		# Update for the next loop
		# Tip likelihoods
		current_tip_relative_probs_of_each_state = new_tip_likelihoods
	
		# Store previous round
		#old_phy_as_it_is_chopped_down = phy_as_it_is_chopped_down
		#old_chainsaw_result = chainsaw_result
		#old_new_tip_likelihoods = new_tip_likelihoods
		
		# Convey the tree to the next iteration
		phy_as_it_is_chopped_down = chainsaw_result$tree_to_chainsaw

		} # End loop through i strata

	# Remove rows that have not been filled (till zero)
	all_condlikes_of_each_state_zero_TF = all_condlikes_of_each_state == 0
	all_condlikes_of_each_state_nonzero_TF = all_condlikes_of_each_state_zero_TF == FALSE
	rows_that_are_numeric_zeros_TF = rowSums(all_condlikes_of_each_state_nonzero_TF) >= 1
	#rowSums(all_condlikes_of_each_state) != 0
	final_all_condlikes_of_each_state = all_condlikes_of_each_state[rows_that_are_numeric_zeros_TF,]
	
	#rowSums(all_relative_probs_of_each_state) != 0
	all_relative_probs_of_each_state = all_relative_probs_of_each_state[rows_that_are_numeric_zeros_TF,]
	
	#all_relative_probs_of_each_state
	# Note: LAGRANGE uses rootedge = TRUE
	# This is not the source of the bug...
	#rootedge=FALSE
	if (rootedge == TRUE)
		{
		grand_total_likelihood = sum(log(rowSums(final_all_condlikes_of_each_state)))
		grand_total_likelihood
		} else {
		# Skip the last row
		grand_total_likelihood = sum(log(rowSums(final_all_condlikes_of_each_state[-nrow(final_all_condlikes_of_each_state),])))
		grand_total_likelihood
		}
	#rootedge=TRUE
	
	# Check for NA -- this can be caused by e.g. dispersal matrix constraints of 0 causing NAs
	# in the calculation
	if (is.na(grand_total_likelihood) == TRUE)
		{
		TF = is.na(all_relative_probs_of_each_state[,1])
		tmpr = (1:nrow(all_relative_probs_of_each_state))[TF]
		
		stoptxt1 = paste("\n\nFATAL ERROR IN calc_loglike_sp_stratified(). grand_total_likelihood=NA.\n",
		"These rows of 'all_relative_probs_of_each_state' had NAs:\n",
		paste(tmpr, collpase=",", sep=""), "\n", 
		"\n",
		"One possible cause of this: your dispersal matrix may be too restrictive; try changing\n",
		"the 0 values to e.g. 0.0000001.  Good luck!", sep="")
		
		cat(stoptxt1)
		
		stop(stoptxt1)
		}
	
	
	if (calc_TTL_loglike_from_condlikes_table == TRUE)
		{
		# Standard LAGRANGE result (exactly)
		TF2 = inputs$master_table$SUBnode.type == "internal"
		TF3 = inputs$master_table$SUBnode.type == "orig_tip"	# These are the original tips likelihoods; doesn't matter for unambiguous tips, but
																# DOES matter if there is a detection model.
		TF4 = inputs$master_table$SUBnode.type == "root"
		TF234 = (TF2 + TF3 + TF4) == 1
		sum(TF234)
		
		TF = TF234 == 1
		nodes_in_original_tree = inputs$master_table[TF,]
		node_order_original = order(nodes_in_original_tree$node)
		
		# Get some output matrices
		condlikes_of_each_state = condlikes_table[TF,][node_order_original,]
		computed_likelihoods_at_each_node = rowSums(condlikes_of_each_state)
		grand_total_likelihood = sum(log(computed_likelihoods_at_each_node))
		
		

		if (calc_ancprobs == TRUE)
			{
			# Downpass relprobs at the branch bottoms BELOW the nodes (just above speciation events)
			# Can't use "root" nodes, they have no DOWNPASS BELOW NODE stored
			TF2 = inputs$master_table$SUBnode.type == "internal"
			TF3 = inputs$master_table$SUBnode.type == "tip"	# These are the original tips likelihoods; doesn't matter for unambiguous tips, but
																	# DOES matter if there is a detection model.
			TF4 = inputs$master_table$piececlass == "subtree"	# Take only subtree tips, or internal nodes; these should have the DOWNPASS below nodes stored
			TF234 = (TF2 + TF3 + TF4) == 2						# ==2, because we need BOTH subtree and subtree-internal-node or subtree-tip-node
			sum(TF234)
			
			TF = TF234 == 1

			nodes_in_original_tree = inputs$master_table[TF,]
			node_order_original = order(nodes_in_original_tree$node)

			tmptable = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_table[TF,][node_order_original,]
			relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS = tmptable / rowSums(tmptable)
			
			# This leaves out the root row, so add that in
			root_row = rep(NA, times=ncol(relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS))
			tmpmat1 = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[1:length(original_phy$tip.label), ]
			tmpmat3_rows = (length(original_phy$tip.label)+1):nrow(relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS)
			tmpmat3 = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[tmpmat3_rows, ]
			relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS = rbind(tmpmat1, root_row, tmpmat3)
			
			# Relative probability of states at nodes, at the branch tops, ON THE DOWNPASS
			# (needs to be recalculated, to be in the right order)
			relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS = NULL
			#tmpmat = matrix(data=computed_likelihoods_at_each_node, ncol=1)
			tmptable = condlikes_of_each_state / rowSums(condlikes_of_each_state)
			relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS = tmptable
			
			
			# Get the relative probabilities at the root
			anc_row_of_master_table_TF = inputs$master_table$node.type=="root"
			#anc_row_of_master_table = (1:nrow(inputs$master_table))[anc_row_of_master_table_TF]
			anc_node_original_tree = inputs$master_table$node[anc_row_of_master_table_TF]
			anc_node_original_tree
			
			# Just always use the root node, not anything below it!
			starting_probs = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[anc_node_original_tree, ]
			}
		}




	#######################################################
	# Now do UPPASS for internal nodes
	#######################################################
	if (calc_ancprobs == TRUE)
		{
		cat("\nUppass started for (STRATIFIED) marginal ancestral states estimation!\n", sep="")





		#######################################################
		#######################################################
		# THIS IS AN UPPASS FROM THE TIPS TO THE ROOT
		#######################################################
		#######################################################

		
		# Setup matrices
		numrows_for_UPPASS = original_phy$Nnode + length(original_phy$tip.label)
		relative_probs_of_each_state_at_branch_top_AT_node_UPPASS = matrix(data=NA, nrow=numrows_for_UPPASS, ncol=ncol(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS))
		relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS = matrix(data=NA, nrow=numrows_for_UPPASS, ncol=ncol(relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS))
		
		
		# Vist edges in reverse order from the downpass
		# This would only work on un-stratified trees
		#edges_to_visit_uppass = seq(from=(num_internal_nodes*2), by=-2, length.out=num_internal_nodes)

		
		# Get the starting probabilities at the root
		# THIS ASSUMES THE STARTING PROBS ARE AT THE ROOT NODE, NOT SOME STUPID BRANCH BELOW THE ROOT NODE
		starting_probs
		
		# Put this starting prob into the node
		relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[anc_node_original_tree,] = starting_probs
		
		# Go through strata in REVERSE order
		for (i in num_iterations:1)
			{
	
			#######################################################
			# Cut down the number of areas, by what is allowed
			# (it would be more efficient to do this once during setup, but probably no biggie)
			#######################################################
			# states_to_use_TF: states to use in Qmat, speciation models, etc.
			# states_allowed_TF: use this to zero out impossible ancestral states according to
			#                    areas_allowed matrix
			# 
			
			if ( (is.null(inputs$list_of_areas_allowed_mats) == FALSE))
				{
				areas_allowed_mat = inputs$list_of_areas_allowed_mats[[i]]
		
				states_allowed_TF = sapply(X=all_states_list, FUN=check_if_state_is_allowed, areas_allowed_mat)
				#states_to_use_TF = all_states_list %in% tmp_states_list
				
				if (null_range_allowed == TRUE)
					{
					states_allowed_TF[1] = TRUE
					}
				# NO; use all areas for this
				# states_to_use_TF = states_allowed_TF
				
				} else {
				# Make no change
				pass = 1
				#states_list = states_list
				states_allowed_TF = rep(TRUE, length(all_states_list))
				}
			# Use this for regular calculations (Qmat, speciation models, etc.)
			states_to_use_TF = rep(TRUE, length(all_states_list))

			
				
			#####################################################
			# Make the dedf matrix for this time period
			#####################################################
			# If there is a distance matrix, use the first one 
			# (non-stratified analysis, here)
	
			# If there is a distance matrix, use the first one 
			# (non-stratified analysis, here)
			if ( (is.null(inputs$list_of_distances_mats) == FALSE))
				{
				distances_mat = inputs$list_of_distances_mats[[i]]
				} else {
				# Default is all areas effectively equidistant
				distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))
				}
		
			# Get the exponent on distance, apply to distances matrix
			dispersal_multipliers_matrix = distances_mat ^ x
		
			# Apply manual dispersal multipliers, if any
			# If there is a manual dispersal multipliers matrix, use the first one 
			# (non-stratified analysis, here)
			if ( (is.null(inputs$list_of_dispersal_multipliers_mats) == FALSE))
				{
				manual_dispersal_multipliers_matrix = as.matrix(inputs$list_of_dispersal_multipliers_mats[[i]])
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
			if ( (is.null(inputs$list_of_area_of_areas) == FALSE))
				{
				area_of_areas = inputs$list_of_area_of_areas[[i]]
				} else {
				# Default is all areas effectively equidistant
				area_of_areas = rep(1, length(areas))
				}
				
			# Get the exponent on extinction, apply to extinction modifiers	
			extinction_modifier_list = area_of_areas ^ (1 * u)
			
			# Apply to extinction rate
			elist = extinction_modifier_list * rep(e, length(areas))
			
			
			
			
			
			# Calculate the Q matrix
			# someday we'll have to put "a" (anagenic range-switching) in here...
	# 		if (is.null(Qmat))
	# 			{
				Qmat_tmp = rcpp_states_list_to_DEmat(areas_list=allareas_list, states_list=all_states_list[states_to_use_TF], 
				dmat=dmat, elist=elist, amat=amat, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)
	# 			} else {
	# 			# If Qmat is pre-specified
	# 			Qmat_tmp = Qmat
	# 			}
			
			if (sparse == TRUE)
				{
				# Sparse matrix exponentiation
				original_Qmat = Qmat_tmp
				
				# number of states in the original matrix
				coo_n = ncol(Qmat_tmp)
				anorm = as.numeric(norm(original_Qmat, type="O"))
				matvec = original_Qmat
				
				# *DO* TRANSPOSE; we want to go FORWARDS in time, NOT BACKWARDS!
				tmatvec = base::t(matvec)
				tmatvec = matvec
				tmpQmat_in_REXPOKIT_coo_fmt = mat2coo(tmatvec)
				}
			
			
			# Now. IF you have a subtree structure, you need to run this with a cladogenesis matrix, 
			# through calc_loglike_sp(), like normal.
			
			# If there's just one tree, store it in the object
			if (is.null(inputs$timeperiods) || length(inputs$timeperiods) == 1)
				{
				tr = read.tree(inputs$trfn)
				
				tree_to_chainsaw = NULL
				tree_to_chainsaw[[1]] = tr
	
				return_pieces_list = NULL
				return_pieces_list[[1]] = tr
				
				return_pieces_basenames = NULL
				
				# Merge THEN split THEN sort!!
				tmp_labels_merge = paste(tr$tip.label, collapse=",", sep="")
				tmp_labels_split = strsplit(tmp_labels_merge, split=",")[[1]]
				return_pieces_basenames[[1]] = paste(sort(tmp_labels_split), collapse=",", sep="")
				
				chainsaw_object = list()
				chainsaw_object$tree_to_chainsaw = tree_to_chainsaw
				chainsaw_object$return_pieces_list = return_pieces_list
				chainsaw_object$return_pieces_basenames = return_pieces_basenames
				attr(chainsaw_object, "class") = "chainsaw_result"
				
				inputs$tree_sections_list[[1]] = chainsaw_object
				}
			
			
			# OK, if you have a tree here, do that
			# if not, exp the branch
			
			#######################################################
			# Cladogenic model 
			#######################################################
			j = BioGeoBEARS_model_object@params_table["j","est"]
			ysv = BioGeoBEARS_model_object@params_table["ysv","est"]
			ys = BioGeoBEARS_model_object@params_table["ys","est"]
			v = BioGeoBEARS_model_object@params_table["v","est"]
			y = BioGeoBEARS_model_object@params_table["y","est"]
			s = BioGeoBEARS_model_object@params_table["s","est"]
			sum_SPweights = y + s + j + v
		
			maxent_constraint_01 = BioGeoBEARS_model_object@params_table["mx01","est"]
			
			# Text version of speciation matrix	
			maxent_constraint_01v = BioGeoBEARS_model_object@params_table["mx01v","est"]
			#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
				
			# Set the parameter controlling the size distribution of 
			# the smaller descendant species
			maxent01s_param = BioGeoBEARS_model_object@params_table["mx01s","est"]
			maxent01v_param = BioGeoBEARS_model_object@params_table["mx01v","est"]
			maxent01j_param = BioGeoBEARS_model_object@params_table["mx01j","est"]
			maxent01y_param = BioGeoBEARS_model_object@params_table["mx01y","est"]
		
		
			# Cladogenesis model inputs
			spPmat_inputs = NULL
	
			# Note that this gets the dispersal multipliers matrix, which is applied to 
			# e.g. the j events, NOT the dmat above which is d*dispersal_multipliers_matrix
			spPmat_inputs$dmat = dispersal_multipliers_matrix
	
			states_indices = all_states_list[states_to_use_TF]
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
			
			# Store the states_list in "l"
			l = spPmat_inputs$l
			
			#######################################################
			# Calculate the speciation model, and put in COO_weights_columnar
			#######################################################
			#Old, seems bogus
			#numareas = max(unlist(spPmat_inputs$l), na.rm=TRUE) + 1
			numareas = max(sapply(X=spPmat_inputs$l, FUN=length), na.rm=TRUE) + 0
			
			maxent01s = relative_probabilities_of_subsets(max_numareas=numareas, maxent_constraint_01=maxent01s_param, NA_val=0)
			maxent01v = relative_probabilities_of_vicariants(max_numareas=numareas, maxent_constraint_01v=maxent01v_param, NA_val=0)
			maxent01j = relative_probabilities_of_subsets(max_numareas=numareas, maxent_constraint_01=maxent01j_param, NA_val=0)
			maxent01y = relative_probabilities_of_subsets(max_numareas=numareas, maxent_constraint_01=maxent01y_param, NA_val=0)

			# Matrix of probs for each ancsize
			maxprob_as_function_of_ancsize_and_decsize = mapply(FUN=max, maxent01s, maxent01v, maxent01j, maxent01y, MoreArgs=list(na.rm=TRUE))
			maxprob_as_function_of_ancsize_and_decsize = matrix(data=maxprob_as_function_of_ancsize_and_decsize, nrow=nrow(maxent01s), ncol=ncol(maxent01s))
			maxprob_as_function_of_ancsize_and_decsize[maxprob_as_function_of_ancsize_and_decsize > 0] = 1
			maxprob_as_function_of_ancsize_and_decsize[maxprob_as_function_of_ancsize_and_decsize <= 0] = 0
			
			# Now, go through, and make a list of the max minsize for each decsize
			max_minsize_as_function_of_ancsize = apply(X=maxprob_as_function_of_ancsize_and_decsize, MARGIN=1, FUN=maxsize)

			
			tmpca_1 = rep(1, sum(states_to_use_TF)-1)  # -1, assumes NULL range is allowed
			tmpcb_1 = rep(1, sum(states_to_use_TF)-1)	# -1, assumes NULL range is allowed
			COO_weights_columnar = rcpp_calc_anclikes_sp_COOweights_faster(Rcpp_leftprobs=tmpca_1, Rcpp_rightprobs=tmpcb_1, l=l, s=s, v=v, j=j, y=y, dmat=dmat, maxent01s=maxent01s, maxent01v=maxent01v, maxent01j=maxent01j, maxent01y=maxent01y, max_minsize_as_function_of_ancsize=max_minsize_as_function_of_ancsize, printmat=FALSE)


			# This gives 15 states
			Rsp_rowsums = rcpp_calc_rowsums_for_COOweights_columnar(COO_weights_columnar=COO_weights_columnar)

			cppSpMethod = 3
			
			# Check to make sure you have the necessary inputs
			if (exists("COO_weights_columnar") == FALSE)
				{
				stop("\nERROR_A: calc_loglike_sp requires 'COO_weights_columnar', 'Rsp_rowsums', and cppSpMethod==3 for marginal ancestral state estimations.\n")
				}
			if (exists("Rsp_rowsums") == FALSE)
				{
				stop("\nERROR_B: calc_loglike_sp requires 'COO_weights_columnar', 'Rsp_rowsums', and cppSpMethod==3 for marginal ancestral state estimations.\n")
				}
			if (cppSpMethod != 3)
				{
				stop("\nERROR_C: calc_loglike_sp requires 'COO_weights_columnar', 'Rsp_rowsums', and cppSpMethod==3 for marginal ancestral state estimations.\n")
				}



				
			#######################################################
			# Go through the tree pieces in this stratum
			#######################################################
			chainsaw_result = inputs$tree_sections_list[[i]]
			
			# Set up a new list item to store uppass tip probs
			inputs$tree_sections_list[[i]]$pieces_relprobs_at_tips = list()

	
			# Go through tree pieces in this stratum (bottom first)
			for (jj in 1:length(chainsaw_result$return_pieces_list))
				{
				treepiece = chainsaw_result$return_pieces_list[[jj]]

				#cat("\ni=", i, "; jj=",jj, "; length(treepiece)=", length(treepiece), sep="")
				
				# If it's just a branch section
				if (is.numeric(treepiece) )
					{
					do_exponentiation = TRUE	# default

					# Also, exclude the case where there is a branch at the bottom below the bottom root node
					if (i == num_iterations)
						{
						errortxt = "ERROR: In stratified analysis, your tree must start with a root node, not a branch below the root node."
						stop(errortxt)
						}
					
					
					# Get the length of this branch
					subbranch_length = treepiece
					
					
					# Get the anc node in the original tree
					TF1 = inputs$master_table$stratum == i
					TF2 = inputs$master_table$piecenum == jj
					TF3 = inputs$master_table$piececlass == "subbranch"
					TF = (TF1 + TF2 + TF3) == 3
					anc_node_original_tree =  inputs$master_table$node[TF]
					
	
					# Check for fossil
					# Find the row
					rownum = (1:nrow(inputs$master_table))[TF]
					tmp_master_table_row = inputs$master_table[rownum, ]
					
					# Error check
					if (nrow(tmp_master_table_row) != 1)
						{
						stoptxt = paste("\n\nFATAL ERROR in stratified loglike UPPASS calculation at i=", i, "; jj=", jj, 
								  "; ", 'inputs$master_table$piececlass == "subbranch"', 
						          "\nnrow(tmp_master_table_row) should =1 but instead =", nrow(tmp_master_table_row), "\n", sep="")
						stop(stoptxt)
						}
					
					# Now check if it's a fossil that appears in this time bin
					master_tip_time_bp = tmp_master_table_row$time_bp
					time_top = tmp_master_table_row$time_top
					time_bot = tmp_master_table_row$time_bot
					is_fossil = tmp_master_table_row$fossils
					
					# If this is TRUE, there's a match and the fossil tip appears in this time period
					if ( (master_tip_time_bp >= time_top) && (master_tip_time_bp < time_bot) && (is_fossil == TRUE))
						{
						# Shorten the branchlength by master_tip_time_bp-time_top
						amount_to_shorten_by = master_tip_time_bp-time_top
						subbranch_length = subbranch_length - amount_to_shorten_by
						do_exponentiation = TRUE
						}

					# If this is TRUE, this fossil hasn't occured yet, and you are looking at the "phantom limb".
					# In this case, DON'T do matrix exponentiation, just copy the likelihoods down!!
					if ( master_tip_time_bp < time_top )
						{
						do_exponentiation = FALSE
						}
					# If FALSE, you're below all this and hopefully don't care
					
					
					
					# Get the uppass probs from the right piece in the previous stratum
					previous_stratum = i + 1
					
					
					# Get the number of the previous treepiece
					previous_stratum_TF = inputs$master_table$stratum == previous_stratum
					node_TF = inputs$master_table$node == anc_node_original_tree
					TF = (previous_stratum_TF + node_TF) == 2
					master_table_row_corresponding_to_anctip = inputs$master_table[TF,]
					
					previous_treepiece_num = master_table_row_corresponding_to_anctip$piecenum


					#print("HEY 2")
					#print(previous_stratum)
					#print(anc_node_original_tree)
					#print(master_table_row_corresponding_to_anctip)
					#print(previous_treepiece_num)

					# The previous treepiece
					previous_treepiece = inputs$tree_sections_list[[previous_stratum]]$return_pieces_list[[previous_treepiece_num]]

					
					# Relprobs from previous treepiece
					relprobs_at_tips_of_anc_treepiece = inputs$tree_sections_list[[previous_stratum]]$pieces_relprobs_at_tips[[previous_treepiece_num]]
					relprobs_at_branch_bottoms_below_tips_from_previous_stratum = inputs$tree_sections_list[[previous_stratum]]$pieces_relprobs_at_bottoms_below_tips[[previous_treepiece_num]]
					
					# If ancestor was a sub-branch
					if (is.numeric(previous_treepiece) == TRUE)
						{
						ancprobs_at_subbranch_bottom = relprobs_at_tips_of_anc_treepiece
						ancprobs_at_bottom_of_total_branch = relprobs_at_branch_bottoms_below_tips_from_previous_stratum
						} else {
						# Ancestor was a  sub-tree
						# Which tip in the previous treepiece?
						tipnum_in_previous_treepiece = master_table_row_corresponding_to_anctip$SUBnode
						
						# Extract those relative probabilities
						ancprobs_at_subbranch_bottom = relprobs_at_tips_of_anc_treepiece[tipnum_in_previous_treepiece, ]
						ancprobs_at_bottom_of_total_branch = relprobs_at_branch_bottoms_below_tips_from_previous_stratum[tipnum_in_previous_treepiece, ]
						}					
					
					
					# Do the exponentiation, unless it's a "phantom limb"!
					if (do_exponentiation == TRUE)
						{
						# Then do a forward matrix exponentiation step
						# Do sparse or dense matrix exponentiation
						if (sparse==FALSE)
							{
							# Dense matrix exponentiation
							# Need to do a forward matrix exponentiation
							actual_probs_after_forward_exponentiation = calc_prob_forward_onebranch_dense(relprobs_branch_bottom=ancprobs_at_subbranch_bottom, branch_length=subbranch_length, Qmat_tmp)
							actual_probs_after_forward_exponentiation[1] = 0 	# NULL range is impossible
							actual_probs_after_forward_exponentiation = actual_probs_after_forward_exponentiation / sum(actual_probs_after_forward_exponentiation)
							} else {
	
							actual_probs_after_forward_exponentiation = calc_prob_forward_onebranch_sparse(relprobs_branch_bottom=ancprobs_at_subbranch_bottom, branch_length=subbranch_length, tmpQmat_in_REXPOKIT_coo_fmt, coo_n=coo_n, anorm=anorm, check_for_0_rows=TRUE)
							actual_probs_after_forward_exponentiation[1] = 0 	# NULL range is impossible
							actual_probs_after_forward_exponentiation = actual_probs_after_forward_exponentiation / sum(actual_probs_after_forward_exponentiation)
							}

						# Zero out impossible states in this zone (but NOT for "phantom limbs")
						if (!is.null(states_allowed_TF))
							{
							actual_probs_after_forward_exponentiation[states_allowed_TF==FALSE] = 0
							actual_probs_after_forward_exponentiation = actual_probs_after_forward_exponentiation / sum(actual_probs_after_forward_exponentiation)
							}	

						} else {
						# Just pass up the ancestral probabilities, without modification, on the "phantom limb"
						actual_probs_after_forward_exponentiation = ancprobs_at_subbranch_bottom
						}
	
	

					##########################################################
					# tip probabilities for next stratum up
					##########################################################
					relprobs_at_tips_for_next_stratum_up = actual_probs_after_forward_exponentiation
					relprobs_at_branch_bottoms_below_tips_for_next_stratum_up = ancprobs_at_bottom_of_total_branch
					
					# Store the relprobs at the tips, so that the next stratum up can access them...
					inputs$tree_sections_list[[i]]$pieces_relprobs_at_tips[[jj]] = relprobs_at_tips_for_next_stratum_up
					inputs$tree_sections_list[[i]]$pieces_relprobs_at_bottoms_below_tips[[jj]] = relprobs_at_branch_bottoms_below_tips_for_next_stratum_up
					
					# If this is a tip in the master tree, store the tmp_relprobs_at_branchtop_AT_node_UPPASS in the master tree
					# (do this after the pass through the whole tree)
					
					
					} else {
					# Otherwise, treepiece is a subtree
					tmp_subtree = treepiece

					# Check for fossils on the UPPASS, shorten branches appropriately if found
					tmp_subtree_tipnums = 1:length(tmp_subtree$tip.label)
					for (iter in 1:length(tmp_subtree_tipnums))
						{
						# Find the row in the master table corresponding to this subtree_tip
						subtree_tip = tmp_subtree_tipnums[iter]
						
						TF1 = inputs$master_table$stratum == i
						TF2 = inputs$master_table$piecenum == jj
						TF3 = inputs$master_table$piececlass == "subtree"
						TF4 = inputs$master_table$SUBnode == subtree_tip
						TF = (TF1 + TF2 + TF3 + TF4) == 4

						# Find the row
						rownum = (1:nrow(inputs$master_table))[TF]
						tmp_master_table_row = inputs$master_table[rownum, ]

						# Error check
						if (nrow(tmp_master_table_row) != 1)
							{
							stoptxt = paste("\n\nFATAL ERROR in stratified loglike UPPASS calculation at i=", i, "; jj=", jj, "; ", 
									  'inputs$master_table$piececlass == "subtree"', "; subtree_tip=", subtree_tip, 
									  "\nnrow(tmp_master_table_row) should =1 but instead =", nrow(tmp_master_table_row), "\n", sep="")
							stop(stoptxt)
							}

						# Now check if it's a fossil that appears in this time bin
						master_tip_time_bp = tmp_master_table_row$time_bp
						time_top = tmp_master_table_row$time_top
						time_bot = tmp_master_table_row$time_bot
						is_fossil = tmp_master_table_row$fossils
						
						# If this is TRUE, there's a match and the fossil tip appears in this time period
						if ( (master_tip_time_bp >= time_top) && (master_tip_time_bp < time_bot) && is_fossil == TRUE)
							{
							# Shorten the branchlength by master_tip_time_bp-time_top
							amount_to_shorten_by = master_tip_time_bp-time_top
							
							# Find the branch of the subtree!
							tmp2_edgeTF = tmp_subtree$edge[,2] == subtree_tip
							tmp2_edgenum = (1:nrow(tmp_subtree$edge))[tmp2_edgeTF]
							
							# Edit the length of the branch on this subtree tip
							tmp_subtree$edge.length[tmp2_edgenum] = tmp_subtree$edge.length[tmp2_edgenum] - amount_to_shorten_by
							# do_exponentiation = TRUE	# not needed here
							}
						} # end forloop through subtree tips
					 	  # End fossils check
						  # (Do on tmp_subtree BEFORE it's converted to phy2!!)


					# Reorder the subtree
					# This is CRUCIAL
					phy2 <- reorder(tmp_subtree, "pruningwise")

					# Get the names of the tips in this subtree
					tipnames = phy2$tip.label
					
					# Get the root node of the subtree (was visit edges in reverse order from the downpass)
					num_internal_nodes = phy2$Nnode
					edges_to_visit_uppass = seq(from=(num_internal_nodes*2), by=-2, length.out=num_internal_nodes)
					tmpj = edges_to_visit_uppass[1]
					tmpi = tmpj - 1
					anc <- phy2$edge[tmpi, 1]				


					# Temporary matrix for UPPASS probabilities
					numnodes = num_internal_nodes + length(tmp_subtree$tip.label)
					tmp_relprobs_at_branchtop_AT_node_UPPASS = matrix(data=NA, nrow=numnodes, length(states_to_use_TF))
					tmp_relprobs_at_branchbot_BELOW_node_UPPASS = matrix(data=NA, nrow=numnodes, length(states_to_use_TF))

					
					# Get the starting uppass probabilities
					#inputs$tree_sections_list[[2]]$return_pieces_list[[3]]$root.edge
					# Otherwise, it's just the ancprobs at the node

					# Get the node number of subtree root node in the original tree
					# Get anc_node_original_tree
					TFi = inputs$master_table$stratum == i
					TFjj = inputs$master_table$piecenum == jj
					TF_SUBnode = inputs$master_table$SUBnode == anc
					TF = ((TFi + TFjj + TF_SUBnode) == 3)
					anc_node_original_tree = inputs$master_table$node[TF]

					# i.e., if i==5 in the Psychotria dataset
					if (i == num_iterations) # You are at the bottom tree piece, just use root node
						{
						#ancprobs_at_subtree_root = starting_probs

						# Anc node of the subtree
						tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, ] = relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[anc_node_original_tree, ]
						
						# None of this, at the root
						# NO: tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc, ] = ancprobs_at_bottom_of_total_branch
						
						} else {
						# You are NOT at the bottom tree piece
						if ((is.numeric(phy2$root.edge) == TRUE) && (!is.null(phy2$root.edge)) && (phy2$root.edge > 0))
							{
							# Get the length of this branch
							root_edge_length = phy2$root.edge
							
							# Get the uppass probs from the right piece in the previous stratum
							previous_stratum = i + 1
							
							
							# Get the number of the previous treepiece
							previous_stratum_TF = inputs$master_table$stratum == previous_stratum
							node_TF = inputs$master_table$node == anc_node_original_tree
							TF = (previous_stratum_TF + node_TF) == 2
							master_table_row_corresponding_to_anctip = inputs$master_table[TF,]
												
							previous_treepiece_num = master_table_row_corresponding_to_anctip$piecenum
							

							# The previous treepiece
							previous_treepiece = inputs$tree_sections_list[[previous_stratum]]$return_pieces_list[[previous_treepiece_num]]
							
							# Relprobs from previous treepiece
							relprobs_at_tips_of_anc_treepiece = inputs$tree_sections_list[[previous_stratum]]$pieces_relprobs_at_tips[[previous_treepiece_num]]
							
							#cat("\ni=", i, "; jj=", jj, "; previous_stratum=", previous_stratum, "; previous_treepiece_num=", previous_treepiece_num, "\n", sep="")
							#print(inputs$tree_sections_list[[previous_stratum]]$pieces_relprobs_at_bottoms_below_tips)
							relprobs_at_branch_bottoms_below_tips_from_previous_stratum = inputs$tree_sections_list[[previous_stratum]]$pieces_relprobs_at_bottoms_below_tips[[previous_treepiece_num]]

							
							
							# If ancestor was a sub-branch
							if (is.numeric(previous_treepiece) == TRUE)
								{
								ancprobs_at_subbranch_bottom = relprobs_at_tips_of_anc_treepiece
								ancprobs_at_bottom_of_total_branch = relprobs_at_branch_bottoms_below_tips_from_previous_stratum
								} else {
								# Ancestor was a  sub-tree
								# Which tip in the previous treepiece?
								tipnum_in_previous_treepiece = master_table_row_corresponding_to_anctip$SUBnode
								
								# Extract those relative probabilities
								ancprobs_at_subbranch_bottom = relprobs_at_tips_of_anc_treepiece[tipnum_in_previous_treepiece, ]
								ancprobs_at_bottom_of_total_branch = relprobs_at_branch_bottoms_below_tips_from_previous_stratum[tipnum_in_previous_treepiece, ]
								}	


							# Then do a forward matrix exponentiation step
							# Do sparse or dense matrix exponentiation
							if (sparse==FALSE)
								{
								# Dense matrix exponentiation
								# Need to do a forward matrix exponentiation
								actual_probs_after_forward_exponentiation = calc_prob_forward_onebranch_dense(relprobs_branch_bottom=ancprobs_at_subbranch_bottom, branch_length=root_edge_length, Qmat_tmp)
								actual_probs_after_forward_exponentiation[1] = 0 	# NULL range is impossible
								actual_probs_after_forward_exponentiation = actual_probs_after_forward_exponentiation / sum(actual_probs_after_forward_exponentiation)
								} else {
								# Sparse matrix exponentiation
								actual_probs_after_forward_exponentiation = calc_prob_forward_onebranch_sparse(relprobs_branch_bottom=ancprobs_at_subbranch_bottom, branch_length=root_edge_length, tmpQmat_in_REXPOKIT_coo_fmt, coo_n=coo_n, anorm=anorm, check_for_0_rows=TRUE)
								actual_probs_after_forward_exponentiation[1] = 0 	# NULL range is impossible
								actual_probs_after_forward_exponentiation = actual_probs_after_forward_exponentiation / sum(actual_probs_after_forward_exponentiation)
								}

							# Zero out impossible states
							if (!is.null(states_allowed_TF))
								{
								actual_probs_after_forward_exponentiation[states_allowed_TF==FALSE] = 0
								actual_probs_after_forward_exponentiation = actual_probs_after_forward_exponentiation / sum(actual_probs_after_forward_exponentiation)
								}	

							ancprobs_at_subtree_root = actual_probs_after_forward_exponentiation
							tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, ] = ancprobs_at_subtree_root
							tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc, ] = ancprobs_at_bottom_of_total_branch
							} else {
							# No root edge; just use probs at anc of subtree
							tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, ] = relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[anc_node_original_tree, ]
							tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc, ] = ancprobs_at_bottom_of_total_branch
							}
						
						}




					# Do dense matrix exponentiation of the subtree branches ahead of time
					if (sparse==FALSE)
							{
							# Get probmats for each branch, put into a big array
							# Create empty array to store results
							#independent_likelihoods_on_each_branch = array(0, dim=c(nrow(Qmat), ncol(Qmat), length(phy2$edge.length)))
							
							independent_likelihoods_on_each_branch = vector("list", length(phy2$edge.length))
							tmpmatrix = matrix(data=0, nrow=nrow(Qmat_tmp), ncol=ncol(Qmat_tmp))
							for (m in 1:length(phy2$edge.length))
								{
								independent_likelihoods_on_each_branch[[m]] = tmpmatrix
								}
							# Calculate the conditional likelihoods for each branch
							# dgexpv NOT ALLOWED when you have a null range state
							# (maybe try very very small values here)
				
							# clusterApply and other multicore stuff (e.g. doMC) are apparently dangerous on R.app
							if (!is.null(cluster_already_open))
								{
								# 
								if (.Platform$GUI == "AQUA")
									{
									cat("In calc_loglike_sp(), cluster_already_open=", cluster_already_open, " which means you want to calculate likelihoods on branches using a multicore option.\n", sep="")
									cat("But .Platform$GUI='AQUA', which means you are running the Mac GUI R.app version of R.  Parallel multicore functions, e.g. as accessed via \n", sep="")
									cat("library(parallel), are apparently dangerous/will crash R.app (google multicore 'R.app').  So, changing to cluster_already_open=NULL.\n", sep="")
									cluster_already_open=NULL
									}
								}
		
							# clusterApply etc. appear to NOT work on R.app
							if (!is.null(cluster_already_open))
								{
								# mcmapply
								#library(parallel)
								#independent_likelihoods_on_each_branch = mcmapply(FUN=expokit_dgpadm_Qmat, Qmat=list(Qmat), t=phy2$edge.length, transpose_needed=TRUE, SIMPLIFY="array", mc.cores=Ncores)
								independent_likelihoods_on_each_branch = clusterApply(cl=cluster_already_open, x=phy2$edge.length, fun=expokit_dgpadm_Qmat2, Qmat=Qmat_tmp, transpose_needed=TRUE)
								} else {
								# Not parallel processing
								#independent_likelihoods_on_each_branch = mapply(FUN=expokit_dgpadm_Qmat, Qmat=list(Qmat), t=phy2$edge.length, transpose_needed=TRUE, SIMPLIFY="array")
								independent_likelihoods_on_each_branch = mapply_likelihoods(Qmat_tmp, phy2, transpose_needed=TRUE)
								#independent_likelihoods_on_each_branch
								}
							}
					
			
			

					#######################################################
					# UPPASS loop here through the subtree
					#######################################################
					# paste here
					for (uj in edges_to_visit_uppass)		# Since we are going backwards
						{
						# First edge visited is ui
						#print(ui)
						
						# Its sister is uj 
						#uj <- ui - 1
						ui <- uj - 1		# Since we are going backwards
						
						# Get the node numbers at the tips of these two edges		
						left_desc_nodenum <- phy2$edge[ui, 2]
						right_desc_nodenum <- phy2$edge[uj, 2]
			
						# And for the ancestor edge (i or j shouldn't matter, should produce the same result!!!)
						anc <- phy2$edge[ui, 1]
						
						# get the correct edges
						left_edge_TF = phy2$edge[,2] == left_desc_nodenum
						right_edge_TF = phy2$edge[,2] == right_desc_nodenum
						
						# Check the branchlength of each edge
						# It's a hook if either branch is super-short
						is_leftbranch_hook_TF = phy2$edge.length[left_edge_TF] < min_branchlength
						is_rightbranch_hook_TF = phy2$edge.length[right_edge_TF] < min_branchlength
						hooknode_TF = (is_leftbranch_hook_TF + is_rightbranch_hook_TF) > 0
			
			
						#cat(i, j, left_desc_nodenum, right_desc_nodenum, hooknode_TF, "\n", sep="	")
						
						
						# You start with these uppass probs, for this node
						tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, ]

						
						# Apply speciation model to get the uppass probs at the base of the two descendant branches
						if (hooknode_TF == TRUE)
							{
							# Just copy the probs up, since a time-continuous model was assumed.
							tmp_relprobs_at_branchbot_BELOW_node_UPPASS[left_desc_nodenum, ] = tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, ]
							tmp_relprobs_at_branchbot_BELOW_node_UPPASS[right_desc_nodenum, ] = tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, ]
							# You're done
							
							} else {
							# Apply regular speciation model, with the weights given in COO_weights_columnar, and the 
							# normalization factor (sum of the weights across each row/ancestral state) in Rsp_rowsums.
							num_nonzero_split_scenarios = length(COO_weights_columnar[[1]])
							
							numstates = ncol(tip_condlikes_of_data_on_each_state)
							relprobs_just_after_speciation_UPPASS_Left = rep(0, numstates)
							relprobs_just_after_speciation_UPPASS_Right = rep(0, numstates)
							
							# Go through each ancestral state
							for (ii in 1:numstates)
								{
								ancprob = tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, ii]
								
								# A certain number of rows in COO_weights_columnar[[1]] will match this ancstate
								# (be sure to convert R 1-based index to C++ 0-based index
								ancstate_matches_TF = ((ii-1) == COO_weights_columnar[[1]])
								
								# For range inheritance scenarios which have this ancestor, get the 
								# Right and Left state indexes (1-based)
								if (null_range_allowed == TRUE)
									{
									# You have to add another 1, since the speciational models EXCLUDE
									# NJM 7/2013 -- no, + 0 works for constrained analysis with UPPASS
									# the first "range", the null range
									# NJM 2013-07-15 -- YES, do +1 or you end up with state #16 (KOMH) getting prob
									# 0 on the UPPASS
									Lstate_1based_indexes = COO_weights_columnar[[2]][ancstate_matches_TF] + 1 + 1
									Rstate_1based_indexes = COO_weights_columnar[[3]][ancstate_matches_TF] + 1 + 1
									} else {
									Lstate_1based_indexes = COO_weights_columnar[[2]][ancstate_matches_TF] + 1					
									Rstate_1based_indexes = COO_weights_columnar[[3]][ancstate_matches_TF] + 1
									}
							
								# And get the probability of each transition
								# AND multiply by the prob of this ancestor
								split_probs_for_this_ancestor = ancprob * COO_weights_columnar[[4]][ancstate_matches_TF]
								
								# Then add to the uppass probs for these branch bottoms
								relprobs_just_after_speciation_UPPASS_Left[Lstate_1based_indexes] = relprobs_just_after_speciation_UPPASS_Left[Lstate_1based_indexes] + split_probs_for_this_ancestor
								relprobs_just_after_speciation_UPPASS_Right[Rstate_1based_indexes] = relprobs_just_after_speciation_UPPASS_Right[Rstate_1based_indexes] + split_probs_for_this_ancestor
								
								} # That should be it for calculating the relative probs. Still have to normalize!

							# Zero out impossible states
							if (!is.null(states_allowed_TF))
								{
								relprobs_just_after_speciation_UPPASS_Left[states_allowed_TF==FALSE] = 0
								relprobs_just_after_speciation_UPPASS_Right[states_allowed_TF==FALSE] = 0
								}	

							
							# Normalize the probs by their sum.
							relprobs_just_after_speciation_UPPASS_Left = relprobs_just_after_speciation_UPPASS_Left / sum(relprobs_just_after_speciation_UPPASS_Left)				
							relprobs_just_after_speciation_UPPASS_Right = relprobs_just_after_speciation_UPPASS_Right / sum(relprobs_just_after_speciation_UPPASS_Right)
							
							# Store these uppass probs for the branch bases
							tmp_relprobs_at_branchbot_BELOW_node_UPPASS[left_desc_nodenum, ] = relprobs_just_after_speciation_UPPASS_Left
							tmp_relprobs_at_branchbot_BELOW_node_UPPASS[right_desc_nodenum, ] = relprobs_just_after_speciation_UPPASS_Right
							
							} # End if hooknode_TF
			
							
						# Finally, we have to retrieve the matrix exponentiations to calculate the probabilities
						# at the branch tops, from the probabilities at the branch bottoms.
				
						# Do sparse or dense matrix exponentiation
						if (sparse==FALSE)
							{
							# Dense matrix exponentiation, which has been done already!
							if (is.null(cluster_already_open))
								{
								# Relative probabilities of states at the top of left branch
								condprobs_Left_branch_top = tmp_relprobs_at_branchbot_BELOW_node_UPPASS[left_desc_nodenum,] %*% independent_likelihoods_on_each_branch[,,ui]
								condprobs_Left_branch_top[1] = 0	# zero out the NULL range, since it is impossible in a survivor
											
								# Relative probabilities of states at the top of right branch
								condprobs_Right_branch_top = tmp_relprobs_at_branchbot_BELOW_node_UPPASS[right_desc_nodenum,] %*% independent_likelihoods_on_each_branch[,,uj]
								condprobs_Right_branch_top[1] = 0	# zero out the NULL range, since it is impossible in a survivor
								} else {
								
								# Here, the independent_likelihoods_on_each_branch are stored in a list of matrices
								# Relative probabilities of states at the top of left branch
								condprobs_Left_branch_top = tmp_relprobs_at_branchbot_BELOW_node_UPPASS[left_desc_nodenum,] %*% independent_likelihoods_on_each_branch[[ui]]
								condprobs_Left_branch_top[1] = 0	# zero out the NULL range, since it is impossible in a survivor
											
								# Relative probabilities of states at the top of right branch
								condprobs_Right_branch_top = tmp_relprobs_at_branchbot_BELOW_node_UPPASS[right_desc_nodenum,] %*% independent_likelihoods_on_each_branch[[uj]]
								condprobs_Right_branch_top[1] = 0	# zero out the NULL range, since it is impossible in a survivor
								}
							
							} else {
							# Sparse matrix exponentiation
							# These are done on the fly, as the transition matrices cannot be stored, really
							
							# Left branch
							relprobs_branch_bottom = tmp_relprobs_at_branchbot_BELOW_node_UPPASS[left_desc_nodenum,]
							branch_length = phy2$edge.length[ui]
							
							condprobs_Left_branch_top = calc_prob_forward_onebranch_sparse(relprobs_branch_bottom, branch_length, tmpQmat_in_REXPOKIT_coo_fmt, coo_n=coo_n, anorm=anorm, check_for_0_rows=TRUE, TRANSPOSE_because_forward=TRUE)
							condprobs_Left_branch_top[1] = 0	# zero out the NULL range, since it is impossible in a survivor
							
							
							# Right branch
							relprobs_branch_bottom = tmp_relprobs_at_branchbot_BELOW_node_UPPASS[right_desc_nodenum,]
							branch_length = phy2$edge.length[uj]
							
							condprobs_Right_branch_top = calc_prob_forward_onebranch_sparse(relprobs_branch_bottom, branch_length, tmpQmat_in_REXPOKIT_coo_fmt, coo_n=coo_n, anorm=anorm, check_for_0_rows=TRUE, TRANSPOSE_because_forward=TRUE)
							condprobs_Right_branch_top[1] = 0	# zero out the NULL range, since it is impossible in a survivor
							}		
						
						# Zero out impossible states
						if (!is.null(states_allowed_TF))
							{
							condprobs_Left_branch_top[states_allowed_TF==FALSE] = 0
							condprobs_Right_branch_top[states_allowed_TF==FALSE] = 0
							
							condprobs_Left_branch_top = condprobs_Left_branch_top / sum(condprobs_Left_branch_top)
							condprobs_Right_branch_top = condprobs_Right_branch_top / sum(condprobs_Right_branch_top)
							}	
						
						
						# Normalize and save these probabilities
						tmp_relprobs_at_branchtop_AT_node_UPPASS[left_desc_nodenum,] = condprobs_Left_branch_top / sum(condprobs_Left_branch_top)
						tmp_relprobs_at_branchtop_AT_node_UPPASS[right_desc_nodenum,] = condprobs_Right_branch_top / sum(condprobs_Right_branch_top)
						
						} # End uppass loop through subtree



					# Store the UPPASS relprobs in the main matrix
					# Temporary matrix for UPPASS probabilities
					for (rownum in 1:nrow(tmp_relprobs_at_branchtop_AT_node_UPPASS))
						{
						# Store the relative probabilities at branch tops, WHEN
						# the internal nodes correspond to the master tree
						tmp_relprobs = tmp_relprobs_at_branchtop_AT_node_UPPASS[rownum,]
						
						subtree_node = rownum
						
						# We could use EITHER
						# (1) Internal nodes and tips of subtrees
						TF1 = inputs$master_table$stratum == i
						TF2 = inputs$master_table$piecenum == jj
						TF3 = inputs$master_table$piececlass == "subtree"
						TF4 = inputs$master_table$SUBnode == subtree_node
						TF5a = inputs$master_table$node.type == "internal"	# Store only if node corresponds to a node or tip in orig_tree
						TF5b = inputs$master_table$node.type == "tip"		# Store only if node corresponds to a node or tip in orig_tree
						TF_subtrees = (TF1 + TF2 + TF3 + TF4 + TF5a + TF5b) == 5
						
						# (2) Tips of subbranches that are also tips of the master tree
						# (see below)
						
						TF = TF_subtrees
						
						
						# Store in the FINAL table
						relative_probs_of_each_state_at_branch_top_AT_node_UPPASS_rownum = inputs$master_table$node[TF]
						relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[relative_probs_of_each_state_at_branch_top_AT_node_UPPASS_rownum, ] = tmp_relprobs
						# the above works
						
						
						#print(relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS)
						#print(dim(relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS))
						#print(relative_probs_of_each_state_at_branch_top_AT_node_UPPASS_rownum)
						#print(tmp_relprobs)
						
						# Storing the UPPass relative probabilities at branch bottoms doesn't work so well, probably
						# since the branch bottom is disconnected from the top sometimes.
						
						# You WANT to pick the subtree tips when they are cut off by a stratum boundary, in the 
						# case where you are doing node bottoms.
						# inputs$master_table[75:100,]
						#TF5a = inputs$master_table$node.type == "internal"	# Store only if node corresponds to a node or tip in orig_tree
						#TF5b = inputs$master_table$node.type == "tip"		# Store only if node corresponds to a node or tip in orig_tree
						#TF = (TF1 + TF2 + TF3 + TF4 + TF5a + TF5b) == 5

						
						tmp_relprobs = tmp_relprobs_at_branchbot_BELOW_node_UPPASS[rownum,]						
						relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS[relative_probs_of_each_state_at_branch_top_AT_node_UPPASS_rownum, ] = tmp_relprobs
						}

					##########################################################
					# tip probabilities for next stratum up
					##########################################################
					relprobs_at_tips_for_next_stratum_up = tmp_relprobs_at_branchtop_AT_node_UPPASS[1:length(phy2$tip.label), ]
					
					# Branch bottoms below tips can also be transferred up
					relprobs_at_branch_bottoms_below_tips_for_next_stratum_up = tmp_relprobs_at_branchbot_BELOW_node_UPPASS[1:length(phy2$tip.label), ]

					# Store the relprobs at the tips, so that the next stratum up can access them...
					inputs$tree_sections_list[[i]]$pieces_relprobs_at_tips[[jj]] = relprobs_at_tips_for_next_stratum_up
					inputs$tree_sections_list[[i]]$pieces_relprobs_at_bottoms_below_tips[[jj]] = relprobs_at_branch_bottoms_below_tips_for_next_stratum_up
					
					} # End if/then on branch vs. subtree

		
				} # End loop through jj tree pieces WITHIN a stratum
	
			} # End loop through i strata

		# (2) Tips of subbranches that are also tips of the master tree
		# (see below)
		# Tips of the master tree
		tipnums_of_master_tree = 1:length(original_phy$tip.labe)
		for (tn in 1:length(tipnums_of_master_tree))
			{
			# Find the row of the master table
			TF1 = inputs$master_table$piececlass == "orig_tip"
			TF2 = inputs$master_table$node == tn
			TF = ((TF1 + TF2) == 2)
			
			tmprow = inputs$master_table[TF,]
			
			# Now find the subbranch or subtree tip that corresponds to this tip (which could be living or extinct)
			TF1 = inputs$master_table$node == tmprow$node
			TF2 = inputs$master_table$time_top == tmprow$time_top
			TF3 = inputs$master_table$piececlass != "orig_tip"
			TF = ((TF1 + TF2 + TF3) == 3)
			tmprow2 = inputs$master_table[TF,]
			
			tmp_stratum = tmprow2$stratum
			tmp_piecenum = tmprow2$piecenum
			
			if (tmprow2$piececlass == "subtree")
				{
				tmp_tipnum = tmprow2$SUBnode
				tmp_tipprobs_at_top_UPPASS = inputs$tree_sections_list[[tmp_stratum]]$pieces_relprobs_at_tips[[tmp_piecenum]][tmp_tipnum, ]
				}
			if (tmprow2$piececlass == "subbranch")
				{
				tmp_tipnum = tmprow2$SUBnode
				tmp_tipprobs_at_top_UPPASS = inputs$tree_sections_list[[tmp_stratum]]$pieces_relprobs_at_tips[[tmp_piecenum]]
				}
			
			# Store these tip uppass probs in the final uppass probs
			relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[tn,] = tmp_tipprobs_at_top_UPPASS
			# The branch bottoms should be fine
			}

		cat("\nUppass completed for (STRATIFIED) marginal ancestral states estimation!\n", sep="")
		} # End UPPASS calculations, if calc_ancprobs == TRUE

	
	#return(grand_total_likelihood)
	# If you are storing ALL of the conditional likelihoods that were calculated
	if ((return_condlikes_table == TRUE) || (calc_TTL_loglike_from_condlikes_table == TRUE))
		{
		# Return an object with the condlikes_table, AND the grand conditional likelihood
		calc_loglike_sp_stratified_results = NULL
		calc_loglike_sp_stratified_results$final_all_condlikes_of_each_state = final_all_condlikes_of_each_state
		calc_loglike_sp_stratified_results$condlikes_table = condlikes_table
		
		if (calc_ancprobs == TRUE)
			{
			calc_loglike_sp_stratified_results$relative_probs_of_each_state_at_branch_top_AT_node_UPPASS = relative_probs_of_each_state_at_branch_top_AT_node_UPPASS
			calc_loglike_sp_stratified_results$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS = relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS
			
			# tmpres$relative_probs_of_each_state_at_branch_top_AT_node_UPPASS
			# tmpres$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS
	
	
			#######################################################
			# For branch bottoms
			#######################################################
			ML_marginal_prob_each_state_at_branch_bottom_below_node = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS * relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS
	
			# print
			ML_marginal_prob_each_state_at_branch_bottom_below_node
			rowSums(ML_marginal_prob_each_state_at_branch_bottom_below_node)
			
			ML_marginal_prob_each_state_at_branch_bottom_below_node = ML_marginal_prob_each_state_at_branch_bottom_below_node / rowSums(ML_marginal_prob_each_state_at_branch_bottom_below_node)
	
			# print
			ML_marginal_prob_each_state_at_branch_bottom_below_node
			sum_MLs_bot = rowSums(ML_marginal_prob_each_state_at_branch_bottom_below_node)
			sum_MLs_bot
			
			# Check for NaN rows
			NaN_TF = is.nan(sum_MLs_bot)
			numNaNs = sum(NaN_TF)
			numNaNs
			
			if (numNaNs > 0)
				{
				nannodenums = (1:length(NaN_TF))[NaN_TF]
				nannodenums_txt = paste(nannodenums, collapse=", ", sep="")
				txt = paste("\n\nWARNING! ML marginal states at branch bottoms produced ", numNaNs, " NaNs for nodes:\n", 
				nannodenums_txt, "\n", 
				"This probably means your downpass probabilities resulted in all 0 probabilities for the node.\n",
				"This might occur in a highly constrained model.\n",
				"As a fix, the downpass probabilities are being used for those nodes.\n", sep="")
				cat(txt)

				ML_marginal_prob_each_state_at_branch_bottom_below_node[NaN_TF,] = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[NaN_TF,]
				}
	
	
	
			#######################################################
			# For branch tops
			#######################################################
			ML_marginal_prob_each_state_at_branch_top_AT_node = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS * relative_probs_of_each_state_at_branch_top_AT_node_UPPASS
			
			# print
			ML_marginal_prob_each_state_at_branch_top_AT_node
			rowSums(ML_marginal_prob_each_state_at_branch_top_AT_node)
			
			ML_marginal_prob_each_state_at_branch_top_AT_node = ML_marginal_prob_each_state_at_branch_top_AT_node / rowSums(ML_marginal_prob_each_state_at_branch_top_AT_node)
			
			# print
			ML_marginal_prob_each_state_at_branch_top_AT_node
			sum_MLs_top = rowSums(ML_marginal_prob_each_state_at_branch_top_AT_node)
			sum_MLs_top
			
			# Check for NaN rows
			NaN_TF = is.nan(sum_MLs_top)
			numNaNs = sum(NaN_TF)
			numNaNs
			
			if (numNaNs > 0)
				{
				nannodenums = (1:length(NaN_TF))[NaN_TF]
				nannodenums_txt = paste(nannodenums, collapse=", ", sep="")
				txt = paste("\n\nWARNING! ML marginal states at branch tops produced ", numNaNs, " NaNs for nodes:\n", 
				nannodenums_txt, "\n", 
				"This probably means your downpass probabilities resulted in all 0 probabilities for the node.\n",
				"This might occur in a highly constrained model.\n",
				"As a fix, the downpass probabilities are being used for those nodes.\n", sep="")
				cat(txt)
				
				ML_marginal_prob_each_state_at_branch_top_AT_node[NaN_TF,] = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[NaN_TF,]
				}
			
			
			# Save them
			calc_loglike_sp_stratified_results$ML_marginal_prob_each_state_at_branch_bottom_below_node = ML_marginal_prob_each_state_at_branch_bottom_below_node
			calc_loglike_sp_stratified_results$ML_marginal_prob_each_state_at_branch_top_AT_node = ML_marginal_prob_each_state_at_branch_top_AT_node
	
			# tmpres$ML_marginal_prob_each_state_at_branch_top_AT_node
			# tmpres$ML_marginal_prob_each_state_at_branch_bottom_below_node
		}

	
		calc_loglike_sp_stratified_results$grand_total_likelihood = grand_total_likelihood

		if (calc_ancprobs == TRUE)
			{
			calc_loglike_sp_stratified_results$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS
			
			calc_loglike_sp_stratified_results$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS
			
			}
		
		#return(calc_loglike_sp_stratified_results)
		#return(condlikes_table)
		}
		
	if (return_what == "loglike")
		{
		return(grand_total_likelihood)
		}
	
	if (return_what == "all")
		{
		return(calc_loglike_sp_stratified_results)
		}

	# Just return the grand_total_likelihood (default)
	return(grand_total_likelihood)
	}








#######################################################
# calc_loglike_for_optim_stratified
#######################################################
#' Take model parameters and the data and calculate the log-likelihood -- stratified version
#' 
#' This is the stratified version of \code{\link{calc_loglike_for_optim}}. This function is an input to optim or optimx, the ML
#' estimation routines.
#' 
#' @param params A vector of parameters for optimization.
#' @param BioGeoBEARS_run_object Object containing the run parameters, and the model.
#' @param phy An ape tree object
#' @param tip_condlikes_of_data_on_each_state Conditional likelihoods at tips. A numeric 
#' matrix with rows representing tips, and columns representing states/geographic ranges.  The cells
#' give the likelihood of the observation data under the assumption that the tip has that 
#' state; typically this means that the known geographic range gets a '1' and all 
#' other states get a 0.
#' @param print_optim If TRUE (default), print the optimization steps as ML estimation progresses.
#' @param areas_list A list of the desired area names/abbreviations/letters (?).
#' @param states_list A list of the possible states/geographic ranges, in 0-based index form.
#' @param force_sparse Should sparse matrix exponentiation be used? Default \code{FALSE}.
#' @param cluster_already_open The cluster object, if it has already been started.
#' @return \code{ttl_loglike} The log-likelihood of the data under the input model and parameters.
#' @export
#' @seealso \code{\link[stats]{convolve}} chainsaw_result
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
#' 
calc_loglike_for_optim_stratified <- function(params, BioGeoBEARS_run_object, phy, tip_condlikes_of_data_on_each_state, print_optim=TRUE, areas_list, states_list, force_sparse=FALSE, cluster_already_open=FALSE)
	{
	# Put the parameters into the BioGeoBEARS_model_object, so that they can be universally read out
	# into any function
	BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
	
	# Put the parameters into the BioGeoBEARS_model_object, so that they can be universally read out
	# into any function
	BioGeoBEARS_model_object = params_into_BioGeoBEARS_model_object(BioGeoBEARS_model_object=BioGeoBEARS_model_object, params=params)
	
	# Update linked parameters
	BioGeoBEARS_model_object = calc_linked_params_BioGeoBEARS_model_object(BioGeoBEARS_model_object)
	
	# Set the dispersal and extinction rate
# 	d = BioGeoBEARS_model_object@params_table["d","est"]
# 	e = BioGeoBEARS_model_object@params_table["e","est"]

	#######################################################
	#######################################################
	# Do branch-length exponentiation if desired
	#######################################################
	#######################################################
	# Don't do this here, it would screw up sectioning of the tree!
# 	b = BioGeoBEARS_model_object@params_table["b","est"]
# 	# Modify the edge.lengths
# 	phy$edge.length = phy$edge.length ^ b


	#######################################################
	#######################################################
	# Do distance-dependence and dispersal multipliers matrix
	#######################################################
	#######################################################
	# Equal dispersal in all directions (unconstrained)
	# Equal extinction probability for all areas
	areas = areas_list
# 	distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

	#dmat = matrix(d, nrow=length(areas), ncol=length(areas))
	#elist = rep(e, length(areas))
	
	# Set up the instantaneous rate matrix (Q matrix)
	#Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)

	# Cladogenic model
# 	j = BioGeoBEARS_model_object@params_table["j","est"]
# 	ysv = BioGeoBEARS_model_object@params_table["ys","est"]
# 	v = BioGeoBEARS_model_object@params_table["v","est"]
# 	ys = BioGeoBEARS_model_object@params_table["ys","est"]
# 	sum_SPweights = ys + j + v
# 	sum_SPweights
# 	
	# Store back in the run object
	BioGeoBEARS_run_object$BioGeoBEARS_model_object = BioGeoBEARS_model_object
	
# 	maxent_constraint_01 = BioGeoBEARS_model_object@params_table["mx01","est"]
# 	
# 	# Text version of speciation matrix	
# 	maxent_constraint_01v = BioGeoBEARS_model_object@params_table["mx01v","est"]
# 	#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
# 		
# 	# Set the parameter controlling the size distribution of 
# 	# the smaller descendant species
# 	maxent01s_param = BioGeoBEARS_model_object@params_table["mx01s","est"]
# 	maxent01v_param = BioGeoBEARS_model_object@params_table["mx01v","est"]
# 	maxent01j_param = BioGeoBEARS_model_object@params_table["mx01j","est"]
# 	maxent01y_param = BioGeoBEARS_model_object@params_table["mx01y","est"]
# 
# 
# 	# Cladogenesis model inputs
# 	spPmat_inputs = NULL
# 	states_indices = states_list
# 	states_indices[1] = NULL	# shorten the states_indices by 1 (cutting the null range state from the speciation matrix)
# 	spPmat_inputs$l = states_indices
# 	spPmat_inputs$s = ys
# 	spPmat_inputs$v = v
# 	spPmat_inputs$j = j
# 	spPmat_inputs$y = ys
# 	spPmat_inputs$dmat = distances_mat
# 	spPmat_inputs$maxent01s_param = maxent01s_param
# 	spPmat_inputs$maxent01v_param = maxent01v_param
# 	spPmat_inputs$maxent01j_param = maxent01j_param
# 	spPmat_inputs$maxent01y_param = maxent01y_param
# 

	if (print_optim == TRUE)
		{
		#outvars = as.data.frame(t(BioGeoBEARS_model_object@params_table$est))
		#names(outvars) = rownames(BioGeoBEARS_model_object@params_table)
		#outvars = c(BioGeoBEARS_model_object@params_table$est)
		
		#cat("\n")
		#cat(outvars, sep="	")
		
		# Before calculating the log likelihood, print it, in case there is e.g. a bug
		#cat("d=", d, "; e=", e, "; j=", j, "; ys=", ys, "; v=", v, "; maxent01=", maxent_constraint_01, "; maxent01v=", maxent_constraint_01v, "; sum=", sum_SPweights, "; LnL=", sep="")
		}

	# Calculate the log-likelihood of the data, given the model parameters during this iteration	
# 	ttl_loglike = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="loglike", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, cluster_already_open=cluster_already_open)
# 	ttl_loglike
	#tmpinputs = BioGeoBEARS_run_object

	fixnode = BioGeoBEARS_run_object$fixnode
	fixlikes = BioGeoBEARS_run_object$fixlikes

	return_condlikes_table = BioGeoBEARS_run_object$return_condlikes_table
	calc_TTL_loglike_from_condlikes_table = BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table

	ttl_loglike = calc_loglike_sp_stratified(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=NULL, spPmat=NULL, min_branchlength=1e-21, return_what="loglike", probs_of_states_at_root=NULL, rootedge=TRUE, sparse=FALSE, printlevel=0, use_cpp=TRUE, input_is_COO=FALSE, spPmat_inputs=NULL, cppSpMethod=3, cluster_already_open=NULL, calc_ancprobs=FALSE, null_range_allowed=TRUE, fixnode=fixnode, fixlikes=fixlikes, inputs=BioGeoBEARS_run_object, allareas=areas, all_states_list=states_list, return_condlikes_table=return_condlikes_table, calc_TTL_loglike_from_condlikes_table=calc_TTL_loglike_from_condlikes_table)

	if (print_optim == TRUE)
		{
		LnL = ttl_loglike
		# If the log likelihood is successful, print it
		outvars = adf(t(c(BioGeoBEARS_model_object@params_table$est, LnL)))
		#outvars = cbind(outvars, LnL)
		
		names(outvars) = c(rownames(BioGeoBEARS_model_object@params_table), "LnL")
		print(round(outvars,3))

		#cat(ttl_loglike, "\n", sep="")
		}
	
	return(ttl_loglike)
	}







