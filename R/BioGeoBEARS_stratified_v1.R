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
#' @param min_dist_between_node_and_stratum_line An error check is run, if any nodes are
#' closer to a stratum boundary than this line, an error is thrown. The easiest 
#' solution is to change the date of your stratum boundary line slightly.
#' @param remove_root_edge Default TRUE, which means the root edge will be removed; 
#' chainsaw2 function will not work if it is present.
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
section_the_tree <- function(inputs, make_master_table=FALSE, plot_pieces=TRUE, cut_fossils=TRUE, fossils_older_than=0.1, min_dist_between_node_and_stratum_line=0.00001, remove_root_edge=TRUE, save_phys_before_they_are_chopped=FALSE)
	{
	runjunk='
	make_master_table=TRUE; plot_pieces=FALSE; cut_fossils=TRUE; fossils_older_than=0.6;
	
	inputs=BioGeoBEARS_run_object
	make_master_table=TRUE
	plot_pieces=TRUE
	cut_fossils=TRUE
	fossils_older_than=0.1
	min_dist_between_node_and_stratum_line=0.00001
	# testing:
	orig_timeperiods = c(0.5, 1.9, 3.7, 5.1)
	save_phys_before_they_are_chopped = FALSE
	' # END runjunk
	
	#	trstr = "((((((((P_hawaiiensis_WaikamoiL1:0.9656850499,P_mauiensis_Eke:0.9656850499):0.7086257935,(P_fauriei2:1.230218511,P_hathewayi_1:1.230218511):0.4440923324):0.1767115552,(P_kaduana_PuuKukuiAS:1.851022399,P_mauiensis_PepeAS:1.851022399):0.0008897862802):0.3347375986,P_kaduana_HawaiiLoa:2.185759997):0.302349378,(P_greenwelliae07:1.131363255,P_greenwelliae907:1.131363255):1.35674612):1.689170274,((((P_mariniana_MauiNui:1.994011054,P_hawaiiensis_Makaopuhi:1.994011054):0.7328279804,P_mariniana_Oahu:2.726839034):0.2574151709,P_mariniana_Kokee2:2.984254205):0.4601084855,P_wawraeDL7428:3.444362691):0.732916959):0.7345185743,(P_grandiflora_Kal2:2.479300491,P_hobdyi_Kuia:2.479300491):2.432497733):0.2873119899,((P_hexandra_K1:2.363984189,P_hexandra_M:2.363984189):0.4630447802,P_hexandra_Oahu:2.826939991):2.372081244);"
	#	tr = read.tree(file="", text=trstr)
	# END runjunk

	if (save_phys_before_they_are_chopped == TRUE)
		{
		phys_before_they_are_chopped = list()
		pcount = 0
		}
	
	
	# Fixing nodes, for marginal local optimum ancestral state reconstruction, is COMPLICATED when you are
	# chopping up an APE tree.  Somehow we would have to keep track of which node.  So, save this for later.
	
	
	orig_timeperiods = inputs$timeperiods
	timeperiods = orig_timeperiods
	#original_tree = read.tree(inputs$trfn)
	original_tree = check_trfn(trfn=inputs$trfn)
	
	
	
	
	# Remove root edge
	if (remove_root_edge == TRUE)
		{
		if ("root.edge" %in% names(original_tree))
			{
			txt = paste0("WARNING in chainsaw2: input tree had a 'root.edge', which crashes chainsaw2. Setting original_tree$root.edge=NULL.")
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			warning(txt)
			original_tree$root.edge = NULL
			} # END if ("root.edge" %in% names(original_tree))
		} # END if (remove_root_edge == TRUE)


	
	phy_as_it_is_chopped_down = original_tree
	
	# Make the tree table for the original tree
	orig_tr_table = prt(original_tree, printflag=FALSE, get_tipnames=TRUE)
	times_older_than_root_node_TF = orig_timeperiods > max(orig_tr_table$node_ht)
	times_younger_than_root_node_TF = orig_timeperiods < max(orig_tr_table$node_ht)
	
	
	# Error check
	if (sum(times_older_than_root_node_TF) >= 2)
		{
		txt = paste0("STOP ERROR in section_the_tree(): the timeperiods file can have ONLY ONE time older than the bottom node in your tree.  Use the function prt() to get a table of node ages for your tree. The oldest node age in your tree is: ", max(orig_tr_table$node_ht), ". The times in your timeperiods file are: ", paste(orig_timeperiods, collapse=" ", sep=""), ".")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		} # END if (sum(times_older_than_root_node_TF) >= 2)

	# Error check
	if (sum(times_younger_than_root_node_TF) == 0)
		{
		txt = paste0("STOP ERROR in section_the_tree(): the timeperiods file *HAS* to have an oldest time that is older than the bottom node in your tree.  Use the function prt() to get a table of node ages for your tree. The oldest node age in your tree is: ", max(orig_tr_table$node_ht), ". The times in your timeperiods file are: ", paste(orig_timeperiods, collapse=" ", sep=""), ".")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		} # END if (sum(times_older_than_root_node_TF) >= 2)
		

	#######################################################
	# Check that only ONE time is older than the root node of the tree
	#######################################################
	
	

	
	
	#######################################################
	# Check for the case where node(s) exactly match stratum boundaries
	# This would be BAD, so disallow it.
	#######################################################
	for (i in 1:length(orig_timeperiods))
		{
		TF = orig_tr_table$time_bp == orig_timeperiods[i]
		if (sum(TF) > 0)
			{
			errortxt = paste("\n\nERROR in section_the_tree(): your tree has ", sum(TF), " nodes with date ", orig_timeperiods[i], ".\nThis is a problem because you have a timeperiod boundary of: ", orig_timeperiods[i], "\nThe function doesn't know how to section a tree exactly at a node boundary.", "\nTo fix: change the timeperiod date, or edit the tree so that all nodes are more than ", min_dist_between_node_and_stratum_line, " time units from a timeperiod boundary\n(specified by 'min_dist_between_node_and_stratum_line', default min_dist_between_node_and_stratum_line=", min_dist_between_node_and_stratum_line, ").", "\n\nIf it makes you feel better, there is no way your dating of either phylogenetic or geological events is all that precise anyway.", sep="")
			cat(errortxt)

			cat("\n\nNodes with this problem:\n\n")
			print(orig_tr_table[TF,])
			
			stop("\n\nStopping on error.\n\n")
			}

		# Or, check for nodes too near to boundaries
		diffs = abs(orig_tr_table$time_bp - orig_timeperiods[i])
		TF = diffs < min_dist_between_node_and_stratum_line
		if (sum(TF) > 0)
			{
			errortxt = paste("\n\nERROR in section_the_tree(): your tree has ", sum(TF), " nodes with a date too close to your timeperiod boundary of: ", orig_timeperiods[i], ".\nThis is a problem because very short branches may cause issues with likelihood calculations, ancestral state estimation, and stochastic mapping.", "\nSee e.g. the min_branchlength option of calc_loglike_sp().", "\nTo fix: change the timeperiod date, or edit the tree so that all nodes are more than ", min_dist_between_node_and_stratum_line, " time units from a timeperiod boundary\n(specified by 'min_dist_between_node_and_stratum_line', default min_dist_between_node_and_stratum_line=", min_dist_between_node_and_stratum_line, ").", "\n\nIf it makes you feel better, there is no way your dating of either phylogenetic or geological events is all that precise anyway.", sep="")
			cat(errortxt)
			
			cat("\n\nNodes with this problem:\n\n")
			print(orig_tr_table[TF,])
			
			stop("\n\nStopping on error.\n\n")
			}
		}



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
							"But you have not turned on fossils by setting 'cut_fossils=FALSE' in section_the_tree().\n", 
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
							"If you actually have that many fossil tips, then everything is fine, and you can ignore this warning. If not, make sure that all fossils are older than whatever you set 'fossils_older_than' to be. If you do *not* have any fossils, then you are probably using an undated tree. This is a Very Bad Idea in general, please see 'BioGeoBEARS Mistakes To Avoid' at PhyloWiki.\n", 
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
			
			
			cat("\n", i, "- top: ", orig_timeperiods[i]-timeperiods[i], ", bot: ", orig_timeperiods[i], ", rel_bot: ", timeperiods[i], "\n", sep="")
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

			if (save_phys_before_they_are_chopped == TRUE)
				{
				phys_before_they_are_chopped[[(pcount=pcount+1)]] = phy_as_it_is_chopped_down
				}
		  
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
			
			# Store the chainsaw result (initial: store again if you change chainsaw_result)
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
				#print(i)
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
						
						
						################################################
						# 2016-02-29 bug fix
						# Fossil branches on sub-trees were not being cut down appropriately,
						# at least for hook nodes
						################################################
						# Identify tips that are fossils *in* the subtree:
						# these branches need to be cut down further
						actual_heights_below_bin_top = tmp_join_table$time_bp - tmp_join_table$time_top
						
						actual_height_lower_than_bin_top_TF = actual_heights_below_bin_top > 1e-10
						subtree_tip_TF = tmp_join_table$SUBnode.type == "tip"
						fossil_in_subtree_TF = (actual_height_lower_than_bin_top_TF + subtree_tip_TF) == 2
						
						# Declare them fossils WITHIN the subtree
						tmp_join_table$SUBfossils[fossil_in_subtree_TF] = TRUE
						tmp_join_table$SUBtime_bp[fossil_in_subtree_TF] = actual_heights_below_bin_top[fossil_in_subtree_TF]
						
						# Adjust the edge lengths in the table, and in the subtree
						# table
						new_subtree_edge_lengths = tmp_join_table$SUBedge.length[fossil_in_subtree_TF]
						tmp_join_table$SUBedge.length[fossil_in_subtree_TF] = tmp_join_table$SUBedge.length[fossil_in_subtree_TF] - actual_heights_below_bin_top[fossil_in_subtree_TF]
						
						# subtree
						# tmp_subtree = chainsaw_result$return_pieces_list[[p]]
						tmp_subtree2 = tmp_subtree
						#print(tmp_subtree2$edge.length)
						
						# Remove the root node edge length, which will NOT be in the 
						# tree object's list of edges
						subtree_tipnums_to_change = tmp_join_table$SUBnode[fossil_in_subtree_TF]
						# Match the subtree tipnums to the subtree's edge table
						edge_table_rownums_to_change = match(x=subtree_tipnums_to_change, table=tmp_subtree2$edge[,2])
						
						tmp_subtree2$edge.length[edge_table_rownums_to_change] = tmp_join_table$SUBedge.length[fossil_in_subtree_TF]
						#print(tmp_subtree2$edge.length)
												
						#print(tmp_join_table)
						#print(write.tree(tmp_subtree2, file=""))
						#plot(tmp_subtree2)
						#axisPhylo()
						
						
						chainsaw_result$return_pieces_list[[p]] = tmp_subtree2
						
						# Store the chainsaw result (again: store again if you change chainsaw_result)
						tree_sections_list[[tnum]] = chainsaw_result
						################################################
						# END 2016-02-29 bug fix
						################################################

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
						
						# END subtree
						# END if (classes_of_pieces[p] == "subtree")
						} else {
						# START sub-branch
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
						} # END if (classes_of_pieces[p] == "subtree")
					
					master_table = rbind(master_table, tmp_join_table)
					} # END for (p in 1:length(classes_of_pieces))
				} # END if (make_master_table == TRUE)
				


	
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


	if (save_phys_before_they_are_chopped == TRUE)
		{
		inputs$phys_before_they_are_chopped = phys_before_they_are_chopped
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
#' @param remove_root_edge Default TRUE, which means the root edge will be removed; 
#' chainsaw2 function will not work if it is present.
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
chainsaw2 <- function(tr, timepoint=10, return_pieces=TRUE, remove_root_edge=TRUE)
	{
	# Take a tree and saw it off evenly across a certain timepoint.
	# This removes any tips above the timepoint, and replaces them 
	# with a single tip representing the lineage crossing
	# the timepoint (with a new tip name).
	
	defaults='
	phy_as_it_is_chopped_down
	timepoint=timepoint
	return_pieces=TRUE
	'
	
	# Remove root edge
	if (remove_root_edge == TRUE)
		{
		if ("root.edge" %in% names(tr))
			{
			txt = paste0("WARNING in chainsaw2: input tree had a 'root.edge', which crashes chainsaw2. Setting tr$root.edge=NULL.")
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			warning(txt)
			tr$root.edge = NULL
			} # END if ("root.edge" %in% names(tr))
		} # END if (remove_root_edge == TRUE)
	
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
	
	
	# NJM 2014-12-11: this assumes your tree has NO root edge length;
	# I'll put in a check for this.Ffunction 
	NA_false = is.not.na(tree_to_chainsaw_table$edge.length)
	
	tree_to_chainsaw$edge.length[parent_branches[NA_false]] = tree_to_chainsaw_table$edge.length[NA_false]
	
	#######################################################
	# Error check
	#######################################################
	tmp_trtable = prt(tree_to_chainsaw, printflag=FALSE)
	brlens = tmp_trtable$edge.length
	TF = brlens <= 0
	TF[is.na(TF)] = FALSE
	if (sum(TF) > 0)
		{
		nodenums = (1:nrow(tmp_trtable))[TF]
		nodenums
		txt = paste0("STOP ERROR in chainsaw2(): the post-chainsaw tree had ", sum(TF), " negative branchlengths. READ THE FOLLOWING ERROR MESSAGE SLOWLY AND CAREFULLY AND YOU MAY FIND A SOLUTION. This error seems to sometimes occur with large cuts on trees with many fossil tips (i.e., non-contemporaneous tips). I'm not sure what causes the bug, except that chainsaw-ing an APE phylo object is quite complex, and it is even more complex for a tree with many non-contemporaneous tips. Imagine a phylogeny made of cardboard, then cutting it at various timepoints, then keeping track of all of the pieces.  Anyway, until I fix it, your best bet is to chainsaw2() in stages, using smaller cuts than the one that caused the error.  Or if you are using BioGeoBEARS and doing a time-stratified analysis, create extra time-strata (perhaps repeating the same settings for several time bins), so that the usage of chainsaw2() within section_the_tree() does not cause an error.")
		
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		} # END if (sum(TF) > 0)
	
	
	
	
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
#' @param include_null_range Does the state space include the null range?
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
calc_loglike_sp_stratified <- function(tip_condlikes_of_data_on_each_state, phy, Qmat=NULL, spPmat=NULL, min_branchlength=0.000001, return_what="loglike", probs_of_states_at_root=NULL, rootedge=TRUE, sparse=FALSE, printlevel=0, use_cpp=TRUE, input_is_COO=FALSE, spPmat_inputs=NULL, cppSpMethod=3, cluster_already_open=NULL, calc_ancprobs=FALSE, include_null_range=TRUE, fixnode=NULL, fixlikes=NULL, inputs=inputs, allareas=allareas, all_states_list=all_states_list, return_condlikes_table=FALSE, calc_TTL_loglike_from_condlikes_table=TRUE, m=NULL)
	{
	defaults='
	Qmat=NULL; spPmat=NULL; min_branchlength=0.000001; return_what="loglike"; probs_of_states_at_root=NULL; rootedge=FALSE; sparse=FALSE; printlevel=1; use_cpp=TRUE; input_is_COO=FALSE; spPmat_inputs=NULL; cppSpMethod=3; cluster_already_open=NULL; calc_ancprobs=FALSE; include_null_range=TRUE; fixnode=NULL; fixlikes=NULL; inputs=inputs; allareas=allareas; all_states_list=all_states_list; return_condlikes_table=FALSE; calc_TTL_loglike_from_condlikes_table=TRUE; m=NULL
	'
	
	defaults='
	maxareas = 4
	include_null_range = TRUE
	phy = read.tree(inputs$trfn)
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=np(inputs$geogfn))
	tip_condlikes_of_data_on_each_state = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges, phy, maxareas=maxareas, include_null_range=include_null_range)
	
	allareas = getareas_from_tipranges_object(tipranges)
	all_states_list = rcpp_areas_list_to_states_list(areas=allareas, include_null_range=TRUE, maxareas=maxareas)
	
	tmpres = calc_loglike_sp_stratified(tip_condlikes_of_data_on_each_state, phy, Qmat=NULL, spPmat=NULL, min_branchlength=0.000001, return_what="all", probs_of_states_at_root=NULL, rootedge=TRUE, sparse=FALSE, printlevel=0, use_cpp=TRUE, input_is_COO=FALSE, spPmat_inputs=NULL, cppSpMethod=3, cluster_already_open=NULL, calc_ancprobs=FALSE, include_null_range=TRUE, fixnode=NULL, fixlikes=NULL, inputs=inputs, allareas=allareas, all_states_list=all_states_list, return_condlikes_table=FALSE)
	tmpres
	
	min_branchlength=0.000001
	include_null_range=TRUE
	printlevel=0
	cppSpMethod=3
	return_condlikes_table=TRUE
	calc_TTL_loglike_from_condlikes_table=TRUE
	calc_ancprobs=TRUE
'
	defaults='
	tmpres = calc_loglike_sp_stratified(tip_condlikes_of_data_on_each_state, phy, Qmat=NULL, spPmat=NULL, min_branchlength=0.000001, return_what="all", probs_of_states_at_root=NULL, rootedge=TRUE, sparse=FALSE, printlevel=0, use_cpp=TRUE, input_is_COO=FALSE, spPmat_inputs=NULL, cppSpMethod=3, cluster_already_open=NULL, calc_ancprobs=TRUE, include_null_range=TRUE, fixnode=NULL, fixlikes=NULL, inputs=inputs, allareas=allareas, all_states_list=all_states_list, return_condlikes_table=TRUE, calc_TTL_loglike_from_condlikes_table=TRUE)
'



# 	defaults='
# 	# STANDARD DEBUGGING HERE
# 	tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state; phy=phy; Qmat=NULL; spPmat=NULL; min_branchlength=0.000001; return_what="loglike"; probs_of_states_at_root=NULL; rootedge=TRUE; sparse=FALSE; printlevel=0; use_cpp=TRUE; input_is_COO=FALSE; spPmat_inputs=NULL; cppSpMethod=3; cluster_already_open=NULL; calc_ancprobs=FALSE; include_null_range=TRUE; fixnode=fixnode; fixlikes=fixlikes; inputs=BioGeoBEARS_run_object; allareas=areas; all_states_list=states_list; return_condlikes_table=TRUE; calc_TTL_loglike_from_condlikes_table=TRUE;
# 	' # end junk


	# START OF FUNCTION
	BioGeoBEARS_run_object = inputs
	if (is.null(inputs$printlevel))
		{
		inputs$printlevel = 1
		}
	printlevel = inputs$printlevel



	# Is this a traits-based analysis?
	traitTF = is.null(BioGeoBEARS_run_object$trait) == FALSE
	if (traitTF == TRUE)
		{
		trait_Pmat_txt = BioGeoBEARS_run_object$trait_Pmat_txt
		num_trait_states = ncol(trait_Pmat_txt)
		} # END if (traitTF == TRUE)
	
	# Initialize m
	m = NULL
	# Initialize jts_matrix, matrix of t12, t23, etc., during a j event
	jts_matrix = NULL


	# Put the parameters into the BioGeoBEARS_model_object, so that they can be universally read out
	# into any function
	BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
	#print(params)
	#print(BioGeoBEARS_model_object)



	######################################################
	# 2016-03-23_NJM: adding rescaling
	# (unscale params, if they were used before)
	######################################################
	if (BioGeoBEARS_run_object$rescale_params == TRUE)
		{
		BioGeoBEARS_model_object@params_table = unscale_BGB_params(scaled_params_table=BioGeoBEARS_model_object@params_table)
		
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table = BioGeoBEARS_model_object@params_table
		}

	# Update linked parameters
	BioGeoBEARS_model_object = calc_linked_params_BioGeoBEARS_model_object(BioGeoBEARS_model_object)
	
	# Update to the run object, just to be SURE
	BioGeoBEARS_run_object$BioGeoBEARS_model_object = BioGeoBEARS_model_object
	inputs$BioGeoBEARS_model_object = BioGeoBEARS_model_object



	#######################################################
	# Error check on fixnode / fixlikes
	#######################################################
	if (!is.null(fixnode))
		{
		if (( is.null(dim(fixlikes)) == TRUE) && (length(fixnode)==1))
			{
			pass_fixlikes = TRUE
			} else {
			if ( (dim(fixlikes)[1]) == length(fixnode) )
				{
				pass_fixlikes = TRUE
				
				# Another error check: Multiple nodes in 'fixnode' MUST be sorted in increasing order
				if ( (all(c(order(fixnode) == 1:length(fixnode)))) == TRUE)
					{
					pass_fixlikes = TRUE
					} else {
					pass_fixlikes = FALSE
					error_msg = "ERROR in calc_loglike_sp_stratified(): \n             Multiple nodes in 'fixnode' MUST be sorted in increasing order.\n"
					cat(error_msg)
					stop(error_msg)
					}
				
				} else {
				pass_fixlikes = FALSE
				error_msg = "ERROR in calc_loglike_sp_stratified(): Either:\n             (1) fixnode must be a single node number, and fixlikes must be a vector, or\n             (2) fixlikes like must be a matrix with the # of rows equal to length(fixnode).\n"
				cat(error_msg)
				stop(error_msg)
				} # end 2nd if()
			} # end 1st if()
		}




	if ((return_condlikes_table == TRUE) || (calc_TTL_loglike_from_condlikes_table == TRUE))
		{
		names_in_inputs = names(inputs)			# can't use exists() on list items; reasons explained here:
												# http://stackoverflow.com/questions/7719741/how-to-test-if-list-element-exists
		if ( ("master_table" %in% names_in_inputs) == TRUE)
			{
			# Old
			#condlikes_table = matrix(data=0, nrow=nrow(inputs$master_table), ncol=length(all_states_list))
			# When traits are possible
			#condlikes_table = matrix(data=0, nrow=nrow(inputs$master_table), ncol=numstates_geogtrait)
			condlikes_table = matrix(data=0, nrow=nrow(inputs$master_table), ncol=ncol(tip_condlikes_of_data_on_each_state))
			
			# Put in the conditional likelihoods at the tips
			tmprownums = nrow(tip_condlikes_of_data_on_each_state)
			condlikes_table[1:tmprownums, ] = tip_condlikes_of_data_on_each_state
			
			if (calc_ancprobs == TRUE)
				{
				if (traitTF == FALSE)
					{
					relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE = matrix(data=0, nrow=nrow(inputs$master_table), ncol=length(all_states_list))
					} else {
					relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE = matrix(data=0, nrow=nrow(inputs$master_table), ncol=length(all_states_list)*num_trait_states)
					} # END if (traitTF == FALSE)
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
		num_timeperiods = 1
		} else {
		# Multiple timeperiods
		timeperiods = inputs$timeperiods
		num_timeperiods = length(timeperiods)
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
	all_relative_probs_of_each_state = matrix(0, ncol=ncol(tip_condlikes_of_data_on_each_state), nrow=(numnodes*length(timeperiods)))
	all_condlikes_of_each_state = matrix(0, ncol=ncol(tip_condlikes_of_data_on_each_state), nrow=(numnodes*length(timeperiods)))
	
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

	########################################################
	# DOWNPASS through the tree pieces
	########################################################
	for (i in 1:num_timeperiods)
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
# 		if (include_null_range == TRUE)
# 			{
# 			TF = all_states_list == "_"
# 			all_states_list[TF] = NA
# 			} else {
# 			TF = all_states_list == "_"
# 			all_states_list[TF] = NA			
# 			}

		#######################################################
		# Cut down the number of areas, by what is allowed
		# (it would be more efficient to do this once during setup, 
		# but probably no biggie)
		#######################################################
		# states_to_use_TF: states to use in Qmat, speciation models, etc.
		# states_allowed_TF: use this to zero out impossible ancestral states according to
		#                    areas_allowed matrix
		#######################################################
		
		# Should we modify the list of allowed states?
		# default: no areas_allowed or areas_adjacency constraints
		user_specified_constraints_on_states_list_TF = FALSE
		states_allowed_TF1 = rep(TRUE, length(all_states_list))
		states_allowed_TF2 = rep(TRUE, length(all_states_list))
		if ( (is.null(inputs$list_of_areas_allowed_mats) == FALSE))
			{
			user_specified_constraints_on_states_list_TF = TRUE
			}
		if ( (is.null(inputs$list_of_areas_adjacency_mats) == FALSE))
			{
			user_specified_constraints_on_states_list_TF = TRUE
			}


		# Get TF for whether each state in the master list is 
		# turned on in this time period.
		# (then edit Qmat etc.)
		if (user_specified_constraints_on_states_list_TF == TRUE)
			{
			# Check that lists_of_states_lists_0based has been specified
			if ( is.null(inputs$lists_of_states_lists_0based) == TRUE )
				{
				errortxt = paste0("STOP ERROR in calc_loglike_sp_stratified(): User has specified areas_allowed or area_adjacency constraints, but 'lists_of_states_lists_0based' has not been added to the BioGeoBEARS_run_object.")
				cat("\n\n")
				cat(errortxt)
				cat("\n\n")
				stop(errortxt)
				}
			
			
			# Areas allowed in this time bin
			if ( (is.null(inputs$list_of_areas_allowed_mats) == FALSE))
				{
				areas_allowed_mat = inputs$list_of_areas_allowed_mats[[i]]
	
				states_allowed_TF1 = sapply(X=all_states_list, FUN=check_if_state_is_allowed, areas_allowed_mat)
				#states_to_use_TF = all_states_list %in% tmp_states_list
			
				if (include_null_range == TRUE)
					{
					states_allowed_TF1[1] = TRUE
					}
				# NO; use all areas for this
				# states_to_use_TF = states_allowed_TF
				} # END if ( (is.null(inputs$list_of_areas_allowed_mats) == FALSE))
			
			# Areas adjacency
			if ( (is.null(inputs$list_of_areas_adjacency_mats) == FALSE))
				{
				areas_adjacency_mat = inputs$list_of_areas_adjacency_mats[[i]]
	
				states_allowed_TF2 = sapply(X=all_states_list, FUN=check_if_state_is_allowed_by_adjacency, areas_adjacency_mat)
				#states_to_use_TF = all_states_list %in% tmp_states_list
			
				if (include_null_range == TRUE)
					{
					states_allowed_TF2[1] = TRUE
					}
				# NO; use all areas for this
				# states_to_use_TF = states_allowed_TF
				} # END if ( (is.null(inputs$list_of_areas_adjacency_mats) == FALSE))
			# Combine the two (areas_allowed and areas_adjacency)
			states_allowed_TF = ((states_allowed_TF1 + states_allowed_TF2) == 2)
			} else {
			# Otherwise, 
			# make no change
			pass = 1
			#states_list = states_list
			states_allowed_TF = rep(TRUE, length(all_states_list))
			} # END if (user_specified_constraints_on_states_list_TF == TRUE)
		# Use this for regular calculations (Qmat, speciation models, etc.)
		states_to_use_TF = rep(TRUE, length(all_states_list))
		
		#print(states_to_use_TF)
		#print(states_allowed_TF)
		
			
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

		# Environmental distances
		if ( (is.null(inputs$list_of_envdistances_mats) == FALSE))
			{
			envdistances_mat = inputs$list_of_envdistances_mats[[i]]
			} else {
			# Default is all areas effectively equidistant
			envdistances_mat = matrix(1, nrow=length(areas), ncol=length(areas))
			}

		# Get the exponent on environmental distance, apply to distances matrix
		n = BioGeoBEARS_model_object@params_table["n","est"]
		dispersal_multipliers_matrix = dispersal_multipliers_matrix * envdistances_mat^n

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

		# Get the exponent on manual dispersal multipliers
		w = BioGeoBEARS_model_object@params_table["w","est"]
		
		#print("manual_dispersal_multipliers_matrix ^ w")
		#print(manual_dispersal_multipliers_matrix ^ w)
				
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
# 		if (is.null(Qmat))
# 			{
			
			# 2014 version
			#Qmat_tmp = rcpp_states_list_to_DEmat(areas_list=allareas_list, states_list=all_states_list[states_to_use_TF], 
			#dmat=dmat_times_d, elist=elist, amat=amat, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)
			# 2015 version
#			Qmat_tmp = rcpp_states_list_to_DEmat(areas_list=allareas_list, states_list=all_states_list[states_allowed_TF], dmat=dmat_times_d, elist=elist, amat=amat, include_null_range=include_null_range, normalize_TF=TRUE, makeCOO_TF=force_sparse)

			# 2018 version
			if (traitTF == FALSE)
				{
				Qmat_tmp = rcpp_states_list_to_DEmat(areas_list=allareas_list, states_list=all_states_list[states_allowed_TF], dmat=dmat_times_d, elist=elist, amat=amat, include_null_range=include_null_range, normalize_TF=TRUE, makeCOO_TF=force_sparse)

				#print(dim(Qmat_tmp))
	# 			} else {
	# 			# If Qmat is pre-specified
	# 			Qmat_tmp = Qmat
	# 			}
				} # END if (traitTF == FALSE)

		# Analysis with a trait modifying dispersal rate
		if (traitTF == TRUE)
			{
			num_geog_states = length(all_states_list[states_allowed_TF])
			numstates_geogtrait = num_trait_states * num_geog_states
			
# 			print("states_allowed_TF")
# 			print(states_allowed_TF)
# 			print("num_geog_states")
# 			print(num_geog_states)
# 			print("num_trait_states")
# 			print(num_trait_states)
# 			print("numstates_geogtrait")
# 			print(numstates_geogtrait)

			# Get the modified Qmatrix (traits + geog)			
			tmpres = modify_Qmat_with_trait(Qmat=NULL, BioGeoBEARS_run_object, numstates_geogtrait=numstates_geogtrait, areas_list=allareas_list, states_list=all_states_list[states_allowed_TF], dispersal_multipliers_matrix=dispersal_multipliers_matrix, elist=elist, force_sparse=force_sparse)
			Qmat_tmp = tmpres$Qmat
			m = tmpres$m
			
			
			# If the trait can change during jump events
			if (is.null(BioGeoBEARS_run_object$jts_txt_matrix) == FALSE)
				{
				jts_txt_matrix = BioGeoBEARS_run_object$jts_txt_matrix
				jts_matrix = matrix(data=0, nrow=nrow(jts_txt_matrix), ncol=ncol(jts_txt_matrix))
				TF_matrix = matrix(data=TRUE, nrow=nrow(jts_txt_matrix), ncol=ncol(jts_txt_matrix))
				diag(TF_matrix) = FALSE
				jts_txt_params = c(jts_txt_matrix[TF_matrix])
				jts_txt_params
			
				# Populate the numeric jts_matrix
				for (jts_i in 1:nrow(jts_txt_matrix))
					{
					diag_val = 1
					for (jts_j in 1:ncol(jts_txt_matrix))
						{
						if (jts_i == jts_j)
							{
							next()
							}
						jts_txt = jts_txt_matrix[jts_i,jts_j]
						newval = as.numeric(BioGeoBEARS_model_object@params_table[jts_txt, "est"])
						jts_matrix[jts_i,jts_j] = newval
						diag_val = 1-newval
						}
					# Populate the diagonal
					jts_matrix[jts_i,jts_i] = diag_val
					} # END for (jts_i in 1:nrow(jts_txt_matrix))
				} # END if (is.null(BioGeoBEARS_run_object$jts_txt_matrix) == FALSE)
			} else {
			num_geog_states = length(all_states_list[states_allowed_TF])
			numstates_geogtrait = num_geog_states
			} # END if (traitTF == TRUE)
		
		
		
		if (force_sparse == TRUE)
			{
			# Convert the COO-formatted trait+geog matrix to CRS format for kexpmv
			tmpQmat_in_REXPOKIT_coo_fmt = Qmat_tmp

			# Make a CRS-formatted matrix, for kexpmv
			# DO THE TRANSPOSE HERE, trait+geog matrices assembled transposed
			tmpQmat_in_kexpmv_crs_fmt = coo2crs(
				ia=tmpQmat_in_REXPOKIT_coo_fmt[,"ia"], 
				ja=tmpQmat_in_REXPOKIT_coo_fmt[,"ja"], 
				a =tmpQmat_in_REXPOKIT_coo_fmt[,"a"],
				n=numstates_geogtrait, transpose_needed=FALSE)
				
			Qmat_tmp = tmpQmat_in_REXPOKIT_coo_fmt
			} # END if (force_sparse == TRUE)

		
		
		# Now. IF you have a subtree structure, you need to run this with a cladogenesis matrix, 
		# through calc_loglike_sp(), like normal.
		
		# If there's just one tree, store it in the object
		if (is.null(inputs$timeperiods) || length(inputs$timeperiods) == 1)
			{
			#tr = read.tree(inputs$trfn)
			tr = check_trfn(trfn=inputs$trfn)
			
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
			} # END if (is.null(inputs$timeperiods) || length(inputs$timeperiods) == 1)
		
		
		# OK, if you have a tree here, do that
		# if not, exp the branch
		
		#######################################################
		# Cladogenic model 
		#######################################################
		spPmat_inputs = get_spPmat_inputs_from_BGB(BioGeoBEARS_run_object=BioGeoBEARS_run_object, states_list=all_states_list[states_allowed_TF], dispersal_multipliers_matrix=dispersal_multipliers_matrix)


	
		########################################################
		# DOWNPASS through the tree pieces
		#######################################################
		# Go through the tree pieces
		#######################################################
		chainsaw_result = inputs$tree_sections_list[[i]]
		
		# You will need the new tip likelihoods of the new tree:
		current_tip_relative_probs_of_each_state
		
		# Old
		#new_tip_likelihoods = matrix(0, nrow=length(chainsaw_result$return_pieces_list), ncol=length(all_states_list))
		# When traits are possible
		new_tip_likelihoods = matrix(0, nrow=length(chainsaw_result$return_pieces_list), ncol=ncol(current_tip_relative_probs_of_each_state))


		# Error check for traits model
		if (traitTF == TRUE)
			{
			# DOWNPASS definition of states_allowed_TF with traits
			wTrait_states_allowed_TF = c(rep(states_allowed_TF, times=num_trait_states))
			#print("wTrait_states_allowed_TF:")
			#print(wTrait_states_allowed_TF)
			if (sum(wTrait_states_allowed_TF) != numstates_geogtrait)
				{
				txt = paste0("STOP ERROR in calc_loglike_sp_stratified(): sum(wTrait_states_allowed_TF)=", sum(wTrait_states_allowed_TF), ", and numstates_geogtrait=", numstates_geogtrait, ". They must be equal to proceed.")
				cat("\n\n")
				cat(txt)
				cat("\n\n")
				stop(txt)
				} # END if (ncol(current_tip_relative_probs_of_each_state) != numstates_geogtrait)
			}


		for (jj in 1:length(chainsaw_result$return_pieces_list))
			{
			#cat("i=", i, " jj=",jj, "\n", sep="")
			treepiece = chainsaw_result$return_pieces_list[[jj]]
			treepiece
		
			############################################
			# DOWNPASS -- process just a branch section (an edge)
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
					this_row_of_master_table_is_being_used = TF
					
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
					
					
					# THIS (FOSSIL_HASNT_OCCURRED_YET) MUST GO *BEFORE* THE 'FOSSIL_TIP_DOES_OCCUR_IN_BIN' IF/THEN
					# If this is TRUE, this fossil hasn't occurred yet, and you are looking at the "phantom limb".
					# In this case, DON'T do matrix exponentiation, just copy the likelihoods down!!
					if (( master_tip_time_bp > time_top) && (is.na(is_fossil) == FALSE) && (is_fossil == TRUE))
						{
						do_exponentiation = FALSE
						}

					# THIS (FOSSIL_TIP_DOES_OCCUR_IN_BIN) MUST GO *AFTER* THE 'FOSSIL_HASNT_OCCURRED_YET' IF/THEN
					# If this is TRUE, there's a match and the fossil tip appears in this time period
					# (THIS IS CRUCIAL TO GETTING STRATIFICATION TO WORK -- YOU NEED THE is_fossil==TRUE ADDED!!)
					if ( (master_tip_time_bp >= time_top) && (master_tip_time_bp < time_bot) && (is.na(is_fossil) == FALSE) && (is_fossil == TRUE))
						{
						# Shorten the branchlength by master_tip_time_bp-time_top
						amount_to_shorten_by = master_tip_time_bp-time_top
						treepiece = treepiece - amount_to_shorten_by
						do_exponentiation = TRUE
						}
					
					# 2016-02-24
					# Also, DON'T do exponentiation if the branch length in the master branch
					# is a direct ancestor, i.e., less than min_branchlength
					if (tmp_master_table_row$edge.length < min_branchlength)
						{
						#print("It's a direct ancestor, so DON'T do matrix exponentiation!")
						# It's a direct ancestor, so DON'T do matrix exponentiation
						do_exponentiation = FALSE
						}
					
						
					# If FALSE, you're below all this and hopefully don't care
					} # END if ((return_condlikes_table == TRUE) || (calc_TTL_loglike_from_condlikes_table == TRUE))



				tipname = chainsaw_result$return_pieces_basenames[[jj]]
				tip_TF = phy_as_it_is_chopped_down$tip.label == tipname
				
				# 22 rather than 17
				relative_probs_of_each_state_at_the_tip_of_this_branch = current_tip_relative_probs_of_each_state[tip_TF, states_to_use_TF]
	
	
				if (do_exponentiation == TRUE)
					{
					# DENSE MATRIX EXPONENTIATION DOWN ONE BRANCH
					if (force_sparse == FALSE)
						{
						# t = treepiece
						independent_likelihoods_at_branch_section_bottom = expokit_dgpadm_Qmat2(times=treepiece, Qmat=Qmat_tmp,  transpose_needed=TRUE)
						#independent_likelihoods_at_branch_section_bottom = expokit_dgpadm_Qmat(Qmat=Qmat_tmp,  t=treepiece, transpose_needed=FALSE)
						
						
						if (traitTF == FALSE)
							{
							# 2014 version
							#conditional_likelihoods_at_branch_section_bottom = matrix(independent_likelihoods_at_branch_section_bottom %*% relative_probs_of_each_state_at_the_tip_of_this_branch, nrow=1)		
					
							# 2015 version
							tmp_conditional_likelihoods_at_branch_section_bottom = matrix(independent_likelihoods_at_branch_section_bottom %*% relative_probs_of_each_state_at_the_tip_of_this_branch[states_allowed_TF], nrow=1)
							} else {
							tmp_conditional_likelihoods_at_branch_section_bottom = matrix(independent_likelihoods_at_branch_section_bottom %*% relative_probs_of_each_state_at_the_tip_of_this_branch[wTrait_states_allowed_TF], nrow=1)
							}
						} else {
						# SPARSE MATRIX EXPONENTIATION DOWN ONE BRANCH
						# CHECK THAT IT'S IN COO FORMAT
						if (class(Qmat_tmp) != "data.frame")
							{
							txt = paste0("ERROR: calc_loglike_sp_stratified is attempting to use a sparse COO-formated Q matrix, to calculated the likelihoods down one branch segment, but you provided a Qmat not in data.frame form")
							cat("\n\n")
							cat(txt)
							cat("\n\n")
							print("class(Qmat_tmp)")
							print(class(Qmat_tmp) )
							print("Qmat_tmp")
							print(Qmat_tmp)
							stop(txt)
							}
			
						if ( (ncol(Qmat_tmp) != 3) )
							{
							txt = paste0("ERROR: calc_loglike_sp_stratified is attempting to use a sparse COO-formated Q matrix, to calculated the likelihoods down one branch segment, but you provided a Qmat that does't have 3 columns")
							cat("\n\n")
							cat(txt)
							cat("\n\n")
							print("class(Qmat_tmp)")
							print(class(Qmat_tmp) )
							print("Qmat_tmp")
							print(Qmat_tmp)
							stop(txt)
							}
						coo_n = numstates_geogtrait
						anorm = 1
						
						#print(relative_probs_of_each_state_at_the_tip_of_this_branch[states_allowed_TF])

						if (traitTF == FALSE)
							{
							try_result_segment = try (
							kexpmv::expokit_dgexpv(mat=tmpQmat_in_kexpmv_crs_fmt, t=treepiece, vector=relative_probs_of_each_state_at_the_tip_of_this_branch[states_allowed_TF], transpose_needed=TRUE, transform_to_crs=FALSE, crs_n=numstates_geogtrait, anorm=NULL, mxstep=10000, tol=as.numeric(1e-10))
							)
							} else {
							try_result_segment = try (
							kexpmv::expokit_dgexpv(mat=tmpQmat_in_kexpmv_crs_fmt, t=treepiece, vector=relative_probs_of_each_state_at_the_tip_of_this_branch[wTrait_states_allowed_TF], transpose_needed=TRUE, transform_to_crs=FALSE, crs_n=numstates_geogtrait, anorm=NULL, mxstep=10000, tol=as.numeric(1e-10))
							)
							}
							
						if (printlevel >=1)
							{	
							txt = "S"
							cat(txt)
							}

						# Error check
						if (class(try_result_segment) == "try-error")
							{
							cat("\n\ntry-error on kexpmv::expokit_dgexpv():\n\n")
							cat("i=", i, "\n")
							cat("t=treepiece==", treepiece, "\n")
							print(tmpQmat_in_kexpmv_crs_fmt)
							print(phy2)
							print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS)
							print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,])
							print(coo_n)
							print(anorm)
							print("BioGeoBEARS_model_object")
							print(BioGeoBEARS_model_object)
				
							stoptxt = "\n\nStopping on error in sparse exponentiation downpass (treepiece, aka a branch segment): NaNs produced in likelihood calculation. This may mean your transition matrix disallows necessary transitions.  E.g., if your ranges are 'A' and 'B', and your model is DEC, then allowing range 'AB' as a possible state is required, so that you can get from 'A' to 'B' via 'AB' as the intermediate. Alternatively, NaNs can be produced sometimes if your Maximum Likelihood (ML) search proposes weird parameter values (such as a negative rate or weight) or a parameter so small that required transitions have a probability that machine precision rounds to zero or negative.  Sometimes this seems to occur because optim, optimx, etc. propose parameters slightly outside the user-specified upper and lower (min/max) boundaries for some reason. One solution is often to narrow the min/max limits. \n\nAnother solution: To have this error report an extremely low log-likelihood,, set BioGeoBEARS_run_object$on_NaN_error to something like -1e50.\n\n"
							
							if (is.null(on_NaN_error))
								{
								stop(stoptxt)
								}
				
							print("print(on_NaN_error):")
							print(on_NaN_error)
							if ( (is.numeric(on_NaN_error)) && (return_what == "loglike") )
								{
								warning(paste0("\n\nWarning  on error in sparse exponentiation downpass (treepiece, aka a branch segment): NaNs produced in likelihood calculation. This may mean your transition matrix disallows necessary transitions.  E.g., if your ranges are 'A' and 'B', and your model is DEC, then allowing range 'AB' as a possible state is required, so that you can get from 'A' to 'B' via 'AB' as the intermediate. Alternatively, NaNs can be produced sometimes if your Maximum Likelihood (ML) search proposes weird parameter values (such as a negative rate or weight) or a parameter so small that required transitions have a probability that machine precision rounds to zero or negative.  Sometimes this seems to occur because optim, optimx, etc. propose parameters slightly outside the user-specified upper and lower (min/max) boundaries for some reason. One solution is often to narrow the min/max limits. \n\nYou are using another solution: Normally, this would be a stop error, but you specified that BioGeoBEARS_run_object$on_NaN_error=", on_NaN_error, "\n\n"))
								return(on_NaN_error)
								} else {
								stop(stoptxt)
								}
							} # END if (any(is.nan(condlikes_Left)))
						
						# If all checks are survived, get the downpass probabilities
						tmp_conditional_likelihoods_at_branch_section_bottom = c(try_result_segment$output_probs)
						tmp_conditional_likelihoods_at_branch_section_bottom[tmp_conditional_likelihoods_at_branch_section_bottom<0] = tmp_conditional_likelihoods_at_branch_section_bottom
						} # END if (force_sparse == FALSE)
					
					
					if (traitTF == FALSE)
						{
						# save to full matrix
						conditional_likelihoods_at_branch_section_bottom = matrix(0, nrow=1, ncol=length(relative_probs_of_each_state_at_the_tip_of_this_branch))
						conditional_likelihoods_at_branch_section_bottom[states_allowed_TF] = tmp_conditional_likelihoods_at_branch_section_bottom
						} else {
						# save to full matrix
						conditional_likelihoods_at_branch_section_bottom = matrix(0, nrow=1, ncol=length(relative_probs_of_each_state_at_the_tip_of_this_branch))
						conditional_likelihoods_at_branch_section_bottom[wTrait_states_allowed_TF] = tmp_conditional_likelihoods_at_branch_section_bottom
						}
					
					
					# Zero out impossible states according to
					# areas_allowed/areas_adjacency
					# keep from 2014 just to double-check
					conditional_likelihoods_at_branch_section_bottom[states_allowed_TF==FALSE] = 0
					
					#print("conditional_likelihoods_at_branch_section_bottom #1a:")
					#print(conditional_likelihoods_at_branch_section_bottom)
					} else {
					# Copying the tip likelihoods down
					conditional_likelihoods_at_branch_section_bottom = matrix(relative_probs_of_each_state_at_the_tip_of_this_branch, nrow=1)
					#print("conditional_likelihoods_at_branch_section_bottom #1b:")
					#print(conditional_likelihoods_at_branch_section_bottom)
					} # END if (do_exponentiation == TRUE)
				
				
				# Test forward exponentiation instead...NO
				# independent_likelihoods_at_branch_section_bottom = expokit_dgpadm_Qmat2(times=treepiece, Qmat=Qmat_tmp,  transpose_needed=FALSE)
# 				#independent_likelihoods_at_branch_section_bottom = expokit_dgpadm_Qmat(Qmat=Qmat_tmp,  t=treepiece, transpose_needed=FALSE)
# 				
# 				
# 				conditional_likelihoods_at_branch_section_bottom = matrix(independent_likelihoods_at_branch_section_bottom %*% relative_probs_of_each_state_at_the_tip_of_this_branch, nrow=1)
#				if (include_null_range == TRUE)
# 					conditional_likelihoods_at_branch_section_bottom[1] = 0
# 				
				
				
				
				# Also, store the conditional likelihoods for all nodes in this subtree
				chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]] = conditional_likelihoods_at_branch_section_bottom
	
	
				# (THIS IS CRUCIAL TO GETTING STRATIFICATION TO WORK -- YOU NEED THE is_fossil==TRUE ADDED!!)
				# (these don't seem essential, whether divided or not, in stratified analysis; what matters is what 
				#  goes into condlikes_table)
				# Relative probabilities -- just the new tip
				chainsaw_result$relative_probs_of_each_state_at_bottom_of_root_branch[[jj]] = conditional_likelihoods_at_branch_section_bottom / sum(conditional_likelihoods_at_branch_section_bottom)
	
				# Relative probabilities -- all nodes plus branch bottom (just branch bottom, here)
				chainsaw_result$relative_probabilities_for_nodes_plus_bottom_in_this_section[[jj]] = conditional_likelihoods_at_branch_section_bottom / sum(conditional_likelihoods_at_branch_section_bottom)
				

				#print("conditional_likelihoods_at_branch_section_bottom #2:")
				#print(conditional_likelihoods_at_branch_section_bottom)

				#print("sum(conditional_likelihoods_at_branch_section_bottom)")
				#print(sum(conditional_likelihoods_at_branch_section_bottom))

				
				#print("conditional_likelihoods_at_branch_section_bottom / sum(conditional_likelihoods_at_branch_section_bottom)")
				#print(conditional_likelihoods_at_branch_section_bottom / sum(conditional_likelihoods_at_branch_section_bottom))
				
				
				
				# If you are storing ALL of the conditional likelihoods that were calculated
				if ((return_condlikes_table == TRUE) || (calc_TTL_loglike_from_condlikes_table == TRUE))
					{
					# Find the row in the big conditional likelihoods table
					TF1 = inputs$master_table$stratum == i
					TF2 = inputs$master_table$piecenum == jj
					
					# 2017-04-06 fix: do BOTH subbranch and orig tip, so that we STORE the 
					# downpass probabilities at branch bottoms, for the original tips
					# This should help with stochastic mapping...
					TF3 = inputs$master_table$piececlass == "subbranch"
					TF4 = inputs$master_table$piececlass == "orig_tip"
					TF5 = (TF3 + TF4) == 1
					TF = (TF1 + TF2 + TF5) == 3
				
					rownum = (1:nrow(condlikes_table))[TF]
					condlikes_table[rownum, ] = conditional_likelihoods_at_branch_section_bottom
					
					# Also store the subbranch downpass relative probabilities at the bottom of each branch
					if (calc_ancprobs == TRUE)
						{
						# Also store the subbranch downpass relative probabalities at the bottom of each branch
						relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE[rownum, ] = conditional_likelihoods_at_branch_section_bottom / sum(conditional_likelihoods_at_branch_section_bottom)
						} # END if (calc_ancprobs == TRUE)
					} # END if ((return_condlikes_table == TRUE) || (calc_TTL_loglike_from_condlikes_table == TRUE))
				
				
				############################################
				# END if (is.numeric(treepiece))
				# It's just a branch section
				############################################
				} else {
				############################################
				# DOWNPASS -- on a subtree
				############################################
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
						this_row_of_master_table_is_being_used = TF
					
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
							# 2016-02-29 -- this adjustment now happens during section_the_tree()
							#
							
							#tmp_subtree$edge.length[tmp2_edgenum] = tmp_subtree$edge.length[tmp2_edgenum] - amount_to_shorten_by
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
				
				if (traitTF == FALSE)
					{
					# 2014 version
					#subtree_tip_relative_probs_of_each_state = current_tip_relative_probs_of_each_state[tips_for_subtree_TF,states_to_use_TF]
					# 2015 version
					subtree_tip_relative_probs_of_each_state = current_tip_relative_probs_of_each_state[tips_for_subtree_TF,][,states_allowed_TF]
					} else {
					subtree_tip_relative_probs_of_each_state = current_tip_relative_probs_of_each_state[tips_for_subtree_TF,][,wTrait_states_allowed_TF]
					}
				
				# 2016-05-28_bug_fix
				# Fix this error, e.g. when DEC* model + areas_allowed means that
				# ranges_list = NULL, Kauai is just
				# ranges_list = Kauai
				# This means that:
				#   subtree_tip_relative_probs_of_each_state 
				# and thus
				#   tip_condlikes_of_data_on_each_state
				# ...are just a list of numbers, not a matrix, thus 
				# rowSums fails in calc_loglike_sp() in that time-stratum.
				# 
				#
				# This was the error message:
				# 
				# Error in rowSums(tip_condlikes_of_data_on_each_state) : 
				#  'x' must be an array of at least two dimensions
				# Calls: bears_optim_run ... calc_loglike_sp_stratified -> calc_loglike_sp -> rowSums
				#
				
				# If there is only 1 geographic state...
				if (sum(states_allowed_TF) == 1)
					{
					if (traitTF == FALSE)
						{
						subtree_tip_relative_probs_of_each_state = matrix(data=subtree_tip_relative_probs_of_each_state, ncol=1)
						} else {
						# If there's a trait, there are at least 2 geographic states
						subtree_tip_relative_probs_of_each_state = matrix(data=subtree_tip_relative_probs_of_each_state, ncol=sum(wTrait_states_allowed_TF))
						} # END if (traitTF == FALSE)
					} # END if (sum(states_allowed_TF) == 1)

				
				
				
				# DOWNPASS: check if this subtree contains fixed internal node(s) on the master tree
				
				# Match the master fixnodes to the fixnodes in *JUST* this subtree
				# We will then pass these fixnodes to the subtree loglike calculation
				# First, we need to get the master node number, iff it's internal
				# 
				tmp_fixnode = NULL		# Default
				tmp_fixlikes = NULL		# Default
				if ((!is.null(fixnode)) && (length(fixnode) > 0))
					{
					# Check for multiple fixnodes
					if (length(fixnode) > 1)
						{
						# If there are multiple fixnodes, 
						# Get the matching nodes in this subtree
						TF1 = inputs$master_table$stratum == i
						TF2 = inputs$master_table$piecenum == jj
						TF3 = inputs$master_table$piececlass == "subtree"
						TF = ((TF1 + TF2 + TF3) == 3)
						tmprows = inputs$master_table[TF,]
						
						# Get the fixnodes found in this subtree
						fixnodes_in_subtree_TF = fixnode %in% tmprows$node
					
						# *IF* the subtree contains fixnodes, do this stuff
						# otherwise, set to NULL
						if (sum(fixnodes_in_subtree_TF) > 0)
							{
							#master_nodes_in_fixnode_TF = inputs$master_table$node %in% fixnode
							#master_nodes_in_fixnode
						
							#TF = (anc == fixnode)	# old
							# we do not use temporary_fixnode, since we need the fixnodes in the subtree numbering (tmprow$SUBnode)
							temporary_fixnodes = fixnode[fixnodes_in_subtree_TF]
							
							if (traitTF == FALSE)
								{
								# But we will use these
								# 2016-03-15_old
								#temporary_fixlikes = fixlikes[fixnodes_in_subtree_TF,]
								# 2016-03-15_new by Torsten
								temporary_fixlikes = fixlikes[fixnodes_in_subtree_TF,states_allowed_TF]
								} else {
								temporary_fixlikes = fixlikes[fixnodes_in_subtree_TF,wTrait_states_allowed_TF]
								}
							
							# The subtree nodenums corresponding to the subset temporary_fixnodes
							# NOTE! THIS SUBSET THING WILL ONLY WORK IF THE NODES ARE SORTED IN ORDER FROM THE START
							subtree_rows_in_fixnodes_TF = tmprows$node %in% fixnode 
							subtree_fixnode_master_nodenums = tmprows$node[subtree_rows_in_fixnodes_TF]
							subtree_fixnode_nums = tmprows$SUBnode[subtree_rows_in_fixnodes_TF]
							
							# We have to order these subtree fixnodes, and order the subtree fixlikes the same way
							order_subtree_fixnode_nums = order(subtree_fixnode_nums)
							subtree_fixnode_nums = subtree_fixnode_nums[order_subtree_fixnode_nums]
							
							# Only reorder if there are 2 or more rows, i.e. if it's a matrix not a vector
							if (length(order_subtree_fixnode_nums) > 1)
								{
								temporary_fixlikes = temporary_fixlikes[order_subtree_fixnode_nums, ]
								}
							
							} else {
							# If *NO* fixnodes in subtree:
							temporary_fixnodes = NULL
							subtree_fixnode_master_nodenums = NULL
							subtree_fixnode_nums = NULL
							temporary_fixlikes = NULL
							} # end if (sum(fixnodes_in_subtree_TF) > 0)


						# Check if we're in the right stratum / piece / piececlass
						# (have account for possible multiple rows)
						TF1 = unique(tmprows$stratum) == i
						TF2 = unique(tmprows$piecenum) == jj
						TF3 = unique(tmprows$piececlass) == "subtree"
						TF = ((TF1 + TF2 + TF3) == 3)
						
						if (TF == TRUE)
							{
							#txt = paste("Master tree node ", fixnode, " matched to i=", i, "; jj=", jj, "; piececlass=", piececlass, "; subtree subnode=", tmprow$SUBnode, sep="")
							#print(txt)
							#print(fixlikes)
				
							# Determine the number of the subnode in the subtree
							if (length(subtree_fixnode_nums) == 0)
								{
								subtree_fixnode_nums = NULL
								temporary_fixlikes = NULL
								}
						
							tmp_fixnode = subtree_fixnode_nums
							tmp_fixlikes = temporary_fixlikes						
							} else {
							tmp_fixnode = NULL
							tmp_fixlikes = NULL
							} # end if (TF == TRUE)


						} else {
						# Only 1 fixnode
						temporary_fixnode = fixnode
						temporary_fixlikes = c(fixlikes)

						# e.g.
						# fixnode=20
						TF1 = inputs$master_table$node == temporary_fixnode
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
							
# 							cat("\n\n")
# 							print(fixnode)
# 							print(temporary_fixnode)
# 							print(tmprow$SUBnode)
# 							print(temporary_fixlikes)
							

							# Determine the number of the subnode in the subtree
							tmp_fixnode = tmprow$SUBnode
							# 2016-03-15_old
							#tmp_fixlikes = temporary_fixlikes
							# 2016-03-15_new by Torsten
							tmp_fixlikes = temporary_fixlikes[states_allowed_TF]
							# end if (TF == TRUE)
							} else {
							tmp_fixnode = NULL
							tmp_fixlikes = NULL
							}

						} # end if (length(fixnode) > 1)
					} # end if (!is.null(fixnode))
				

								
				# Calculate the likelihoods for this subtree
				#prt(tmp_subtree)
				#print("subtree_tip_relative_probs_of_each_state")
				#print(subtree_tip_relative_probs_of_each_state)
				#print("min_branchlength")
				#print(min_branchlength)
				
# 				print("subtree_tip_relative_probs_of_each_state:")
# 				print(subtree_tip_relative_probs_of_each_state)
# 				print("spPmat_inputs:")
# 				print(spPmat_inputs)
# 				print("Qmat_tmp:")
# 				print(Qmat_tmp)
				
				calc_loglike_sp_results = calc_loglike_sp(
					tip_condlikes_of_data_on_each_state=subtree_tip_relative_probs_of_each_state, 
					phy=tmp_subtree, 
					Qmat=Qmat_tmp, 
					spPmat=NULL,
					min_branchlength=min_branchlength,
					return_what="all",
					probs_of_states_at_root=NULL,
					rootedge=TRUE,
					sparse=force_sparse, 
					printlevel=printlevel,
					use_cpp=TRUE,
					input_is_COO=force_sparse,
					spPmat_inputs=spPmat_inputs,
					cppSpMethod=cppSpMethod,
					cluster_already_open=cluster_already_open,
					calc_ancprobs=calc_ancprobs,	 # If TRUE, get e.g. relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS
					include_null_range=include_null_range,
					fixnode=tmp_fixnode,
					fixlikes=tmp_fixlikes,
					stratified=TRUE,		# This makes calc_loglike_sp skip UPPASS probs, which are irrelevant inside stratified analyses
					# 2014:states_allowed_TF=states_allowed_TF
					states_allowed_TF=rep(TRUE, times=ncol(subtree_tip_relative_probs_of_each_state)),
					m=m,
					jts_matrix=jts_matrix,
					BioGeoBEARS_model_object=BioGeoBEARS_model_object,
					on_NaN_error=BioGeoBEARS_run_object$on_NaN_error
					)

				# Slot these likelihoods into a bigger object, if needed due to
				# lists_of_states_lists_0based
				if (!is.null(inputs$lists_of_states_lists_0based))
					{
					#print("Here!!")
					#print(calc_loglike_sp_results)
					#print("Here!!")
					names_of_calc_loglike_sp_results_objects = names(calc_loglike_sp_results)
					for (name_i in 1:length(calc_loglike_sp_results))
						{
						# If it's a matrix, slot it inside a new matrix
						oldmat = calc_loglike_sp_results[[name_i]]
						
						TF1 = names_of_calc_loglike_sp_results_objects[name_i] == "relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS"
						TF2 = names_of_calc_loglike_sp_results_objects[name_i] == "condlikes_of_each_state"
						TF3 = names_of_calc_loglike_sp_results_objects[name_i] == "relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS"
						TF4 = names_of_calc_loglike_sp_results_objects[name_i] == "relative_probs_of_each_state_at_bottom_of_root_branch"
						
						if (TF1 || TF2 || TF3)
							{
							# Old
							if (traitTF == FALSE)
								{
								newmat = matrix(0, nrow=nrow(oldmat), ncol=length(states_allowed_TF))
								newmat[,states_allowed_TF] = oldmat
								}
							
							# New, with traits possible
							if (traitTF == TRUE)
								{
								full_matrix_ncols = length(states_allowed_TF) * num_trait_states
								newmat = matrix(0, nrow=nrow(oldmat), ncol=full_matrix_ncols)
								wTrait_states_allowed_TF = c(rep(states_allowed_TF, times=num_trait_states))
								newmat[,wTrait_states_allowed_TF] = oldmat
								}
							
							calc_loglike_sp_results[[name_i]] = newmat
							#print("oldmat")
							#print(oldmat)
							#print("newmat")
							#print(newmat)
							} # END if (is.matrix(oldmat))

						if (TF4)
							{
							# Old
							if (traitTF == FALSE)
								{
								newmat = matrix(0, nrow=1, ncol=length(states_allowed_TF))
								newmat[,states_allowed_TF] = oldmat
								}

							# New, with traits possible
							if (traitTF == TRUE)
								{
								full_matrix_ncols = length(states_allowed_TF) * num_trait_states
								newmat = matrix(0, nrow=1, ncol=full_matrix_ncols)
								wTrait_states_allowed_TF = c(rep(states_allowed_TF, times=num_trait_states))
								newmat[,wTrait_states_allowed_TF] = oldmat
								}

							calc_loglike_sp_results[[name_i]] = newmat
							#print("oldmat")
							#print(oldmat)
							#print("newmat")
							#print(newmat)
							} # END if (is.matrix(oldmat))

						} # END for (name_i in 1:length(calc_loglike_sp_results))
					} # END if (!is.null(inputs$lists_of_states_lists_0based))

				
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
							# Store the state probabilities at the branch bottoms below nodes
							# NOTE: calc_loglikes_sp() returns NA for the bottom of the root branch, and 
							# stores that instead in 
							if (rownum <= nrow(calc_loglike_sp_results$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS))
								{
								# Check if you are the subtree root or not
								if (inputs$master_table$SUBnode.type[condlikes_table_rownum] != "root")
									{
									# Get relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS
									# For subtree tip and internal nodes
									tmp = calc_loglike_sp_results$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[rownum,]
																		relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE[condlikes_table_rownum,] = tmp
									} else {
									# For subtree root node
									tmp = calc_loglike_sp_results$relative_probs_of_each_state_at_bottom_of_root_branch
																		relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE[condlikes_table_rownum,] = tmp
									}
								} # END check of subtree internal/tip vs. subtree root
							# Store the state probabilities at the branch bottom below the root node of the subtree	
							#if (rownum == nrow(calc_loglike_sp_results$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS))
							#	{
								# Get relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS
								#tmp = NA
								# Save the relative probabilities of each state at the BOTTOM of the branch
								# BELOW the subtree root
							# 	tmp = calc_loglike_sp_results$relative_probs_of_each_state_at_bottom_of_root_branch
# 								relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE[condlikes_table_rownum,] = tmp
# 								
# 								cat("\n\n")
# 								print(i)
# 								print(jj)
# 								print(rownum)
# 								print(calc_loglike_sp_results$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS)
# 								print(relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE[condlikes_table_rownum,])
# 								print(condlikes_table_rownum)
# 								}

							} # END if (calc_ancprobs == TRUE)
						} # END for (rownum in 1:nrow(calc_loglike_sp_results$condlikes_of_each_state))
					} # END if ((return_condlikes_table == TRUE) || (calc_TTL_loglike_from_condlikes_table == TRUE))

				
				
				chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]] = matrix(data=calc_loglike_sp_results$condlikes_of_each_state[-tmp_tipnums, ], ncol=ncol(calc_loglike_sp_results$condlikes_of_each_state))
				
		
				# Matrix of tip likelihoods to delete so you don't repeat using them in the total
				# loglike
				#tiplikes_to_delete[[jj]] = calc_loglike_sp_results$condlikes_of_each_state[tmp_tipnums, ]
				
				# Relative probabilities -- all nodes plus branch bottom (just branch bottom, here)
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
				#if (i != num_timeperiods)
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

		} # END for (i in 1:num_timeperiods)
		# END loop through i strata

		##################################################################	
		##################################################################
		##################################################################
		# ENDING DOWNPASS
		##################################################################	
		##################################################################
		##################################################################


	# Remove rows that have not been filled (till zero)
	all_condlikes_of_each_state_zero_TF = all_condlikes_of_each_state == 0
	all_condlikes_of_each_state_nonzero_TF = all_condlikes_of_each_state_zero_TF == FALSE
	rows_that_are_NOT_numeric_zeros_TF = rowSums(all_condlikes_of_each_state_nonzero_TF) >= 1
	#rowSums(all_condlikes_of_each_state) != 0
	final_all_condlikes_of_each_state = all_condlikes_of_each_state[rows_that_are_NOT_numeric_zeros_TF,]
	
	#rowSums(all_relative_probs_of_each_state) != 0
	all_relative_probs_of_each_state = all_relative_probs_of_each_state[rows_that_are_NOT_numeric_zeros_TF,]
	
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
		"e.g. the 0 values to e.g. 0.0000001.  Good luck!", sep="")
		
		cat(stoptxt1)
		
		stop(stoptxt1)
		}
	
	
	if (calc_TTL_loglike_from_condlikes_table == TRUE)
		{
		#print ("HEY!")
		
		# Standard LAGRANGE result (exactly)
		TF2 = inputs$master_table$SUBnode.type == "internal"
		TF3 = inputs$master_table$SUBnode.type == "orig_tip"	# These are the original tips likelihoods; doesn't matter for unambiguous tips, but
																# DOES matter if there is a detection model.
		TF4 = inputs$master_table$SUBnode.type == "root"
		TF234 = (TF2 + TF3 + TF4) == 1
		sum(TF234)
		
		#TF5 = inputs$master_table$piececlass == "subbranch"
		#TF234 = (TF2 + TF3 + TF4 + TF5) == 1

		
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

			tmptable = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE[TF,][node_order_original,]
			relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS = tmptable / rowSums(tmptable)
			
			# This leaves out the master tree root row, so add that in
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
			anc_node_original_tree = inputs$master_table$node[anc_row_of_master_table_TF]
			anc_node_original_tree
			
			# Just always use the root node, not anything below it!
			starting_probs = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[anc_node_original_tree, ]
			} # END if (calc_ancprobs == TRUE)
		} # END if (calc_TTL_loglike_from_condlikes_table == TRUE)




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
		relative_probs_of_each_state_at_branch_top_AT_node_UPPASS = matrix(data=0, nrow=numrows_for_UPPASS, ncol=ncol(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS))
		relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS = matrix(data=0, nrow=numrows_for_UPPASS, ncol=ncol(relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS))
		
		
		# Vist edges in reverse order from the downpass
		# This would only work on un-stratified trees
		#edges_to_visit_uppass = seq(from=(num_internal_nodes*2), by=-2, length.out=num_internal_nodes)

		
		# Get the starting probabilities at the root
		# THIS ASSUMES THE STARTING PROBS ARE AT THE ROOT NODE, NOT SOME STUPID BRANCH BELOW THE ROOT NODE
		starting_probs
		
		# Put this starting prob into the root node
		relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[anc_node_original_tree,] = 1/length(starting_probs)
		
		#print("starting_probs")
		#print(starting_probs)
		#print(relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[(anc_node_original_tree-3):(anc_node_original_tree+3),])
		#print(dim(relative_probs_of_each_state_at_branch_top_AT_node_UPPASS))
		
		
		# Go through strata in REVERSE order
		#print("num_timeperiods")
		#print(num_timeperiods)
		#print("inputs$list_of_areas_adjacency_mats")
		#print(inputs$list_of_areas_adjacency_mats)
		for (i in num_timeperiods:1)
			{
	
			#######################################################
			# Cut down the number of areas, by what is allowed
			# (it would be more efficient to do this once during setup, 
			#  but probably no biggie)
			#######################################################
			# states_to_use_TF: states to use in Qmat, speciation models, etc.
			# states_allowed_TF: use this to zero out impossible ancestral 
			# states according to areas_allowed matrix/areas_adjacency matrix
			# 
			
			# Should we modify the list of allowed states?
			# default: no areas_allowed or areas_adjacency constraints
			user_specified_constraints_on_states_list_TF = FALSE
			states_allowed_TF1 = rep(TRUE, length(all_states_list))
			states_allowed_TF2 = rep(TRUE, length(all_states_list))
			if ( (is.null(inputs$list_of_areas_allowed_mats) == FALSE))
				{
				user_specified_constraints_on_states_list_TF = TRUE
				}
			if ( (is.null(inputs$list_of_areas_adjacency_mats) == FALSE))
				{
				user_specified_constraints_on_states_list_TF = TRUE
				}


		
			if (user_specified_constraints_on_states_list_TF == TRUE)
				{
				# Areas allowed
				if ( (is.null(inputs$list_of_areas_allowed_mats) == FALSE))
					{
					areas_allowed_mat = inputs$list_of_areas_allowed_mats[[i]]
					
					cat("\ni=", i, "\n", sep="")
					cat("areas_allowed_mat: ", sep="")
					print(areas_allowed_mat)
					
					states_allowed_TF1 = sapply(X=all_states_list, FUN=check_if_state_is_allowed, areas_allowed_mat)
					#states_to_use_TF = all_states_list %in% tmp_states_list
			
					if (include_null_range == TRUE)
						{
						states_allowed_TF1[1] = TRUE
						}
					# NO; use all areas for this
					# states_to_use_TF = states_allowed_TF
					} # END if ( (is.null(inputs$list_of_areas_allowed_mats) == FALSE))
			
				# Areas adjacency
				if ( (is.null(inputs$list_of_areas_adjacency_mats) == FALSE))
					{
					areas_adjacency_mat = inputs$list_of_areas_adjacency_mats[[i]]
	
					states_allowed_TF2 = sapply(X=all_states_list, FUN=check_if_state_is_allowed_by_adjacency, areas_adjacency_mat)
					#states_to_use_TF = all_states_list %in% tmp_states_list
			
					if (include_null_range == TRUE)
						{
						states_allowed_TF2[1] = TRUE
						}
					# NO; use all areas for this
					# states_to_use_TF = states_allowed_TF
					} # END if ( (is.null(inputs$list_of_areas_adjacency_mats) == FALSE))
				# Combine the two (areas_allowed and areas_adjacency)
				states_allowed_TF = ((states_allowed_TF1 + states_allowed_TF2) == 2)
				} else {
				# Otherwise, 
				# make no change
				pass = 1
				#states_list = states_list
				states_allowed_TF = rep(TRUE, length(all_states_list))
				} # END if (user_specified_constraints_on_states_list_TF == TRUE)
			# Use this for regular calculations (Qmat, speciation models, etc.)
			states_to_use_TF = rep(TRUE, length(all_states_list))
			
			#print("states_allowed_TF")
			#print(states_allowed_TF)
				
			#####################################################
			# Make the dedf matrix for this time period
			#####################################################
			# If there is a distance matrix, use the first one 
			# (non-stratified analysis, here)
	
			# If there is a distance matrix, take the ith one... 
			# (stratified analysis, here)
			if ( (is.null(inputs$list_of_distances_mats) == FALSE))
				{
				distances_mat = inputs$list_of_distances_mats[[i]]
				} else {
				# Default is all areas effectively equidistant
				distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))
				}
		
			# Get the exponent on distance, apply to distances matrix
			dispersal_multipliers_matrix = distances_mat ^ x

			# Environmental distances
			if ( (is.null(inputs$list_of_envdistances_mats) == FALSE))
				{
				envdistances_mat = inputs$list_of_envdistances_mats[[1]]
				} else {
				# Default is all areas effectively equidistant
				envdistances_mat = matrix(1, nrow=length(areas), ncol=length(areas))
				}

			# Get the exponent on environmental distance, apply to distances matrix
			n = BioGeoBEARS_model_object@params_table["n","est"]
			dispersal_multipliers_matrix = dispersal_multipliers_matrix * envdistances_mat^n

		
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
			
			
			
			
			
			# 2018 version
			if (traitTF == FALSE)
				{
				Qmat_tmp = rcpp_states_list_to_DEmat(areas_list=allareas_list, states_list=all_states_list[states_allowed_TF], dmat=dmat_times_d, elist=elist, amat=amat, include_null_range=include_null_range, normalize_TF=TRUE, makeCOO_TF=force_sparse)

				#print(dim(Qmat_tmp))
	# 			} else {
	# 			# If Qmat is pre-specified
	# 			Qmat_tmp = Qmat
	# 			}
				} # END if (traitTF == FALSE)

			# Analysis with a trait modifying dispersal rate
			if (traitTF == TRUE)
				{
				num_geog_states = length(all_states_list[states_allowed_TF])
				numstates_geogtrait = num_trait_states * num_geog_states
			
				# UPPASS definition of states_allowed_TF with traits
				wTrait_states_allowed_TF =c(rep(states_allowed_TF, times=num_trait_states))
			
# 				print("num_geog_states")
# 				print(num_geog_states)
# 			
# 				print("num_trait_states")
# 				print(num_trait_states)
# 			
# 				print("numstates_geogtrait")
# 				print(numstates_geogtrait)
			
				if (ncol(tip_condlikes_of_data_on_each_state[,wTrait_states_allowed_TF]) != numstates_geogtrait)
					{
					txt = paste0("STOP ERROR in calc_loglike_sp_stratified(): ncol(tip_condlikes_of_data_on_each_state)=", ncol(tip_condlikes_of_data_on_each_state), ", and numstates_geogtrait=", numstates_geogtrait, ". They must be equal to proceed.")
					cat("\n\n")
					cat(txt)
					cat("\n\n")
					stop(txt)
					} # END if (ncol(tip_condlikes_of_data_on_each_state) != numstates_geogtrait)

				# Get the modified Qmatrix (traits + geog)			
				tmpres = modify_Qmat_with_trait(Qmat=NULL, BioGeoBEARS_run_object, numstates_geogtrait=numstates_geogtrait, areas_list=allareas_list, states_list=all_states_list[states_allowed_TF], dispersal_multipliers_matrix=dispersal_multipliers_matrix, elist=elist, force_sparse=force_sparse)
				Qmat_tmp = tmpres$Qmat
				m = tmpres$m
			
			
				# If the trait can change during jump events
				if (is.null(BioGeoBEARS_run_object$jts_txt_matrix) == FALSE)
					{
					jts_txt_matrix = BioGeoBEARS_run_object$jts_txt_matrix
					jts_matrix = matrix(data=0, nrow=nrow(jts_txt_matrix), ncol=ncol(jts_txt_matrix))
					TF_matrix = matrix(data=TRUE, nrow=nrow(jts_txt_matrix), ncol=ncol(jts_txt_matrix))
					diag(TF_matrix) = FALSE
					jts_txt_params = c(jts_txt_matrix[TF_matrix])
					jts_txt_params
			
					# Populate the numeric jts_matrix
					for (jts_i in 1:nrow(jts_txt_matrix))
						{
						diag_val = 1
						for (jts_j in 1:ncol(jts_txt_matrix))
							{
							if (jts_i == jts_j)
								{
								next()
								}
							jts_txt = jts_txt_matrix[jts_i,jts_j]
							newval = as.numeric(BioGeoBEARS_model_object@params_table[jts_txt, "est"])
							jts_matrix[jts_i,jts_j] = newval
							diag_val = 1-newval
							}
						# Populate the diagonal
						jts_matrix[jts_i,jts_i] = diag_val
						} # END for (jts_i in 1:nrow(jts_txt_matrix))
					} # END if (is.null(BioGeoBEARS_run_object$jts_txt_matrix) == FALSE)
				} else {
				num_geog_states = length(all_states_list[states_allowed_TF])
				numstates_geogtrait = num_geog_states
				} # END if (traitTF == TRUE)


		if (force_sparse == TRUE)
			{
			tmpQmat_in_REXPOKIT_coo_fmt = Qmat_tmp

			# Make a CRS-formatted matrix, for kexpmv
			# DO THE TRANSPOSE HERE, trait+geog matrices assembled transposed
			tmpQmat_in_kexpmv_crs_fmt = coo2crs(
				ia=tmpQmat_in_REXPOKIT_coo_fmt[,"ia"], 
				ja=tmpQmat_in_REXPOKIT_coo_fmt[,"ja"], 
				a =tmpQmat_in_REXPOKIT_coo_fmt[,"a"],
				n=numstates_geogtrait, transpose_needed=FALSE)
			} # END if (force_sparse == TRUE)

					
# 			if (sparse == TRUE)
# 				{
# 				# Sparse matrix exponentiation
# 				original_Qmat = Qmat_tmp
# 				
# 				# number of states in the original matrix
# 				coo_n = ncol(Qmat_tmp)
# 				anorm = as.numeric(norm(original_Qmat, type="O"))
# 				matvec = original_Qmat
# 				
# 				# *DO* TRANSPOSE; we want to go FORWARDS in time, NOT BACKWARDS!
# 				tmatvec = base::t(matvec)
# 				tmatvec = matvec
# 				tmpQmat_in_REXPOKIT_coo_fmt = mat2coo(tmatvec)
# 				}
			
			
			# Now. IF you have a subtree structure, you need to run this with a cladogenesis matrix, 
			# through calc_loglike_sp(), like normal.
			
			# If there's just one tree, store it in the object
			if (is.null(inputs$timeperiods) || length(inputs$timeperiods) == 1)
				{
				#tr = read.tree(inputs$trfn)
				tr = check_trfn(trfn=inputs$trfn)

				
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
			spPmat_inputs = get_spPmat_inputs_from_BGB(BioGeoBEARS_run_object=BioGeoBEARS_run_object, states_list=all_states_list[states_allowed_TF], dispersal_multipliers_matrix=dispersal_multipliers_matrix)
		
			dmat = dispersal_multipliers_matrix
		
			maxent01s_param = spPmat_inputs$maxent01s_param
			maxent01v_param = spPmat_inputs$maxent01v_param
			maxent01j_param = spPmat_inputs$maxent01j_param
			maxent01y_param = spPmat_inputs$maxent01y_param
			
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



			# -1 for null range
			if (include_null_range == TRUE)
				{
				state_space_size_Qmat_to_cladoMat = -1
				} else {
				state_space_size_Qmat_to_cladoMat = 0
				}

			# -1, assumes NULL range is allowed
# 			tmpca_1 = rep(1, sum(states_allowed_TF)-1)
# 			tmpcb_1 = rep(1, sum(states_allowed_TF)-1)
			tmpca_1 = rep(1, sum(states_allowed_TF)+state_space_size_Qmat_to_cladoMat)
			tmpcb_1 = rep(1, sum(states_allowed_TF)+state_space_size_Qmat_to_cladoMat)
			
			
			COO_weights_columnar = rcpp_calc_anclikes_sp_COOweights_faster(Rcpp_leftprobs=tmpca_1, Rcpp_rightprobs=tmpcb_1, l=l, s=maxent01s_param, v=maxent01v_param, j=maxent01j_param, y=maxent01y_param, dmat=dmat, maxent01s=maxent01s, maxent01v=maxent01v, maxent01j=maxent01j, maxent01y=maxent01y, max_minsize_as_function_of_ancsize=max_minsize_as_function_of_ancsize, printmat=FALSE, m=m)


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
			# UPPASS THROUGH THE TREE PIECES - CALCULATIONS
			# Go through the tree pieces in this stratum
			#######################################################
			chainsaw_result = inputs$tree_sections_list[[i]]
			
			# Set up a new list item to store uppass tip probs
			inputs$tree_sections_list[[i]]$pieces_relprobs_at_tips = list()

	
			# UPPASS: Go through tree pieces in this stratum (bottom first)
			for (jj in 1:length(chainsaw_result$return_pieces_list))
				{
				treepiece = chainsaw_result$return_pieces_list[[jj]]

				#cat("\ni=", i, "; jj=",jj, "; length(treepiece)=", length(treepiece), sep="")
				
				# If it's just a branch section
				if (is.numeric(treepiece) )
					{
					do_exponentiation = TRUE	# default

					# Also, exclude the case where there is a branch at the bottom below the bottom root node
					if (i == num_timeperiods)
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
					
					# 2016-02-29: rearranged
					# UPPASS:
					# THIS (FOSSIL_HASNT_OCCURRED_YET) MUST GO *BEFORE* THE 'FOSSIL_TIP_DOES_OCCUR_IN_BIN' IF/THEN
					# If this is TRUE, this fossil hasn't occured yet, and you are looking at the "phantom limb".
					# In this case, DON'T do matrix exponentiation, just copy the probabilities up!!
					if ( master_tip_time_bp < time_top )
						{
						do_exponentiation = FALSE
						}

					# UPPASS:
					# THIS (FOSSIL_TIP_DOES_OCCUR_IN_BIN) MUST GO *AFTER* THE 'FOSSIL_HASNT_OCCURRED_YET' IF/THEN
					# If this is TRUE, there's a match and the fossil tip appears in this time period
					if ( (master_tip_time_bp >= time_top) && (master_tip_time_bp < time_bot) && (is_fossil == TRUE))
						{
						# Shorten the branchlength by master_tip_time_bp-time_top
						amount_to_shorten_by = master_tip_time_bp-time_top
						subbranch_length = subbranch_length - amount_to_shorten_by
						do_exponentiation = TRUE
						}
					# If FALSE, you're below all this and hopefully don't care
					
					
					# 2016-02-29
					# UPPASS:
					# Also, DON'T do exponentiation if the branch length in the master branch
					# is a direct ancestor, i.e., less than min_branchlength
					if (tmp_master_table_row$edge.length < min_branchlength)
						{
						#print("It's a direct ancestor, so DON'T do matrix exponentiation!")
						# It's a direct ancestor, so DON'T do matrix exponentiation
						do_exponentiation = FALSE
						}


					
					# Get the uppass probs from the correct piece in the previous (below, older) stratum
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
					relprobs_at_branch_bottoms_below_tips_from_previous_stratum = inputs$tree_sections_list[[previous_stratum]]$pieces_relprobs_at_bottoms_below_tips[[previous_treepiece_num]]
					
					# If ancestor was a sub-branch
					if (is.numeric(previous_treepiece) == TRUE)
						{
						ancprobs_at_subbranch_bottom = relprobs_at_tips_of_anc_treepiece
						ancprobs_at_bottom_of_total_branch = relprobs_at_branch_bottoms_below_tips_from_previous_stratum
						} else {
						# Ancestor was a sub-tree
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
							
							if (traitTF == FALSE)
								{
								actual_probs_after_forward_exponentiation = calc_prob_forward_onebranch_dense(relprobs_branch_bottom=ancprobs_at_subbranch_bottom[states_allowed_TF], branch_length=subbranch_length, Qmat_tmp)
								} else {
								actual_probs_after_forward_exponentiation = calc_prob_forward_onebranch_dense(relprobs_branch_bottom=ancprobs_at_subbranch_bottom[wTrait_states_allowed_TF], branch_length=subbranch_length, Qmat_tmp)
								}
							
							if (include_null_range == TRUE)
								{
								# NULL range is impossible
								actual_probs_after_forward_exponentiation[1] = 0
								} # END if (include_null_range == TRUE)
								
							actual_probs_after_forward_exponentiation = actual_probs_after_forward_exponentiation / sum(actual_probs_after_forward_exponentiation)
							} else {
							# Sparse matrix exponentiation
#							print("this1")
							
							if (traitTF == FALSE)
								{
								actual_probs_after_forward_exponentiation = calc_prob_forward_onebranch_sparse(relprobs_branch_bottom=ancprobs_at_subbranch_bottom[states_allowed_TF], branch_length=subbranch_length, tmpQmat_in_REXPOKIT_coo_fmt=tmpQmat_in_REXPOKIT_coo_fmt, coo_n=numstates_geogtrait, anorm=NULL, check_for_0_rows=TRUE)
								} else {
								actual_probs_after_forward_exponentiation = calc_prob_forward_onebranch_sparse(relprobs_branch_bottom=ancprobs_at_subbranch_bottom[wTrait_states_allowed_TF], branch_length=subbranch_length, tmpQmat_in_REXPOKIT_coo_fmt=tmpQmat_in_REXPOKIT_coo_fmt, coo_n=numstates_geogtrait, anorm=NULL, check_for_0_rows=TRUE)
								}

							if (include_null_range == TRUE)
								{
								# NULL range is impossible
								actual_probs_after_forward_exponentiation[1] = 0
								} # END if (include_null_range == TRUE)

							actual_probs_after_forward_exponentiation = actual_probs_after_forward_exponentiation / sum(actual_probs_after_forward_exponentiation)
							}

						if (traitTF == FALSE)
							{
							# 2015 fix (and states_allowed_TF above)
							actual_probs_after_forward_exponentiation_new = rep(0, length(states_allowed_TF))
							actual_probs_after_forward_exponentiation_new[states_allowed_TF] = actual_probs_after_forward_exponentiation
							} else {
							actual_probs_after_forward_exponentiation_new = rep(0, length(wTrait_states_allowed_TF))
							actual_probs_after_forward_exponentiation_new[wTrait_states_allowed_TF] = actual_probs_after_forward_exponentiation							
							} # END if (traitTF == FALSE)
						actual_probs_after_forward_exponentiation = actual_probs_after_forward_exponentiation_new
					
						# Zero out impossible states in this zone (but NOT for "phantom limbs")
						# This CAN work, since they've been reset to main state space
						if (!is.null(states_allowed_TF))
							{
							actual_probs_after_forward_exponentiation[states_allowed_TF==FALSE] = 0
							actual_probs_after_forward_exponentiation = actual_probs_after_forward_exponentiation / sum(actual_probs_after_forward_exponentiation)
							}	
					
						} else {
						# Just pass up the ancestral probabilities, without modification, on the "phantom limb"
						# Don't do 2015 fix here, since this is a fossil branch and
						# we are just passing up the probabilities
						actual_probs_after_forward_exponentiation = ancprobs_at_subbranch_bottom
						}

					
					if (any(is.na(actual_probs_after_forward_exponentiation)))
						{
						print("i, jj, anc")
						print(i)
						print(jj)
						print(anc)
						print("actual_probs_after_forward_exponentiation")
						print(actual_probs_after_forward_exponentiation)
						print("ancprobs_at_subbranch_bottom")
						print(ancprobs_at_subbranch_bottom)
						print("ancprobs_at_bottom_of_total_branch")
						print(ancprobs_at_bottom_of_total_branch)
						stop("ERROR #1: see stratified code")
						}
						
					##########################################################
					# tip probabilities for next stratum up
					##########################################################
					relprobs_at_tips_for_next_stratum_up = actual_probs_after_forward_exponentiation
					relprobs_at_branch_bottoms_below_tips_for_next_stratum_up = ancprobs_at_bottom_of_total_branch
					
					# Store the relprobs at the tips, so that the next stratum up can access them...
					inputs$tree_sections_list[[i]]$pieces_relprobs_at_tips[[jj]] = relprobs_at_tips_for_next_stratum_up
					
					relprobs_at_tips_for_next_stratum_up
					
					inputs$tree_sections_list[[i]]$pieces_relprobs_at_bottoms_below_tips[[jj]] = relprobs_at_branch_bottoms_below_tips_for_next_stratum_up
					
					# If this is a tip in the master tree, store the tmp_relprobs_at_branchtop_AT_node_UPPASS in the master tree
					# (do this after the pass through the whole tree)
					
					# END if (is.numeric(treepiece) )
					} else {
					######################################################
					# UPPASS: Treepiece is a subtree!!
					tmp_subtree = treepiece
					######################################################
					
					######################################################
					# Check if this subtree contains a fixed internal node on the master tree
					######################################################
					# 2014-03-20_NJM
					# FIXNODES FOR UPPASS
					# NON-BUG BUG BUG (I don't think anyone has ever addressed UPPASS 
					# calculations with fixed ancestral states)
			
					use_fixnodes_on_uppass = TRUE	# turn off if desired/buggy
					if (use_fixnodes_on_uppass) {

					# Match the master fixnodes to the fixnodes in *JUST* this subtree
					# We will then pass these fixnodes to the subtree loglike calculation
					# First, we need to get the master node number, iff it's internal
					# 
					tmp_fixnode = NULL		# Default
					tmp_fixlikes = NULL		# Default
					if ((!is.null(fixnode)) && (length(fixnode) > 0))
						{
						# Check for multiple fixnodes
						if (length(fixnode) > 1)
							{
							# If there are multiple fixnodes, 
							# Get the matching nodes in this subtree
							TF1 = inputs$master_table$stratum == i
							TF2 = inputs$master_table$piecenum == jj
							TF3 = inputs$master_table$piececlass == "subtree"
							TF = ((TF1 + TF2 + TF3) == 3)
							tmprows = inputs$master_table[TF,]
						
							# Get the fixnodes found in this subtree
							fixnodes_in_subtree_TF = fixnode %in% tmprows$node
					
							# *IF* the subtree contains fixnodes, do this stuff
							# otherwise, set to NULL
							if (sum(fixnodes_in_subtree_TF) > 0)
								{
								#master_nodes_in_fixnode_TF = inputs$master_table$node %in% fixnode
								#master_nodes_in_fixnode
						
								#TF = (anc == fixnode)	# old
								# we do not use temporary_fixnode, since we need the fixnodes in the subtree numbering (tmprow$SUBnode)
								temporary_fixnodes = fixnode[fixnodes_in_subtree_TF]
								# But we will use these
								temporary_fixlikes = fixlikes[fixnodes_in_subtree_TF,]
						
								# The subtree nodenums corresponding to the subset temporary_fixnodes
								# NOTE! THIS SUBSET THING WILL ONLY WORK IF THE NODES ARE SORTED IN ORDER FROM THE START
								subtree_rows_in_fixnodes_TF = tmprows$node %in% fixnode 
								subtree_fixnode_master_nodenums = tmprows$node[subtree_rows_in_fixnodes_TF]
								subtree_fixnode_nums = tmprows$SUBnode[subtree_rows_in_fixnodes_TF]
							
								# We have to order these subtree fixnodes, and order the subtree fixlikes the same way
								order_subtree_fixnode_nums = order(subtree_fixnode_nums)
								subtree_fixnode_nums = subtree_fixnode_nums[order_subtree_fixnode_nums]
							
								# Only reorder if there are 2 or more rows, i.e. if it's a matrix not a vector
								if (length(order_subtree_fixnode_nums) > 1)
									{
									temporary_fixlikes = temporary_fixlikes[order_subtree_fixnode_nums, ]
									}
							
								} else {
								# If *NO* fixnodes in subtree:
								temporary_fixnodes = NULL
								subtree_fixnode_master_nodenums = NULL
								subtree_fixnode_nums = NULL
								temporary_fixlikes = NULL
								} # end if (sum(fixnodes_in_subtree_TF) > 0)


							# Check if we're in the right stratum / piece / piececlass
							# (have account for possible multiple rows)
							TF1 = unique(tmprows$stratum) == i
							TF2 = unique(tmprows$piecenum) == jj
							TF3 = unique(tmprows$piececlass) == "subtree"
							TF = ((TF1 + TF2 + TF3) == 3)
						
							if (TF == TRUE)
								{
								#txt = paste("Master tree node ", fixnode, " matched to i=", i, "; jj=", jj, "; piececlass=", piececlass, "; subtree subnode=", tmprow$SUBnode, sep="")
								#print(txt)
								#print(fixlikes)
				
								# Determine the number of the subnode in the subtree
								if (length(subtree_fixnode_nums) == 0)
									{
									subtree_fixnode_nums = NULL
									temporary_fixlikes = NULL
									}
						
								tmp_fixnode = subtree_fixnode_nums
								tmp_fixlikes = temporary_fixlikes
								} else {
								tmp_fixnode = NULL
								tmp_fixlikes = NULL
								} # end if (TF == TRUE)

							} else {
							# Only 1 fixnode
							temporary_fixnode = fixnode
							temporary_fixlikes = c(fixlikes)

							# e.g.
							# fixnode=20
							TF1 = inputs$master_table$node == temporary_fixnode
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
							
	# 							cat("\n\n")
	# 							print(fixnode)
	# 							print(temporary_fixnode)
	# 							print(tmprow$SUBnode)
	# 							print(temporary_fixlikes)
							

								# Determine the number of the subnode in the subtree
								tmp_fixnode = tmprow$SUBnode
								tmp_fixlikes = temporary_fixlikes
								# end if (TF == TRUE)
								} else {
								tmp_fixnode = NULL
								tmp_fixlikes = NULL
								}

							} # end if (length(fixnode) > 1)
						} # end if (!is.null(fixnode))
					} # end if (use_fixnodes_on_uppass)

					
					
					
					####################################################
					# Check for fossils on the UPPASS, shorten branches 
					# appropriately if found
					####################################################
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
							
							# 2016-02-29: this is now done in section_the_tree()
							# Edit the length of the branch on this subtree tip
							#tmp_subtree$edge.length[tmp2_edgenum] = tmp_subtree$edge.length[tmp2_edgenum] - amount_to_shorten_by
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
					
					
					# Get the SUBnode numbers of the subtree, and their corresponding nodenums in the 
					# master table
					
					if (traitTF == FALSE)
						{
						# 2014 version
						#tmp_relprobs_at_branchtop_AT_node_UPPASS = matrix(data=NA, nrow=numnodes, length(states_to_use_TF))
						#tmp_relprobs_at_branchbot_BELOW_node_UPPASS = matrix(data=NA, nrow=numnodes, length(states_to_use_TF))
						# 2015 version
						tmp_relprobs_at_branchtop_AT_node_UPPASS = matrix(data=0, nrow=numnodes, sum(states_allowed_TF))
						tmp_relprobs_at_branchbot_BELOW_node_UPPASS = matrix(data=0, nrow=numnodes, sum(states_allowed_TF))
						} else {
						tmp_relprobs_at_branchtop_AT_node_UPPASS = matrix(data=0, nrow=numnodes, sum(states_allowed_TF)*num_trait_states)
						tmp_relprobs_at_branchbot_BELOW_node_UPPASS = matrix(data=0, nrow=numnodes, sum(states_allowed_TF)*num_trait_states)
						} # END if (traitTF == FALSE)
					
					master_tree_nodenums = NULL
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
						TF5c = inputs$master_table$node.type == "root"		# Store only if node corresponds to a node or tip in orig_tree
						TF_subtrees = (TF1 + TF2 + TF3 + TF4 + TF5a + TF5b + TF5c) == 5
						
						master_tree_nodenums = c(master_tree_nodenums, inputs$master_table$node[TF_subtrees])
						} # END for (rownum in 1:nrow(tmp_relprobs_at_branchtop_AT_node_UPPASS))



					# Also, we need to get the saved DOWNPASS stuff at branch BOTTOMS, for the subtree uppass
#					print("dim(relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS)")
#					print(dim(relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS))
# 					print(master_tree_nodenums)
# 					print(states_allowed_TF)
# 					print(relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[master_tree_nodenums,])
# 					print(relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[,states_allowed_TF])
# 					print(relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[master_tree_nodenums,][,states_allowed_TF])
					if (traitTF == FALSE)
						{
						tmp_relprobs_at_branchbot_BELOW_node_DOWNPASS = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[master_tree_nodenums,][,states_allowed_TF]
						} else {
						tmp_relprobs_at_branchbot_BELOW_node_DOWNPASS = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[master_tree_nodenums,][,wTrait_states_allowed_TF]
						}
					
					# 2016-05-28_bug_fix
					# Fix this error, e.g. when DEC* model + areas_allowed means that
					# ranges_list = NULL, Kauai is just
					# ranges_list = Kauai
					# This means that:
					#   subtree_tip_relative_probs_of_each_state 
					# and thus
					#   tip_condlikes_of_data_on_each_state
					# ...are just a list of numbers, not a matrix, thus 
					# rowSums fails in calc_loglike_sp() in that time-stratum.
					# 
					#
					# This was the error message:
					# 
					# Error in rowSums(tip_condlikes_of_data_on_each_state) : 
					#  'x' must be an array of at least two dimensions
					# Calls: bears_optim_run ... calc_loglike_sp_stratified -> calc_loglike_sp -> rowSums
					#
					if (sum(states_allowed_TF) == 1)
						{
						if (traitTF == FALSE)
							{
							tmp_relprobs_at_branchbot_BELOW_node_DOWNPASS = matrix(data=tmp_relprobs_at_branchbot_BELOW_node_DOWNPASS, ncol=1)
							} else {
							tmp_relprobs_at_branchbot_BELOW_node_DOWNPASS = matrix(data=tmp_relprobs_at_branchbot_BELOW_node_DOWNPASS, ncol=sum(wTrait_states_allowed_TF))
							} # END if (traitTF == FALSE)
						} # END if (sum(states_allowed_TF) == 1)
					

					# i.e., if i==5 in the Psychotria dataset
					if (i == num_timeperiods) # You are at the bottom tree piece, just use root node
						{
						#ancprobs_at_subtree_root = starting_probs

						# Anc node of the subtree
						if (traitTF == FALSE)
							{
							tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, ] = relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[anc_node_original_tree, ][states_allowed_TF]
							} else {
							#print("Here1")
							#print(tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, ])
							#print(relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[anc_node_original_tree, ][wTrait_states_allowed_TF])
							#print(wTrait_states_allowed_TF)
							tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, ] = relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[anc_node_original_tree, ][wTrait_states_allowed_TF]
							} # END if (traitTF == FALSE)
						
						# None of this, at the root
						# NO: tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc, ] = ancprobs_at_bottom_of_total_branch
						
						} else {
						# You are NOT at the bottom tree piece
						if ((is.numeric(phy2$root.edge) == TRUE) && (!is.null(phy2$root.edge)) && (phy2$root.edge > 0))
							{
							# Get the length of this branch
							root_edge_length = phy2$root.edge
							
							# Get the uppass probs from the correct piece in the previous stratum
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
								if (traitTF == FALSE)
									{
									actual_probs_after_forward_exponentiation = calc_prob_forward_onebranch_dense(relprobs_branch_bottom=ancprobs_at_subbranch_bottom[states_allowed_TF], branch_length=root_edge_length, Qmat_tmp)
									} else {
									actual_probs_after_forward_exponentiation = calc_prob_forward_onebranch_dense(relprobs_branch_bottom=ancprobs_at_subbranch_bottom[wTrait_states_allowed_TF], branch_length=root_edge_length, Qmat_tmp)
									} # END if (traitTF == FALSE)
								} else {
								# Sparse matrix exponentiation
# 								print(tmpQmat_in_REXPOKIT_coo_fmt)
# 								print("this2")
# 								print("numstates_geogtrait")
# 								print(numstates_geogtrait)
# 								print(ancprobs_at_subbranch_bottom[states_allowed_TF])
								
								if (traitTF == FALSE)
									{
									actual_probs_after_forward_exponentiation = calc_prob_forward_onebranch_sparse(relprobs_branch_bottom=ancprobs_at_subbranch_bottom[states_allowed_TF], branch_length=root_edge_length, tmpQmat_in_REXPOKIT_coo_fmt=tmpQmat_in_REXPOKIT_coo_fmt, coo_n=numstates_geogtrait, anorm=NULL, check_for_0_rows=TRUE)
									} else {
									actual_probs_after_forward_exponentiation = calc_prob_forward_onebranch_sparse(relprobs_branch_bottom=ancprobs_at_subbranch_bottom[wTrait_states_allowed_TF], branch_length=root_edge_length, tmpQmat_in_REXPOKIT_coo_fmt=tmpQmat_in_REXPOKIT_coo_fmt, coo_n=numstates_geogtrait, anorm=NULL, check_for_0_rows=TRUE)
									} # END if (traitTF == FALSE)
								} # END if (sparse==FALSE)

							
							if (traitTF == FALSE)
								{
								# 2015 fix
								actual_probs_after_forward_exponentiation_new = rep(0, length(states_allowed_TF))
								actual_probs_after_forward_exponentiation_new[states_allowed_TF] = actual_probs_after_forward_exponentiation
								actual_probs_after_forward_exponentiation = actual_probs_after_forward_exponentiation_new
								} else {
								actual_probs_after_forward_exponentiation_new = rep(0, length(wTrait_states_allowed_TF))
								actual_probs_after_forward_exponentiation_new[wTrait_states_allowed_TF] = actual_probs_after_forward_exponentiation
								actual_probs_after_forward_exponentiation = actual_probs_after_forward_exponentiation_new								
								}
							
							if (include_null_range == TRUE)
								{
								# NULL range is impossible
								actual_probs_after_forward_exponentiation[1] = 0
								} # END if (include_null_range == TRUE)

							actual_probs_after_forward_exponentiation = actual_probs_after_forward_exponentiation / sum(actual_probs_after_forward_exponentiation)

							# This CAN work, since they've been reset to main state space
							# Zero out impossible states
							if (!is.null(states_allowed_TF))
								{
								actual_probs_after_forward_exponentiation[states_allowed_TF==FALSE] = 0
								actual_probs_after_forward_exponentiation = actual_probs_after_forward_exponentiation / sum(actual_probs_after_forward_exponentiation)
								}	

							# 2015 fix:ancprobs_at_subtree_root is for FULL state space
							ancprobs_at_subtree_root = actual_probs_after_forward_exponentiation
							
							if (traitTF == FALSE)
								{
								# But this is for reduced state space
								tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, ] = actual_probs_after_forward_exponentiation[states_allowed_TF]
								# This is also for reduced state space
								tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc, ] = ancprobs_at_bottom_of_total_branch[states_allowed_TF]
								} else {
								# But this is for reduced state space
								tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, ] = actual_probs_after_forward_exponentiation[wTrait_states_allowed_TF]
								# This is also for reduced state space
								tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc, ] = ancprobs_at_bottom_of_total_branch[wTrait_states_allowed_TF]								
								}
							} else {

							if (traitTF == FALSE)
								{
								# No trait
								# No root edge; just use probs at anc of subtree
								# 2015 reduced state space
								tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, ] = relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[anc_node_original_tree, ][states_allowed_TF]
								# 2015 reduced state space
								tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc, ] = ancprobs_at_bottom_of_total_branch[states_allowed_TF]
								} else {
								# With trait
								#print("Here2")
								#print(tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, ])
								#print(relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[anc_node_original_tree, ][wTrait_states_allowed_TF])
								#print(wTrait_states_allowed_TF)
								tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, ] = relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[anc_node_original_tree, ][wTrait_states_allowed_TF]
								tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc, ] = ancprobs_at_bottom_of_total_branch[wTrait_states_allowed_TF]
								} # END if (traitTF == FALSE)

							} # END if ((is.numeric(phy2$root.edge) == TRUE) && (!is.null(phy2$root.edge)) && (phy2$root.edge > 0))
							# END check for root edge
						
						
						if (any(is.na(tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, ])))
							{
							print("i, jj, anc")
							print(i)
							print(jj)
							print(anc)
							print("actual_probs_after_forward_exponentiation")
							print(actual_probs_after_forward_exponentiation)
							print("tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, ]")
							print(tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, ])		
							print("anc_node_original_tree")
							print(anc_node_original_tree)
							print("states_allowed_TF")
							print(states_allowed_TF)
							if (traitTF == FALSE)
								{
								print("ancprobs_at_subbranch_bottom[states_allowed_TF]")
								print(ancprobs_at_subbranch_bottom[states_allowed_TF])
								} else {
								print("wTrait_states_allowed_TF")
								print(wTrait_states_allowed_TF)
								print("ancprobs_at_subbranch_bottom[wTrait_states_allowed_TF]")
								print(ancprobs_at_subbranch_bottom[wTrait_states_allowed_TF])
								}
							stop("ERROR #2: see stratified code")
							}
						if (any(is.na(tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc, ])))
							{
							print("i, jj, anc")
							print(i)
							print(jj)
							print(anc)
							print("actual_probs_after_forward_exponentiation")
							print(actual_probs_after_forward_exponentiation)
							print("ancprobs_at_bottom_of_total_branch")
							print(ancprobs_at_bottom_of_total_branch)
							print("tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc, ]")
							print(tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc, ])
							stop("ERROR #3: see stratified code")
							}
						
						} # END if (i == num_timeperiods) # You are at the bottom tree piece, just use root node




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
					# The root node of this subtree
					rootnode = length(phy2$tip.label) + 1
					
					for (uj in edges_to_visit_uppass)		# Since we are going backwards
						{
						# First edge visited is ui
						#print(ui)

						#print("Qmat_tmp")						
						#print(Qmat_tmp)
						#print(dim(Qmat_tmp))
						
						# Its sister is uj 
						#uj <- ui - 1
						ui <- uj - 1		# Since we are going backwards
						
						# Get the node numbers at the tips of these two edges		
						left_desc_nodenum <- phy2$edge[ui, 2]
						right_desc_nodenum <- phy2$edge[uj, 2]
			
						# And for the ancestor edge (i or j shouldn't matter, should produce the same result!!!)
						anc <- phy2$edge[ui, 1]
						# ancedge
						anc_edgenum_TF = phy2$edge[,2] == anc
						anc_edgenum = (1:length(anc_edgenum_TF))[anc_edgenum_TF]
			
						# For the marginal state probability, uppass calculations
						# Get the mother and sister of "anc" (which is the focal node)
						mother_of_anc_TF = phy2$edge[,2] == anc
						mother_of_anc = phy2$edge[mother_of_anc_TF,1]

						sister_of_anc_TF = phy2$edge[,1] == mother_of_anc
						sister_of_anc_TF2 = (sister_of_anc_TF + mother_of_anc_TF) == 1
						sister_of_anc = phy2$edge[sister_of_anc_TF2,2]
						mother_of_anc
						sister_of_anc
			
						# Is the sister left or right?
						# (note: these are reversed from what you would get with:
						#  plot(tr); nodelabels()
						sister_is_LR = "rootnode"
						if (anc != rootnode)
							{
							if (sister_of_anc > anc)
								{
								sister_is_LR = "right"
								} else {
								sister_is_LR = "left"
								} # END if (sister_of_anc > anc)
							} # END if (anc != rootnode)
						
						# get the correct edges
						left_edge_TF = phy2$edge[,2] == left_desc_nodenum
						right_edge_TF = phy2$edge[,2] == right_desc_nodenum
						left_edgenum = (1:length(left_edge_TF))[left_edge_TF]
						right_edgenum = (1:length(right_edge_TF))[right_edge_TF]

						# Check the branchlength of each edge
						# It's a hook if either branch is super-short
						is_leftbranch_hook_TF = phy2$edge.length[left_edge_TF] < min_branchlength
						is_rightbranch_hook_TF = phy2$edge.length[right_edge_TF] < min_branchlength
						hooknode_TF = (is_leftbranch_hook_TF + is_rightbranch_hook_TF) > 0
			
			
						#cat(i, j, left_desc_nodenum, right_desc_nodenum, hooknode_TF, "\n", sep="	")
						
						
						# You start with these uppass probs, for this node
						tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, ]
						#print("tmp_relprobs_at_branchtop_AT_node_UPPASS")
						#print(dim(tmp_relprobs_at_branchtop_AT_node_UPPASS))
						#print("tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, ]")


						# 2014 version
						#numstates = ncol(tip_condlikes_of_data_on_each_state)
						# 2015 version
						sum_states_allowed = sum(states_allowed_TF)
						if (traitTF == TRUE)
							{
							wTrait_sum_states_allowed = sum(states_allowed_TF) * num_trait_states
							} # END if (traitTF == TRUE)
						#print("states_allowed_TF")
						#print(states_allowed_TF)
						
						# Apply speciation model to get the uppass probs at the base of the two descendant branches
						if (hooknode_TF == TRUE)
							{
							# Just copy the probs up, since a time-continuous model was assumed.
							# If you have a "hooknode" (short branch = direct ancestor), for
							# the uppass, it is simpler to convert the cladogenesis model
							# to an all-1s model
							temp_COO_weights_columnar = COO_weights_columnar

							# NJM 2016-02-24 -- see:
							# /drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/AAAB_M3_ancestor_check
							# ...for an example that traces the -2 issue in detail. 
							# Or:
							# http://phylo.wikidot.com/fossil-data-in-biogeographical-analysis-in-biogeobears#toc3
							#
							# Basically, this converts 1-based state numbers (e.g., 1-16, with
							# 1: null range
							# 2: A
							# 3: B
							# ...etc..
							#
							# ...to the 0-based state names in a cladogenesis matrix, where the
							# null range is automatically excluded (if used in the first place).
							#
							# Then the states in the cladogenesis matrix are numbered, starting from 0.
							# E.g.:
							# 0: A
							# 1: B
							# 2: C
							# ...etc...
				
							if (include_null_range == TRUE)
								{
								#highest_clado_state_0based_considering_null_range = numstates - 1
								highest_clado_state_0based_considering_null_range = sum_states_allowed - 2
								} else {
								#highest_clado_state_0based_considering_null_range = numstates
								highest_clado_state_0based_considering_null_range = sum_states_allowed - 1
								} # if (include_null_range == TRUE)

							# If you have a "hooknode" (short branch = direct ancestor), for
							# the uppass, it is simpler to convert the cladogenesis model
							# to an all-1s model

							# Ancestral, left, and right states all the same
							temp_COO_weights_columnar[[1]] = 0:highest_clado_state_0based_considering_null_range
							temp_COO_weights_columnar[[2]] = 0:highest_clado_state_0based_considering_null_range
							temp_COO_weights_columnar[[3]] = 0:highest_clado_state_0based_considering_null_range
							temp_COO_weights_columnar[[4]] = rep(1, highest_clado_state_0based_considering_null_range+1)
							} else {
							temp_COO_weights_columnar = COO_weights_columnar
							} # END if (hooknode_TF == TRUE)
						
						#print("temp_COO_weights_columnar")
						#print(temp_COO_weights_columnar)
						
						##############################################################################################
						# Apply regular speciation model, with the weights given in COO_weights_columnar, and the 
						# normalization factor (sum of the weights across each row/ancestral state) in Rsp_rowsums.
						##############################################################################################
						num_nonzero_split_scenarios = length(COO_weights_columnar[[1]])
						

						# Probs at the mother have been predetermined, in the uppass
						# 1. Get uppass probabilities at the base of the branch below the 
						#    focal (anc) node, including the probabilities coming down
						#    from the sister, and up from the mother.
						
						# Check if you are at the global root
						TFi = inputs$master_table$stratum == i
						TFjj = inputs$master_table$piecenum == jj
						TF_SUBnode = inputs$master_table$SUBnode == anc
						TF = ((TFi + TFjj + TF_SUBnode) == 3)
						anc_node_original_tree = inputs$master_table$node[TF]
						global_root_TF = inputs$master_table$node.type[TF]

						if ((anc == rootnode) && (global_root_TF == TRUE))
						#if (anc == rootnode)
							{
							# You ARE at the global ancestor (root) node
							probs_at_mother = 1/length(starting_probs)
							likes_at_sister = 1/length(starting_probs)
							left_branch_downpass_likes = NULL
							right_branch_downpass_likes = NULL
							#tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc,] = NA		
							probs_of_mother_and_sister_uppass_to_anc = 1/length(starting_probs)
							tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc,] = probs_of_mother_and_sister_uppass_to_anc
							} else {
							# You ARE NOT at the global ancestor (root) node
							if (anc == rootnode)
								{
								probs_at_mother = 1/length(starting_probs)
								probs_of_mother_and_sister_uppass_to_anc = tmp_relprobs_at_branchtop_AT_node_UPPASS[anc,]
								} else {
								probs_at_mother = tmp_relprobs_at_branchtop_AT_node_UPPASS[mother_of_anc,]
								} # END if (anc == rootnode)
							likes_at_sister_branch_bottom = tmp_relprobs_at_branchbot_BELOW_node_DOWNPASS[sister_of_anc,]

							if (sister_is_LR == "left")
								{	
								left_branch_downpass_likes = likes_at_sister_branch_bottom
								right_branch_downpass_likes = NULL
								}
							if (sister_is_LR == "right")
								{	
								left_branch_downpass_likes = NULL
								right_branch_downpass_likes = likes_at_sister_branch_bottom
								}
				
							# Calculate the uppass probs at the branch 	
							#print("calculation of uppass at split")
							#print(probs_at_mother)
							#print(left_branch_downpass_likes)
							#print(right_branch_downpass_likes)
							
							if (anc != rootnode)
								{
								#print("Here123")
								#print("probs_at_mother")
								#print(probs_at_mother)
								#print("left_branch_downpass_likes")
								#print(left_branch_downpass_likes)
								#print("right_branch_downpass_likes")
								#print(right_branch_downpass_likes)
								
								if (traitTF == FALSE)
									{
									uppass_probs_at_bottom_below_anc_results = calc_uppass_probs_new2(probs_ancstate=probs_at_mother, COO_weights_columnar=temp_COO_weights_columnar, numstates=sum_states_allowed, include_null_range=include_null_range, left_branch_downpass_likes=left_branch_downpass_likes, right_branch_downpass_likes=right_branch_downpass_likes, Rsp_rowsums=NULL)
									} else {
									uppass_probs_at_bottom_below_anc_results = calc_uppass_probs_new2(probs_ancstate=probs_at_mother, COO_weights_columnar=temp_COO_weights_columnar, numstates=wTrait_sum_states_allowed, include_null_range=include_null_range, left_branch_downpass_likes=left_branch_downpass_likes, right_branch_downpass_likes=right_branch_downpass_likes, Rsp_rowsums=NULL)
									} # END if (traitTF == FALSE)
								#print("...finished uppass at split")

								# Store
								if (sister_is_LR == "left")
									{	
									Rprobs_brbot_below_anc = uppass_probs_at_bottom_below_anc_results$relprobs_just_after_speciation_UPPASS_Right
									tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc,] = Rprobs_brbot_below_anc
									}
								if (sister_is_LR == "right")
									{
									Lprobs_brbot_below_anc = uppass_probs_at_bottom_below_anc_results$relprobs_just_after_speciation_UPPASS_Left
									tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc,] = Lprobs_brbot_below_anc
									}
				
								# 2. Exponentiate up from the mother to the focal/anc node
								probs_at_branch_bot = tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc,]
								
								if (force_sparse == FALSE)
									{
									probs_of_mother_and_sister_uppass_to_anc = probs_at_branch_bot %*% expokit_dgpadm_Qmat2(times=phy2$edge.length[anc_edgenum], Qmat=Qmat_tmp, transpose_needed=TRUE)
									} else {
									probs_of_mother_and_sister_uppass_to_anc = calc_prob_forward_onebranch_sparse(relprobs_branch_bottom=probs_at_branch_bot, branch_length=phy2$edge.length[anc_edgenum], tmpQmat_in_REXPOKIT_coo_fmt=tmpQmat_in_REXPOKIT_coo_fmt, coo_n=length(probs_at_branch_bot), anorm=NULL, check_for_0_rows=TRUE)
									}
								} else {
								# Subtree rootnode, so, the uppass probability was already
								# determined in processing the root branch of the subtree
								probs_of_mother_and_sister_uppass_to_anc
								} # END if (anc != rootnode)
				
							# Store in uppass
							tmp_relprobs_at_branchtop_AT_node_UPPASS[anc,] = probs_of_mother_and_sister_uppass_to_anc
							} # END if (anc == rootnode)

						##################################################################
						# Finish uppass to tips
						##################################################################
						# Check if either the left or right descendant nodes are tips;
						# if so, do the exponentiation here, so as to completely fill
						# in the UPPASS table
						##################################################################
						# If Left descendant is a tip
						if (left_desc_nodenum <= length(phy2$tip.label))
							{
							probs_at_anc = tmp_relprobs_at_branchtop_AT_node_UPPASS[anc,]
							left_branch_downpass_likes = NULL
							right_branch_downpass_likes = tmp_relprobs_at_branchbot_BELOW_node_DOWNPASS[right_desc_nodenum,]
							
							if (traitTF == FALSE)
								{
								uppass_probs_at_bottom_below_tip_results = calc_uppass_probs_new2(probs_ancstate=probs_at_anc, COO_weights_columnar=temp_COO_weights_columnar, numstates=sum_states_allowed, include_null_range=include_null_range, left_branch_downpass_likes=left_branch_downpass_likes, right_branch_downpass_likes=right_branch_downpass_likes, Rsp_rowsums=NULL)
								} else {
								uppass_probs_at_bottom_below_tip_results = calc_uppass_probs_new2(probs_ancstate=probs_at_anc, COO_weights_columnar=temp_COO_weights_columnar, numstates=wTrait_sum_states_allowed, include_null_range=include_null_range, left_branch_downpass_likes=left_branch_downpass_likes, right_branch_downpass_likes=right_branch_downpass_likes, Rsp_rowsums=NULL)
								}
							# The UPPASS probabilities below the tip:
							Lprobs_brbot_below_tip = uppass_probs_at_bottom_below_tip_results$relprobs_just_after_speciation_UPPASS_Left
							
							#print(dim(Qmat_tmp))
							#print(Lprobs_brbot_below_tip)
							#print(dim(Lprobs_brbot_below_tip))
							#print(length(Lprobs_brbot_below_tip))
							
							# The UPPASS probabilities AT the tip:
							if (force_sparse == FALSE)
								{
								Lprobs_brtop_AT_tip = Lprobs_brbot_below_tip %*% expokit_dgpadm_Qmat2(times=phy2$edge.length[left_edgenum], Qmat=Qmat_tmp, transpose_needed=TRUE)
								} else {
								# Sparse matrix exponentiation
								Lprobs_brtop_AT_tip = calc_prob_forward_onebranch_sparse(relprobs_branch_bottom=Lprobs_brbot_below_tip, branch_length=phy2$edge.length[left_edgenum], tmpQmat_in_REXPOKIT_coo_fmt=tmpQmat_in_REXPOKIT_coo_fmt, coo_n=numstates_geogtrait, anorm=NULL, check_for_0_rows=TRUE)
								}

				
							# Store: branch bottoms
							tmp_relprobs_at_branchbot_BELOW_node_UPPASS[left_desc_nodenum,] = Lprobs_brbot_below_tip
							# Store: branch tops
							tmp_relprobs_at_branchtop_AT_node_UPPASS[left_desc_nodenum,] = Lprobs_brtop_AT_tip
							} # END if (left_desc_nodenum <= length(phy2$tip.label))
			
						# If Right descendant is a tip
						if (right_desc_nodenum <= length(phy2$tip.label))
							{
							probs_at_anc = tmp_relprobs_at_branchtop_AT_node_UPPASS[anc,]
							right_branch_downpass_likes = NULL
							
# 							print("left_desc_nodenum:")
# 							print(left_desc_nodenum)
# 							
# 							print("tmp_relprobs_at_branchbot_BELOW_node_DOWNPASS:")
# 							print(tmp_relprobs_at_branchbot_BELOW_node_DOWNPASS)
# 							
# 							print("dim(tmp_relprobs_at_branchbot_BELOW_node_DOWNPASS):")
# 							print(dim(tmp_relprobs_at_branchbot_BELOW_node_DOWNPASS))
# 							
# 							print("states_allowed_TF:")
# 							print(states_allowed_TF)
# 							print("length(states_allowed_TF):")
# 							print(length(states_allowed_TF))
# 							
# 							print("c(uj, ui, jj, i):")
# 							print(c(uj, ui, jj, i))

							left_branch_downpass_likes = tmp_relprobs_at_branchbot_BELOW_node_DOWNPASS[left_desc_nodenum,]
							
							if (traitTF == FALSE)
								{
								uppass_probs_at_bottom_below_tip_results = calc_uppass_probs_new2(probs_ancstate=probs_at_anc, COO_weights_columnar=temp_COO_weights_columnar, numstates=sum_states_allowed, include_null_range=include_null_range, right_branch_downpass_likes=right_branch_downpass_likes, left_branch_downpass_likes=left_branch_downpass_likes, Rsp_rowsums=NULL)
								} else {
								uppass_probs_at_bottom_below_tip_results = calc_uppass_probs_new2(probs_ancstate=probs_at_anc, COO_weights_columnar=temp_COO_weights_columnar, numstates=wTrait_sum_states_allowed, include_null_range=include_null_range, right_branch_downpass_likes=right_branch_downpass_likes, left_branch_downpass_likes=left_branch_downpass_likes, Rsp_rowsums=NULL)
								}
				
							# The UPPASS probabilities below the tip:
							Rprobs_brbot_below_tip = uppass_probs_at_bottom_below_tip_results$relprobs_just_after_speciation_UPPASS_Right

							#print(dim(Qmat_tmp))
							#print(Rprobs_brbot_below_tip)
							#print(dim(Rprobs_brbot_below_tip))
							#print(length(Rprobs_brbot_below_tip))

							# The UPPASS probabilities AT the tip:
							if (force_sparse == FALSE)
								{
								Rprobs_brtop_AT_tip = Rprobs_brbot_below_tip %*% expokit_dgpadm_Qmat2(times=phy2$edge.length[right_edgenum], Qmat=Qmat_tmp, transpose_needed=TRUE)
								} else {
								# Sparse matrix exponentiation
								Rprobs_brtop_AT_tip = calc_prob_forward_onebranch_sparse(relprobs_branch_bottom=Rprobs_brbot_below_tip, branch_length=phy2$edge.length[right_edgenum], tmpQmat_in_REXPOKIT_coo_fmt=tmpQmat_in_REXPOKIT_coo_fmt, coo_n=numstates_geogtrait, anorm=NULL, check_for_0_rows=TRUE)
								}
				
							# Store: branch bottoms
							tmp_relprobs_at_branchbot_BELOW_node_UPPASS[right_desc_nodenum,] = Rprobs_brbot_below_tip
							# Store: branch tops
							tmp_relprobs_at_branchtop_AT_node_UPPASS[right_desc_nodenum,] = Rprobs_brtop_AT_tip
							} # END if (right_desc_nodenum <= length(phy2$tip.label))				

						##################################################################
						# END finish uppass to tips
						##################################################################
			
						# Zero out impossible states
						# Do NOT do this, since we are in the subset state space
# 						if (!is.null(states_allowed_TF))
# 							{
# 			tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc,][states_allowed_TF==FALSE] = 0
# 			tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc,][states_allowed_TF==FALSE] = 0
# 
# 			tmp_relprobs_at_branchbot_BELOW_node_UPPASS[left_desc_nodenum,][states_allowed_TF==FALSE] = 0
# 			tmp_relprobs_at_branchbot_BELOW_node_UPPASS[right_desc_nodenum,][states_allowed_TF==FALSE] = 0
# 
# 			tmp_relprobs_at_branchbot_BELOW_node_UPPASS[left_desc_nodenum,][states_allowed_TF==FALSE] = 0				
# 			tmp_relprobs_at_branchbot_BELOW_node_UPPASS[right_desc_nodenum,][states_allowed_TF==FALSE] = 0
# 
# 							} # END if (!is.null(states_allowed_TF))	
			
			
						# 2014-03-20_NJM
						# FIXNODES FOR UPPASS
						# NON-BUG BUG BUG (I don't think anyone has ever addressed UPPASS 
						# calculations with fixed ancestral states)
			
						use_fixnodes_on_uppass = TRUE	# turn off if desired/buggy
						if (use_fixnodes_on_uppass)
							{
							#######################################################
							# If the states/likelihoods have been fixed at a particular node
							# (check top of anc branch)
							#######################################################
							if (!is.null(fixnode))
								{
								# For multiple fixnodes
								# 2016-03-15_old
								# if (length(fixnode) > 1)
								# 2016-03-15_new by Torsten
								if (length(tmp_fixnode) > 1)
									{
									# Get the matching node
									TF = (anc == tmp_fixnode)
									temporary_fixnode = tmp_fixnode[TF]
									temporary_fixlikes = c(tmp_fixlikes[TF,])
									} else {
									temporary_fixnode = tmp_fixnode
									temporary_fixlikes = c(tmp_fixlikes)
									}

				
								if ((length(temporary_fixnode) > 0) && (anc == temporary_fixnode))
									{
									# If the node is fixed, ignore the calculation for this node, and
									# instead use the fixed likelihoods (i.e., the "known" state) for
									# this node.
									# fix the likelihoods of the (NON-NULL) states
# 2016-03-15_old
# 									tmp_relprobs_at_branchtop_AT_node_UPPASS[anc,] = tmp_relprobs_at_branchtop_AT_node_UPPASS[anc,] * temporary_fixlikes
# 2016-03-15_new by Torsten
 									if (traitTF == FALSE)
 										{
	 									tmp_relprobs_at_branchtop_AT_node_UPPASS[anc,] =  tmp_relprobs_at_branchtop_AT_node_UPPASS[anc,] * temporary_fixlikes[states_allowed_TF]
	 									} else {
	 									tmp_relprobs_at_branchtop_AT_node_UPPASS[anc,] =  tmp_relprobs_at_branchtop_AT_node_UPPASS[anc,] * temporary_fixlikes[wTrait_states_allowed_TF]
	 									}
									}
								} # end if (!is.null(fixnode))
							} # end if (use_fixnodes_on_uppass)

						# Normalize and save these probabilities
						#tmp_relprobs_at_branchtop_AT_node_UPPASS[left_desc_nodenum,] = condprobs_Left_branch_top / sum(condprobs_Left_branch_top)
						#tmp_relprobs_at_branchtop_AT_node_UPPASS[right_desc_nodenum,] = condprobs_Right_branch_top / sum(condprobs_Right_branch_top)

						#######################################################
						# End of UPPASS loop for this ancnode.  Move to next ancnode.
						#######################################################
						} # End uppass loop



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
						TF5c = inputs$master_table$node.type == "root"		# Store only if node corresponds to a node or tip in orig_tree
						TF_subtrees = (TF1 + TF2 + TF3 + TF4 + TF5a + TF5b + TF5c) == 5
						
						# (2) Tips of subbranches that are also tips of the master tree
						# (see below)
						
						TF = TF_subtrees
						
						
						# Store in the FINAL table
						relative_probs_of_each_state_at_branch_top_AT_node_UPPASS_rownum = inputs$master_table$node[TF]
						
						if (traitTF == FALSE)
							{	relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[relative_probs_of_each_state_at_branch_top_AT_node_UPPASS_rownum, ][states_allowed_TF] = tmp_relprobs
							} else {
							relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[relative_probs_of_each_state_at_branch_top_AT_node_UPPASS_rownum, ][wTrait_states_allowed_TF] = tmp_relprobs
							} # END if (traitTF == FALSE)
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
						
						if (traitTF == FALSE)
							{	relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS[relative_probs_of_each_state_at_branch_top_AT_node_UPPASS_rownum, ][states_allowed_TF] = tmp_relprobs
							} else {
relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS[relative_probs_of_each_state_at_branch_top_AT_node_UPPASS_rownum, ][wTrait_states_allowed_TF] = tmp_relprobs
							}

						} # END Store the UPPASS relprobs in the main matrix

					##########################################################
					# tip probabilities for next stratum up
					##########################################################
					if (traitTF == FALSE)
						{
						relprobs_at_tips_for_next_stratum_up = matrix(0, nrow=length(phy2$tip.label), ncol=length(states_allowed_TF))
						relprobs_at_tips_for_next_stratum_up[,states_allowed_TF] = tmp_relprobs_at_branchtop_AT_node_UPPASS[1:length(phy2$tip.label), ]
					
						# Branch bottoms below tips can also be transferred up
						relprobs_at_branch_bottoms_below_tips_for_next_stratum_up = matrix(0, nrow=length(phy2$tip.label), ncol=length(states_allowed_TF))
						relprobs_at_branch_bottoms_below_tips_for_next_stratum_up[,states_allowed_TF] = tmp_relprobs_at_branchbot_BELOW_node_UPPASS[1:length(phy2$tip.label), ]
						} else {
						relprobs_at_tips_for_next_stratum_up = matrix(0, nrow=length(phy2$tip.label), ncol=length(wTrait_states_allowed_TF))
						relprobs_at_tips_for_next_stratum_up[,wTrait_states_allowed_TF] = tmp_relprobs_at_branchtop_AT_node_UPPASS[1:length(phy2$tip.label), ]
					
						# Branch bottoms below tips can also be transferred up
						relprobs_at_branch_bottoms_below_tips_for_next_stratum_up = matrix(0, nrow=length(phy2$tip.label), ncol=length(wTrait_states_allowed_TF))
						relprobs_at_branch_bottoms_below_tips_for_next_stratum_up[,wTrait_states_allowed_TF] = tmp_relprobs_at_branchbot_BELOW_node_UPPASS[1:length(phy2$tip.label), ]
						} # END if (traitTF == FALSE)
					
					
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
			
			if (any(is.na(relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[tn,])))
				{
				print("i, jj, anc, tn")
				print(i)
				print(jj)
				print(anc)
				print(tn)
				print("relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[tn,]")
				print(relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[tn,])
				print("tmp_tipprobs_at_top_UPPASS")
				print(tmp_tipprobs_at_top_UPPASS)	
				print("tmprow2$piececlass")
				print(tmprow2$piececlass)
				stop("ERROR #4: see stratified code")
				}
			
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

		# Downpass conditional likelihoods at nodes and subtree tips
		# (master tree tip likelihoods in stratum 0 of master_table)
		calc_loglike_sp_stratified_results$condlikes_table = condlikes_table
		
		if (calc_ancprobs == TRUE)
			{
			# Relative downpass probabilities (rescaled conditional likelihoods)
			# at the BOTTOMS of branches of subtrees
			# 2014-05-26_NJM: now INCLUDES the probabilities at the BOTTOM of the subtree ROOT BRANCH
			# 2014-05-26_NJM: AND the probabilities at the bottom of sub-branches (except for the master_tree tip nodes,
			# which have tip likelihoods)
			calc_loglike_sp_stratified_results$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE

			# Uppass probabilities
			calc_loglike_sp_stratified_results$relative_probs_of_each_state_at_branch_top_AT_node_UPPASS = relative_probs_of_each_state_at_branch_top_AT_node_UPPASS
			calc_loglike_sp_stratified_results$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS = relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS
			
			# tmpres$relative_probs_of_each_state_at_branch_top_AT_node_UPPASS
			# tmpres$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS
	
	
			#######################################################
			# For branch bottoms
			#######################################################
			ML_marginal_prob_each_state_at_branch_bottom_below_node = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS * relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS
	
			# NJM - 2015-01-06: But, the root probabilities 
			# should NOT be multiplied, that would
			# be downpass * downpass, which 
			# results in focusing probability on the 
			# most-probable downpass state
			
			# Get the root node
			anc_row_of_master_table_TF = inputs$master_table$node.type=="root"
			anc_node_original_tree = inputs$master_table$node[anc_row_of_master_table_TF]
			anc_node_original_tree
			
			ML_marginal_prob_each_state_at_branch_bottom_below_node[anc_node_original_tree,] = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[anc_node_original_tree, ]	
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
				"This might occur in a highly constrained model, or if your data strongly contradicts your manual fixed\n",
				"likelihoods ('fixlikes') at some node(s) ('fixnode').\n",
				"As a 'fix', the downpass probabilities are being used for those nodes. But this is NOT RECOMMENDED!\n",
				"You should instead figure out what is causing the problem.", sep="")
				cat(txt)

				cat("\n\nPrinting (partial) downpass, uppass, and probability matrices to screen:\n\n", sep="") 
				print(relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS)				
				print(relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS)
				print(ML_marginal_prob_each_state_at_branch_bottom_below_node)

				ML_marginal_prob_each_state_at_branch_bottom_below_node[NaN_TF,] = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[NaN_TF,]
				}
	
	
	
			#######################################################
			# For state probabilities at branch tops
			#######################################################
			# State estimates under specified model (usually global ML)
			# are downpass * uppass
			ML_marginal_prob_each_state_at_branch_top_AT_node = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS * relative_probs_of_each_state_at_branch_top_AT_node_UPPASS
			
			# NJM - 2015-01-06: But, the root probabilities 
			# should NOT be multiplied, that would
			# be downpass * downpass, which 
			# results in focusing probability on the 
			# most-probable downpass state
			# Instead, the probabilities are just the downpass probabilities
			# The root node just gets the downpass probabilities
			ML_marginal_prob_each_state_at_branch_top_AT_node[anc_node_original_tree, ] = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[anc_node_original_tree, ]
			
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
				"This might occur in a highly constrained model, or if your data strongly contradicts your manual fixed\n",
				"likelihoods ('fixlikes') at some node(s) ('fixnode').\n",
				"As a 'fix', the downpass probabilities are being used for those nodes. But this is NOT RECOMMENDED!\n",
				"You should instead figure out what is causing the problem.", sep="")
				cat(txt)

				cat("\n\nPrinting (partial) downpass, uppass, and probability matrices to screen:\n\n", sep="") 
				print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS)				
				print(relative_probs_of_each_state_at_branch_top_AT_node_UPPASS)
				print(ML_marginal_prob_each_state_at_branch_top_AT_node)
				
				
				ML_marginal_prob_each_state_at_branch_top_AT_node[NaN_TF,] = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[NaN_TF,]
				}
			
			
			# Save them
			calc_loglike_sp_stratified_results$ML_marginal_prob_each_state_at_branch_bottom_below_node = ML_marginal_prob_each_state_at_branch_bottom_below_node
			calc_loglike_sp_stratified_results$ML_marginal_prob_each_state_at_branch_top_AT_node = ML_marginal_prob_each_state_at_branch_top_AT_node
	
			# tmpres$ML_marginal_prob_each_state_at_branch_top_AT_node
			# tmpres$ML_marginal_prob_each_state_at_branch_bottom_below_node
		}

	
		calc_loglike_sp_stratified_results$grand_total_likelihood = grand_total_likelihood
		
		# 2014-02-05_NJM fix
		calc_loglike_sp_stratified_results$total_loglikelihood = grand_total_likelihood
		
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
	} # END calc_loglike_sp_stratified



# Negative version of calc_loglike_for_optim_stratified() 
# (for e.g. minimization with GenSA)
calc_loglike_for_optim_stratified_neg <- function(params, BioGeoBEARS_run_object, phy, tip_condlikes_of_data_on_each_state, print_optim=TRUE, areas_list, states_list, force_sparse=FALSE, cluster_already_open=FALSE, min_branchlength=0.000001)
	{
	logLike = calc_loglike_for_optim_stratified(params, BioGeoBEARS_run_object, phy, tip_condlikes_of_data_on_each_state, print_optim=print_optim, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, min_branchlength=min_branchlength)
	
	neg_logLike = -1 * logLike
	return(neg_logLike)
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
#' @param min_branchlength Nodes with branches below this branchlength will not be treated as cladogenesis events; instead, they will be treated as 
#' if an OTU had been sampled from an anagenetic lineage, i.e. as if you had a direct ancestor.  This is useful for putting fossils into the biogeography analysis,
#' when you have fossil species that range through time. (Note: the proper way to obtain such trees, given that most phylogenetic methods force all OTUs to be tips 
#' rather than direct ancestors, is another question subject to active research.  However, one method might be to just set a branch-length cutoff, and treat any
#' branches sufficiently small as direct ancestors.)
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
calc_loglike_for_optim_stratified <- function(params, BioGeoBEARS_run_object, phy, tip_condlikes_of_data_on_each_state, print_optim=TRUE, areas_list, states_list, force_sparse=FALSE, cluster_already_open=FALSE, min_branchlength=0.000001)
	{
	# Put the parameters into the BioGeoBEARS_model_object, so that they can be universally read out
	# into any function
	BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
	
	# Put the parameters into the BioGeoBEARS_model_object, so that they can be universally read out
	# into any function
	BioGeoBEARS_model_object = params_into_BioGeoBEARS_model_object(BioGeoBEARS_model_object=BioGeoBEARS_model_object, params=params)
	
		
	######################################################
	# 2016-03-23_NJM: adding rescaling
	# (unscale params, if they were used before)
	######################################################
	if (BioGeoBEARS_run_object$rescale_params == TRUE)
		{
		#print("Before unscaling:")
		#print(BioGeoBEARS_model_object@params_table)
		BioGeoBEARS_model_object@params_table = unscale_BGB_params(scaled_params_table=BioGeoBEARS_model_object@params_table)
		
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table = BioGeoBEARS_model_object@params_table
		#print("After unscaling:")
		#print(BioGeoBEARS_model_object@params_table)
		}
	
	
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

	#dmat_times_d = matrix(d, nrow=length(areas), ncol=length(areas))
	#elist = rep(e, length(areas))
	
	# Set up the instantaneous rate matrix (Q matrix)
	#Qmat = rcpp_states_list_to_DEmat(areas_list=areas_list, states_list=states_list, dmat=dmat_times_d, elist=elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)

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
	# shorten the states_indices by 1 (cutting the 
	# null range state from the speciation matrix)
# 	if (include_null_range == TRUE)
# 		{
# 		states_indices[1] = NULL
# 		} # END if (include_null_range == TRUE)
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





	#######################################################
	# Get the detection model
	#######################################################
	if (BioGeoBEARS_run_object$use_detection_model == TRUE)
		{
		mean_frequency = BioGeoBEARS_model_object@params_table["mf","est"]
		dp = BioGeoBEARS_model_object@params_table["dp","est"]
		fdp = BioGeoBEARS_model_object@params_table["fdp","est"]

		# Calculate the initial tip likelihoods, using the detection model
		# Assumes correct order, double-check this
		numareas = length(areas)
		detects_df = BioGeoBEARS_run_object$detects_df
		controls_df = BioGeoBEARS_run_object$controls_df
		
		# return_LnLs=TRUE ensures no under-flow
		tip_condlikes_of_data_on_each_state = tiplikes_wDetectionModel(states_list_0based_index=states_list, phy=phy, numareas=numareas, detects_df=detects_df, controls_df=controls_df, mean_frequency=mean_frequency, dp=dp, fdp=fdp, null_range_gets_0_like=TRUE, return_LnLs=TRUE, relative_LnLs=TRUE, exp_LnLs=TRUE, error_check=TRUE)
		}





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
	include_null_range=BioGeoBEARS_run_object$include_null_range

	ttl_loglike = try(calc_loglike_sp_stratified(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=NULL, spPmat=NULL, min_branchlength=min_branchlength, return_what="loglike", probs_of_states_at_root=NULL, rootedge=TRUE, sparse=force_sparse, printlevel=printlevel, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=NULL, cppSpMethod=3, cluster_already_open=cluster_already_open, calc_ancprobs=FALSE, include_null_range=include_null_range, fixnode=fixnode, fixlikes=fixlikes, inputs=BioGeoBEARS_run_object, allareas=areas, all_states_list=states_list, return_condlikes_table=return_condlikes_table, calc_TTL_loglike_from_condlikes_table=calc_TTL_loglike_from_condlikes_table))

	if (("try-error" %in% class(ttl_loglike)) == TRUE)
		{
		ttl_loglike = BioGeoBEARS_run_object$on_NaN_error
		print_optim = FALSE
		}

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
	} # END calc_loglike_for_optim_stratified







