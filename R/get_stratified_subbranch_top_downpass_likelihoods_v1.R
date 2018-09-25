get_tip_likelihoods_of_subbranch_from_resCondlikes_given_master_table <- function(stratum, piecenum, master_table, condlikes, relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE)
	{
	# The output of calc_loglike_sp_strat:
	
	# $relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_table -- the relative probabilities at the BOTTOMS
	# of branches -- subbranch pieces will have NAs
	
	# $condlikes_table -- the regular likelihoods at the tips and nodes of subtrees
	# $condlikes_table -- for subbranches, the likelihoods at the BOTTOM of the subbranch, i.e. at the stratum boundary
	# $condlikes_table -- the tip likelihoods of the master tree are in stratum 0
	
	
	
	# Convert variable to match stochastic mapping setup
	timeperiod_i = stratum
	timeperiod_i_up = timeperiod_i - 1


	# Get the list of tree pieces for the stratum
	#stochastic_mapping_inputs = stochastic_mapping_inputs_list[[timeperiod_i]]
	#tree_sections = stochastic_mapping_tree_sections_list
	#tree_pieces = tree_sections$return_pieces_list
	#chainsaw_result = tree_pieces
	#num_pieces = length(tree_pieces)

	# Get the specific current treepiece (a subbranch)
	jj = piecenum
	#treepiece = tree_pieces[[jj]]
	
	# Find the row in the master_table
	TF1 = master_table$stratum == timeperiod_i
	TF2 = master_table$piecenum == jj
	TF3 = master_table$piececlass == "subbranch"
	TF = (TF1 + TF2 + TF3) == 3
	this_row_of_master_table_is_being_used = TF
	
	# Find the row
	#rownum_master_table = (1:nrow(condlikes_table))[TF]
	rownum_master_table = (1:nrow(master_table))[TF]
	tmp_master_table_row = master_table[this_row_of_master_table_is_being_used, ]
	
	# Error check -- you should get only one row for the current treepiece
	if (nrow(tmp_master_table_row) != 1)
		{
		errortxt = paste("\n\nSTOP ERROR in\this_row_of_master_table_is_the_node_in_the_next_stratum_up():\n\nnrow(tmp_master_table_row) != 1\n\n", sep="")
		cat(errortxt)
		print("tmp_master_table_row:")
		cat("\n\n")
		print(tmp_master_table_row)
		stop("Stopping on error")
		}

	# Fossil check
	# Now check if it's a fossil that appears in this time bin
	master_tip_time_bp = tmp_master_table_row$time_bp
	time_top = tmp_master_table_row$time_top
	time_bot = tmp_master_table_row$time_bot
	is_fossil = tmp_master_table_row$fossils

	# If this is TRUE, there's a match and the fossil tip appears in this time period
	if ( (master_tip_time_bp >= time_top) && (master_tip_time_bp < time_bot) && (is_fossil == TRUE))
		{
		# Shorten the branchlength by master_tip_time_bp-time_top
		#amount_to_shorten_by = master_tip_time_bp-time_top
		#subbranch_length = subbranch_length - amount_to_shorten_by
		#do_exponentiation = TRUE
		
		# If it's a fossil tip, it's tip likelihoods appear in the orig_tip part of the master_table
		node_at_top_of_subbranch_in_master_tree = tmp_master_table_row$node
	
		# Find the row in the master_table
		#TF1 = master_table$stratum == 0		# 0, because it's a tip likelihood in the master_tree table
		TF1 = master_table$stratum == timeperiod_i		# 2017-04-06: NO, not zero, use stratum
		
		# 2017-04-06 bug fix: added TF2 and new TF3 (original tips downpass were 0s; remove other fix if this works)
		TF2 = master_table$piecenum == jj
		#TF3 = master_table$piececlass == "orig_tip" # This forces us to pick from the orig_tips
		TF3 = master_table$SUBnode.type == "tip" # This forces us NOT to pick from the orig_tips
		TF4 = master_table$node == node_at_top_of_subbranch_in_master_tree
		TF = (TF1 + TF2 + TF3 + TF4) == 4
		this_row_of_master_table_is_the_node_of_the_master_tree_tip_fossil = (1:nrow(master_table))[TF]

		# Error check -- you should get only one row for the current treepiece
		if (length(this_row_of_master_table_is_the_node_of_the_master_tree_tip_fossil) != 1)
			{
			errortxt = paste("\n\nSTOP ERROR in\this_row_of_master_table_is_the_node_in_the_next_stratum_up():\n\nnrow(this_row_of_master_table_is_the_node_of_the_master_tree_tip_fossil) != 1\n\n", sep="")
			cat(errortxt)
			print("print(this_row_of_master_table_is_the_node_of_the_master_tree_tip_fossil):")
			print(this_row_of_master_table_is_the_node_of_the_master_tree_tip_fossil)

			print("TF1 = master_table$stratum == timeperiod_i ... print(timeperiod_i):")
			print(timeperiod_i)
			print("sum(TF1):")
			print(sum(TF1))
			print("TF2 = master_table$piecenum == jj ... print(jj / piecenum):")
			print(jj)
			print("sum(TF2):")
			print(sum(TF2))
			print("TF3 = master_table$SUBnode.type == 'tip' ... print(jj / piecenum):")
			print("sum(TF3):")
			print(sum(TF3))
			print("TF4 = master_table$node == node_at_top_of_subbranch_in_master_tree ... print(node_at_top_of_subbranch_in_master_tree):")
			print(node_at_top_of_subbranch_in_master_tree)
			print("sum(TF4):")
			print(sum(TF4))
			stop("Stopping on error")
			}
		
		# These should be the correct downpass conditional likelihoods at the top of the input piecenum (subbranch)
		condlikes_at_tip_in_master_tree = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE[this_row_of_master_table_is_the_node_of_the_master_tree_tip_fossil,]
		
		
		if (sum(condlikes_at_tip_in_master_tree) == 0)
			{
			txt = paste0("STOP ERROR_5A in get_tip_likelihoods_of_subbranch_from_resCondlikes_given_master_table():\nAt timeperiod_i=", timeperiod_i, ", sum(condlikes_at_tip_in_master_tree)==0, this_row_of_master_table_is_the_node_of_the_master_tree_tip_fossil=", this_row_of_master_table_is_the_node_of_the_master_tree_tip_fossil)
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			stop(txt)
			}		
		
		return(condlikes_at_tip_in_master_tree)
		}


	
	if (timeperiod_i > 1)
		{
		# Find the node in the stratum above (whether or not it's an actual node, or a sub-branch)
		# It's identified by just "node", since "node" is the node at the top of the branches in the master_tree
		node_at_top_of_subbranch_in_master_tree = tmp_master_table_row$node
	
		# Find the row in the master_table
		TF1 = master_table$stratum == timeperiod_i_up
		#TF2 = master_table$piecenum == jj
		TF3 = master_table$piececlass != "orig_tip" # This excludes the orig_tips
		TF4 = master_table$node == node_at_top_of_subbranch_in_master_tree
		TF = (TF1 + TF3 + TF4) == 3
		this_row_of_master_table_is_the_node_in_the_next_stratum_up = (1:nrow(master_table))[TF]

		# Error check -- you should get only one row for the current treepiece
		if (length(this_row_of_master_table_is_the_node_in_the_next_stratum_up) != 1)
			{
			errortxt = paste("\n\nSTOP ERROR in\this_row_of_master_table_is_the_node_in_the_next_stratum_up():\n\nnrow(this_row_of_master_table_is_the_node_in_the_next_stratum_up) != 1\n\n", sep="")
			cat(errortxt)
			print("this_row_of_master_table_is_the_node_in_the_next_stratum_up:")
			cat("\n\n")
			print(this_row_of_master_table_is_the_node_in_the_next_stratum_up)
			stop("Stopping on error")
			}

		# These should be the correct downpass conditional likelihoods at the top of the input piecenum (subbranch)
		condlikes_at_bottom_of_branch_below_node_in_stratum_above = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE[this_row_of_master_table_is_the_node_in_the_next_stratum_up,]
		
		if (sum(condlikes_at_bottom_of_branch_below_node_in_stratum_above) == 0)
			{
			txt = paste0("STOP ERROR_5B in get_tip_likelihoods_of_subbranch_from_resCondlikes_given_master_table():\nAt timeperiod_i=", timeperiod_i, ", sum(condlikes_at_bottom_of_branch_below_node_in_stratum_above)==0, this_row_of_master_table_is_the_node_in_the_next_stratum_up=", this_row_of_master_table_is_the_node_in_the_next_stratum_up)
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			stop(txt)
			}
		
		return(condlikes_at_bottom_of_branch_below_node_in_stratum_above)
		} else {
		# For getting the tip likelihoods at the tips of the master tree, when stratum == 0
		
		# Find the node in the stratum above (whether or not it's an actual node, or a sub-branch)
		# It's identified by just "node", since "node" is the node at the top of the branches in the master_tree
		node_at_top_of_subbranch_in_master_tree = tmp_master_table_row$node
	
		# Find the row in the master_table
		TF1 = master_table$stratum == timeperiod_i_up
		#TF2 = master_table$piecenum == jj
		TF3 = master_table$piececlass == "orig_tip" # This excludes the orig_tips
		TF4 = master_table$node == node_at_top_of_subbranch_in_master_tree
		TF = (TF1 + TF3 + TF4) == 3
		TF = (TF1 + TF3 + TF4) == 3
		this_row_of_master_table_is_the_node_of_the_master_tree_tip = (1:nrow(master_table))[TF]

		# Error check -- you should get only one row for the current treepiece
		if (length(this_row_of_master_table_is_the_node_of_the_master_tree_tip) != 1)
			{
			errortxt = paste("\n\nSTOP ERROR in\this_row_of_master_table_is_the_node_in_the_next_stratum_up():\n\nnrow(this_row_of_master_table_is_the_node_of_the_master_tree_tip) != 1\n\n", sep="")
			cat(errortxt)
			print("this_row_of_master_table_is_the_node_of_the_master_tree_tip:")
			cat("\n\n")
			print(this_row_of_master_table_is_the_node_of_the_master_tree_tip)
			stop("Stopping on error")
			}
		
		# These should be the correct downpass conditional likelihoods at the top of the input piecenum (subbranch)
		condlikes_at_tip_in_master_tree = condlikes[this_row_of_master_table_is_the_node_of_the_master_tree_tip,]

		if (sum(condlikes_at_tip_in_master_tree) == 0)
			{
			txt = paste0("STOP ERROR_5C in get_tip_likelihoods_of_subbranch_from_resCondlikes_given_master_table():\nAt timeperiod_i=", timeperiod_i, ", sum(condlikes_at_tip_in_master_tree)==0. this_row_of_master_table_is_the_node_of_the_master_tree_tip=", this_row_of_master_table_is_the_node_of_the_master_tree_tip)
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			stop(txt)
			}

		return(condlikes_at_tip_in_master_tree)
		}
	
	stop("Stopping: you shouldn't get here.")
	return(NULL)
	}