# Runs biogeographic stochastic mapping 
# (non-stratified, or stratified, versions, read from the 
#  input res object)
# 
# The 'stochastic_mapping_inputs_list' object comes from:
# 
# Non-stratified:
# 
# 
# stochastic_mapping_inputs_list2 = get_inputs_for_stochastic_mapping_stratified(res=resDIVALIKEj)

# wait_before_save = number of iterations of a meaningless for-loop, to 
# space-out the saving operations to avoid memory faults in R.app

# 
runBSM <- function(res, stochastic_mapping_inputs_list, maxnum_maps_to_try=1, nummaps_goal=maxnum_maps_to_try, maxtries_per_branch=40000, save_after_every_try=TRUE, savedir=getwd(), seedval=as.numeric(Sys.time()), wait_before_save=0.1, master_nodenum_toPrint=0)
	{
	defaults='
	maxnum_maps_to_try=1
	nummaps_goal=maxnum_maps_to_try
	maxtries_per_branch=40000
	save_after_every_try=TRUE
	savedir=getwd()
	seedval=12345
	'
	
	
	# Necessary default setting to avoid numbers instead of states
	options(stringsAsFactors=FALSE)
	
	# Starter
	clado_events_tables = list()
	ana_events_tables = list()
	
	# 2016-05-05_bug fix
	lnum = 0
	
	# Reconcile # of maps to try with nummaps goal
	if (nummaps_goal > maxnum_maps_to_try)
		{
		maxnum_maps_to_try = 2 * nummaps_goal
		} # END if (nummaps_goal > maxnum_maps_to_try)
	

	# Loop through the number of attempts
	for (m in 1:maxnum_maps_to_try)
		{
		# Initialize
		stochastic_mapping_results = NA
		clado_events_table = NA
		ana_events_table = NA
		
		# Skip loop once you get to nummaps_goal
		# (the desired number of successes)
		if (lnum >= nummaps_goal)
			{
			next()
			}
		
		startseed = seedval + m

		# Check if this is stratified or not
		strat_TF = FALSE
		if (is.null(stochastic_mapping_inputs_list$stratified) == FALSE)
			{
			if (stochastic_mapping_inputs_list$stratified == TRUE)
				{
				strat_TF = TRUE
				}
			} else {
			strat_TF = TRUE
			} # END if (is.null(stochastic_mapping_inputs_list$stratified) == FALSE)
		
		#print(paste0("runBSM startseed = ", startseed))
		
		# Run stratified or non-stratified stochastic mapping
		if (strat_TF == FALSE)
			{
			cmdstr = paste0("stochastic_mapping_results = stochastic_map_given_inputs(stochastic_mapping_inputs=stochastic_mapping_inputs_list, piecenum=NULL, maxtries=", maxtries_per_branch, ", seedval=", startseed, ", include_null_range=res$inputs$include_null_range, master_nodenum_toPrint=master_nodenum_toPrint)")
			} else {
			# Stratified
			#maxtries=40000; seedval=521; master_nodenum_toPrint=0
			cmdstr = paste0("stochastic_mapping_results = stochastic_mapping_on_stratified(res=res, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxtries=", maxtries_per_branch, ", seedval=", startseed, ", master_nodenum_toPrint=master_nodenum_toPrint)")
			} # END if (strat_TF == FALSE)
		#print("AA_stochastic_mapping_on_stratified")
		#print(cmdstr)
		
		# Run either of the options
		#cat(" ", m, sep="")
		#print("try_result:")
		try_result = try(expr=eval(parse(text=cmdstr)))
		#eval(parse(text=cmdstr))
		#print(try_result)
		#cat("\n")
		#print(warnings())
		#print(class(try_result))
		#cat("\n")
		#cat(",", sep="")
		
		# 2016-05-05
		# Now that stochastic_mapping_results is initialized to NA each loop, should
		# have no weird blank results
		TF2 = ( (length(stochastic_mapping_results)==1) && (is.na(stochastic_mapping_results) == TRUE) )
		
		if ( (class(try_result) != "try-error") && TF2==FALSE )
			{
			if (strat_TF == FALSE)
				{
				##########################################
				# For NON-stratified stochastic mapping
				##########################################
				clado_events_table = stochastic_mapping_results
				
				#print(clado_events_table[20,])
				
				# Error check, for no anagenetic events
				TF1 = stochastic_mapping_results$anagenetic_events_txt_below_node == "none"
				TF2 = is.na(stochastic_mapping_results$anagenetic_events_txt_below_node)
				TF3 = stochastic_mapping_results$anagenetic_events_txt_below_node == "NA"
				TF4 = stochastic_mapping_results$anagenetic_events_txt_below_node == ""
				TFall = (TF1 + TF2 + TF3 + TF4) 
				TF = TFall == 0
				
				if (sum(TF, na.rm=TRUE) > 0)
					{
					# Some anagenetic events, make table
					#print("here1")
					ana_events_table = events_txt_list_into_events_table(events_txt_list=clado_events_table$anagenetic_events_txt_below_node, trtable=clado_events_table, recalc_abs_ages=TRUE)
					#print("here2")
					} else {
					# No anagenetic events
					ana_events_table = NA
					}
				} else {
				####################################daughter_nds
				# For STRATIFIED stochastic mapping
				####################################
				#clado_colnames = c("node", "node.type", "edge.length", "time_bp", "stratum", "piecenum", "piececlass", "SUBnode", "SUBnode.type", "SUBedge.length", "clado_event_type", "clado_event_txt", "clado_dispersal_to", "sampled_states_AT_nodes", "sampled_states_AT_brbots", "left_desc_nodes", "right_desc_nodes", "samp_LEFT_dcorner", "samp_RIGHT_dcorner", "anagenetic_events_txt_below_node")
				clado_colnames = c("node", "node.type", "parent_br", "edge.length", "ancestor", "daughter_nds", "time_bp", "fossils", "label", "stratum", "time_top", "time_bot", "reltimept", "piecenum", "piececlass", "SUBnode", "SUBnode.type", "SUBparent_br", "SUBedge.length", "SUBancestor", "SUBdaughter_nds" , "SUBtime_bp", "SUBfossils", "SUBlabel", "clado_event_type", "clado_event_txt", "clado_dispersal_to", "sampled_states_AT_nodes", "sampled_states_AT_brbots", "left_desc_nodes", "right_desc_nodes", "samp_LEFT_dcorner", "samp_RIGHT_dcorner", "anagenetic_events_txt_below_node")
				#rowsTF = stochastic_mapping_results$master_table_cladogenetic_events$node.type != "tip"
				clado_events_table = stochastic_mapping_results$master_table_cladogenetic_events[, clado_colnames]
				ana_events_table = stochastic_mapping_results$table_w_anagenetic_events
				} # END if (strat_TF == FALSE)
			
			# Print update
			lnum = lnum+1
			cat("\nSuccess #", lnum, "/", nummaps_goal, " on stochastic mapping attempt #", m, "/", maxnum_maps_to_try, " tries.", sep="")
			
			# Wait a second or two, every 100 successes
			# -- attempting to avoid segmentation faults during
			# rapid stochastic mapping
			if (FALSE)
				{
				Sys.sleep(wait_before_save)
				}
			
			if (TRUE)
				{
				if (lnum/50 == round(lnum/50))
					{
					Sys.sleep(wait_before_save * 3)
					} # END if (lnum/100 == round(lnum/100))
				}
			
			# Store results
			clado_events_tables[[lnum]] = clado_events_table
			# Error trap -- if no anagenetic events, put NA
			# 2016-05-05 bug: is.na(ana_events_table) is flagging yes when whole 
			# row is NAs
			if ((length(ana_events_table) > 1) && (is.na(ana_events_table)==TRUE))
				{
				error_msg = "STOP ERROR in runBSM(): your anagenetic events table seems to have row(s) that have multiple NAs. This could be due to an out of date 'BioGeoBEARS_stochastic_mapping_v1.R' file, missing this bug fix: '# NJM 2016-05-05 bug fix: add 'as.numeric'   rownums_in_trtable = as.numeric(tmptable_rows$nodenum_at_top_of_branch)'"
				cat("\n\n")
				cat(error_msg)
				cat("\n\n")
				stop(error_msg)
				}
			
			
			if ( is.null(ana_events_table) || (length(ana_events_table)==1 && is.na(ana_events_table)) || (length(ana_events_table) < 1) )
				{
				cat("\nWARNING: 'ana_events_table' was null, NA or length <1. Replacing with 'NA'\n")
				ana_events_table = NA
				}
			ana_events_tables[[lnum]] = ana_events_table
			
			if (save_after_every_try == TRUE)
				{
				# Sometimes it looks like saving happens too fast and R crashes
				if (is.null(wait_before_save) == FALSE)
					{
					# Wait this long before saving
					Sys.sleep(wait_before_save)
					} # END if (is.null(wait_before_save) == FALSE)
				
				# Archive every round
				RES_clado_events_tables = clado_events_tables
				RES_ana_events_tables = ana_events_tables

				RES_clado_events_tables_fn = slashslash(paste0(savedir, "/", "RES_clado_events_tables_PARTIAL.Rdata"))
				RES_ana_events_tables_fn = slashslash(paste0(savedir, "/", "RES_ana_events_tables_PARTIAL.Rdata"))
				
				# Message about saving
				cat("...saving to RES_clado_events_tables_PARTIAL.Rdata")
				cmdstr = "save(RES_clado_events_tables, file=RES_clado_events_tables_fn)"
				savetry = try(expr=cmdstr)
				if (class(savetry) == "try-error")
					{
					cat("...save FAILED for 'RES_clado_events_tables_PARTIAL.Rdata' for some reason, sometimes rapid saves seem to crash on some systems.\n")
					} else {
					cat("...saved 'RES_clado_events_tables_PARTIAL.Rdata'\n")
					}
				#print(RES_ana_events_tables)
				
				cmdstr = "save(RES_ana_events_tables, file=RES_ana_events_tables_fn)"
				savetry = try(expr=eval(parse(text=cmdstr)))
				if (class(savetry) == "try-error")
					{
					cat("...save FAILED for 'RES_ana_events_tables_PARTIAL.Rdata' for some reason, sometimes rapid saves seem to crash on some systems.\n")
					} else {
					cat("...saved 'RES_ana_events_tables_PARTIAL.Rdata'")
					}
				} # END if (save_after_every_try == TRUE)
			
			} else {
			cat("\nFailure on stochastic mapping attempt #", m, "/", maxnum_maps_to_try, " tries. Holding at success #", lnum, "/", nummaps_goal, ".", sep="")
			cat("\n")
			print(try_result)
			} # END if (class(try_result) != "try-error")

		} # END for (m in 1:maxnum_maps_to_try)
	
	# Archive at the end
	RES_clado_events_tables = clado_events_tables
	RES_ana_events_tables = ana_events_tables
	RES_clado_events_tables_fn = slashslash(paste0(savedir, "/", "RES_clado_events_tables.Rdata"))
	RES_ana_events_tables_fn = slashslash(paste0(savedir, "/", "RES_ana_events_tables.Rdata"))
	
	# Print results
	cat("\n\n\n=============================================\nBiogeographic Stochastic Mapping is finished.\n\n")
	txt = paste0("Saving cladogenetic and anagenetic events tables. ")
	cat(txt)
	txt2 = paste0("The commands:\n\nload(file=", RES_clado_events_tables_fn, ")\n\nload(file=", RES_ana_events_tables_fn, ")\n\n...will load them to objects RES_clado_events_tables and RES_ana_events_tables.")
	cat(txt2)
	cat("\n\n")
		
	save(RES_clado_events_tables, file=RES_clado_events_tables_fn)
	save(RES_ana_events_tables, file=RES_ana_events_tables_fn)
	
	cat("Now returning to console 'BSM_output', which contains:\n")
	cat("BSM_output$RES_clado_events_tables")
	cat("\n")
	cat("BSM_output$RES_ana_events_tables")
	cat("\n=============================================")
	cat("\n")
	cat("\n")
	
	BSM_output = NULL
	BSM_output$RES_clado_events_tables = RES_clado_events_tables
	BSM_output$RES_ana_events_tables = RES_ana_events_tables
	
	extract='
	clado_events_tables = BSM_output$RES_clado_events_tables
	ana_events_tables = BSM_output$RES_ana_events_tables
	'
	
	return(BSM_output)
	} # END runBSM()


