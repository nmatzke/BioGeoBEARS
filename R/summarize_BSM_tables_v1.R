#######################################################
# Intrinsically, the transition matrix used in BioGeoBEARS, Lagrange, etc.,
# does not specify the exact source area when a transition occurs
# from a widespread ancestor. However, this can be simulated given
# knowledge of dmat (dmat = all dispersal multipliers multiplied)
#######################################################

simulate_source_areas_ana_clado <- function(res, clado_events_tables, ana_events_tables, areanames)
	{
	# Get the dmat and times
	dmat_times = get_dmat_times_from_res(res=res, numstates=NULL)
	dmat = dmat_times$dmat
	times = dmat_times$times
	
	# Error check
	if (length(clado_events_tables) != length(ana_events_tables))
		{
		txt = paste0("STOP ERROR in simulate_source_areas_ana_clado(): ana_events_tables and clado_events_tables should have the same length. Instead, length(ana_events_tables)=", length(ana_events_tables), ", length(clado_events_tables)=", (clado_events_tables), ".")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		} # END if (length(clado_events_tables) != length(clado_events_tables))
	
	# Simulate 
	# - the anagenetic dispersal event source areas
	# - the cladogenetic dispersal event (jump dispersal) source areas
	numBSMs = length(clado_events_tables)
	txt = paste0("Simulating the ancestral areas for cases where the ancestral area occupied 2 or more areas (number of stochastic maps=", numBSMs, "):\n")
	cat("\n\n")
	cat(txt)
	for (i in 1:numBSMs)
		{
		cat(i, " ", sep="")
		
		# Simulate the source areas for the jump-dispersal events;
		# modify clado_event_tables
		clado_events_tables[[i]] = simulate_source_area_clado(clado_events_tables[[i]], areanames, dmat=dmat, times=times)

		# Simulate the source areas for the anagenetic-dispersal events;
		# modify ana_event_tables
		ana_events_tables[[i]] = simulate_source_area_ana(ana_events_tables[[i]], areanames, dmat=dmat, times=times)
		} # END for (i in 1:numBSMs)
	cat("\n\n")	
		
	BSMs_w_sourceAreas = NULL
	BSMs_w_sourceAreas$clado_events_tables = clado_events_tables
	BSMs_w_sourceAreas$ana_events_tables = ana_events_tables
	
	extract='
	clado_events_tables = BSMs_w_sourceAreas$clado_events_tables
	ana_events_tables = BSMs_w_sourceAreas$ana_events_tables
	'
	
	return(BSMs_w_sourceAreas)
	}



# areanames = names(tipranges@df)
# actual_names = areanames
# actual_names
# 
# # Get the dmat and times (if any)
# dmat_times = get_dmat_times_from_res(res=res, numstates=NULL)
# dmat_times
# 
# # Simulate the source areas for the jump-dispersal events; add to table
# clado_events_table = clado_events_tables[[1]]
# names(clado_events_table)
# clado_events_table = simulate_source_area_clado(clado_events_table, areanames, dmat=dmat_times$dmat)
# names(clado_events_table)
# 
# # Simulate the source areas for the jump-dispersal events; add to table
# ana_events_table = ana_events_tables[[1]]
# names(ana_events_table)
# ana_events_table = simulate_source_area_ana(ana_events_table, areanames, dmat=NULL)
# names(ana_events_table)


simulate_source_area_ana <- function(ana_events_table, areanames, dmat=NULL, times=NULL)
	{
	# NA if no events; return NULL
	if ( (length(ana_events_table)==1) && (is.na(ana_events_table) == TRUE) )
		{
		return(ana_events_table)
		} # END if (is.na(ana_events_table) == TRUE)
	
	ana_events_table_OLD = ana_events_table
	
	# Check for times (different dmats through time)
	strat_TF = FALSE
	if (is.list(dmat) && is.null(times))
		{
		txt = "STOP ERROR in simulate_source_area_ana(): if 'dmat' is a list of dmats, 'times' cannot be NULL."
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		}
	if (is.list(dmat) && !is.null(times))
		{
		strat_TF = TRUE
		
		# Get the top/bottom of the time bins
		tops = c(0, times[-length(times)])
		bots = times
		
		if (length(dmat) != length(times))
			{
			txt = "STOP ERROR in simulate_source_area_ana(): the length of the 'dmat' list must equal the length of 'times' in a time-stratified analysis."
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			stop(txt)
			} # END if (length(dmat) != length(times))
		
		if ( all(times == sort(times))==FALSE || (0 %in% times) == TRUE )
			{
			txt = "STOP ERROR in simulate_source_area_ana(): 'times' in a time-stratified analysis must be sorted from youngest to oldest. Also, time '0' should *not* be included."
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			stop(txt)
			} # END if ( all( times == sort(times) ) )
		} # END if (is.list(dmat) && is.null(times))
	
	# If dmat is NULL, set to all 1s
	# This is the dispersal_multipliers_matrix, which includes the influence 
	# of manual dispersal multipliers, distance, environmental distance,
	# etc. on the dispersal (d) and jump dispersal (j) processes.
	# (and anagenetic range-switching, a)
	if (is.null(dmat) == TRUE)
		{
		dmat = matrix(1, nrow=length(areanames), ncol=length(areanames))
		} # END if (is.null(dmat) == TRUE)

	# Actual source area (single area)
	ana_dispersal_from = rep("", nrow(ana_events_table))
	for (ii in 1:length(areanames))
		{
		areaname = areanames[ii]
		TF = ana_events_table$dispersal_to == areaname
		TF[is.na(TF)] = FALSE
	
		if ((length(TF) > 0) && (sum(TF) > 0))
			{
			# Split events on "->"
			tmp_table = ana_events_table[TF,]
			events_txt = tmp_table$event_txt

			events_df = as.data.frame(matrix(unlist(sapply(X=events_txt, FUN=strsplit, split="->")), ncol=2, byrow=TRUE), stringsAsFactors=FALSE)
			names(events_df) = c("ancrange", "decrange")
			orig_rownum = 1:nrow(events_df)

			# OLD: nchars is labeled as the weight
			#nchars = nchar(events_df$ancrange)

			# NEW: weight by prob of being the source area,
			# then pick a source area
			tmp_ana_dispersal_from = rep(NA, nrow(tmp_table))
			for (j in 1:nrow(tmp_table))
				{
				# Destination area was set by i through areanames
				# (assumes they are in the same order); filter
				# dmat accordingly
				if (strat_TF == FALSE)
					{
					source_area_relprob1 = dmat[,ii]
					} else {
					timebp_of_event = tmp_table$abs_event_time[j]
					TF1 = timebp_of_event >= tops
					TF2 = timebp_of_event < bots
					inbin_TF = (TF1 + TF2) == 2
					binnum = (1:length(times))[inbin_TF]
					
					source_area_relprob1 = dmat[[binnum]][,ii]
					} # END if (strat_TF == FALSE)

				# source range
				source_range = events_df$ancrange[j]
				source_areas = strsplit(x=source_range, split="")[[1]]
				source_indices = match(x=source_areas, table=areanames)
				source_area_relprob2 = source_area_relprob1[source_indices]
				
				# Check for zeros (might happen in "manual" histories)
				if (sum(source_area_relprob2) <= 0)
					{
					stoptxt = paste0("simulate_source_area_ana() says: WARNING: the sum of source_area_relprob2 is 0. This means that all possible source areas have probability zero. In stochastic mapping, this could happen if, e.g. in a time-stratified stochastic map, the simulated state at the top of the previous time-bin is an impossible starting range for the manual history in the current time bin. In this situation, the 'manual history' option probably forced a history, just to stop you from running millions of tries pointlessly. This is part of the overall problem with very-improbable histories in stochastic mapping with huge biogeographical state spaces and complex modifications and lots of constratints.  If you get this warning, you should (a) be grateful that this code works at all, programming all of this is ridiculously complex, especially in R, and (b) you should probably delete these histories from your collection, if you can, although the overall effect should be minor.  Here, we are re-setting source_area_relprob2 to source_area_relprob2=rep(1, times=length(source_area_relprob2)). The events_txt being tried here was: '", events_txt, "'.")
					cat("\n\n")
					cat(stoptxt)
					cat("\n\n")
					warning(stoptxt)
					
					source_area_relprob2=rep(1, times=length(source_area_relprob2))
					}
				
				source_area_relprob = source_area_relprob2 / sum(source_area_relprob2)
				
				# Randomly sample a source area
				source_area = sample(x=source_areas, size=1, replace=FALSE, prob=source_area_relprob)
				tmp_ana_dispersal_from[j] = source_area
				} # END for (j in 1:nrow(tmp_table))
			# Save in the overall ana_dispersal_from
			ana_dispersal_from[TF] = tmp_ana_dispersal_from
			} # END if ((length(TF) > 0) && (sum(TF) > 0))
		} # END for (ii in 1:length(areanames))

	# Edit the input ana_events_table
	dispersal_to_colnum_TF = names(ana_events_table_OLD) == "dispersal_to"
	dispersal_to_colnum = (1:ncol(ana_events_table_OLD))[dispersal_to_colnum_TF]
	
	end1 = (dispersal_to_colnum-1)
	start2 = (dispersal_to_colnum)
	end2 = ncol(ana_events_table_OLD)
	if ( (dispersal_to_colnum+1) <= ncol(ana_events_table_OLD))
		{
		ana_events_table = cbind(ana_events_table_OLD[1:end1], ana_dispersal_from, ana_events_table_OLD[start2:end2])
		} else {
		ana_events_table = cbind(ana_events_table_OLD[1:end1], ana_dispersal_from)
		} # END if ( (dispersal_to_colnum+1) <= ncol(ana_events_table_OLD))
	
	return(ana_events_table)
	} # END simulate_source_area_ana <- function(ana_events_tables, ana_events_tables, areanames, dmat=NULL)



simulate_source_area_clado <- function(clado_events_table, areanames, dmat=NULL, times=NULL)
	{
	clado_events_table_OLD = clado_events_table

	# Check for times (different dmats through time)
	strat_TF = FALSE
	if (is.list(dmat) && is.null(times))
		{
		txt = "STOP ERROR in simulate_source_area_ana(): if 'dmat' is a list of dmats, 'times' cannot be NULL."
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		}
	if (is.list(dmat) && !is.null(times))
		{
		strat_TF = TRUE
		
		# Get the top/bottom of the time bins
		tops = c(0, times[-length(times)])
		bots = times
		
		if (length(dmat) != length(times))
			{
			txt = "STOP ERROR in simulate_source_area_ana(): the length of the 'dmat' list must equal the length of 'times' in a time-stratified analysis."
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			stop(txt)
			} # END if (length(dmat) != length(times))
		
		if ( all(times == sort(times))==FALSE || (0 %in% times) == TRUE )
			{
			txt = "STOP ERROR in simulate_source_area_ana(): 'times' in a time-stratified analysis must be sorted from youngest to oldest. Also, time '0' should *not* be included."
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			stop(txt)
			} # END if ( all( times == sort(times) ) )
		} # END if (is.list(dmat) && is.null(times))
	
	
	# If dmat is NULL, set to all 1s
	# This is the dispersal_multipliers_matrix, which includes the influence 
	# of manual dispersal multipliers, distance, environmental distance,
	# etc. on the dispersal (d) and jump dispersal (j) processes.
	# (and anagenetic range-switching, a)
	if (is.null(dmat) == TRUE)
		{
		dmat = matrix(1, nrow=length(areanames), ncol=length(areanames))
		} # END if (is.null(dmat) == TRUE)

	# Actual source area (single area)
	clado_dispersal_from = rep("", nrow(clado_events_table))
	for (i in 1:length(areanames))
		{
		areaname = areanames[i]
		TF = clado_events_table$clado_dispersal_to == areaname
		TF[is.na(TF)] = FALSE
		
		if ((length(TF) > 0) && (sum(TF) > 0))
			{
			# Split events on "->"
			tmp_table = clado_events_table[TF,]
			events_txt = tmp_table$clado_event_txt
	
			events_df = as.data.frame(matrix(unlist(sapply(X=events_txt, FUN=strsplit, split="->")), ncol=2, byrow=TRUE), stringsAsFactors=FALSE)
			names(events_df) = c("ancrange", "decrange")
			orig_rownum = 1:nrow(events_df)

			# OLD: nchars is labeled as the weight
			#nchars = nchar(events_df$ancrange)

			# NEW: weight by prob of being the source area,
			# then pick a source area
			tmp_clado_dispersal_from = rep(NA, nrow(tmp_table))
			for (j in 1:nrow(tmp_table))
				{

				# Destination area was set by i through areanames
				# (assumes they are in the same order); filter
				# dmat accordingly
				if (strat_TF == FALSE)
					{
					source_area_relprob1 = dmat[,i]
					} else {
					timebp_of_event = tmp_table$time_bp[j]
					TF1 = timebp_of_event >= tops
					TF2 = timebp_of_event < bots
					inbin_TF = (TF1 + TF2) == 2
					binnum = (1:length(times))[inbin_TF]
					source_area_relprob1 = dmat[[binnum]][,i]
					} # END if (strat_TF == FALSE)

				# source range
				source_range = events_df$ancrange[j]
				source_areas = strsplit(x=source_range, split="")[[1]]
				source_indices = match(x=source_areas, table=areanames)
				source_area_relprob2 = source_area_relprob1[source_indices]
				source_area_relprob = source_area_relprob2 / sum(source_area_relprob2)

				# Randomly sample a source area
				source_area = sample(x=source_areas, size=1, replace=FALSE, prob=source_area_relprob)
				tmp_clado_dispersal_from[j] = source_area
				} # END for (j in 1:nrow(tmp_table))
			# Save in the overall clado_dispersal_from
			clado_dispersal_from[TF] = tmp_clado_dispersal_from
			} # END if ((length(TF) > 0) && (sum(TF) > 0))
		} # END for (i in 1:length(areanames))
	
	# Edit the input clado_events_table
	clado_to_colnum_TF = names(clado_events_table_OLD) == "clado_dispersal_to"
	clado_to_colnum = (1:ncol(clado_events_table_OLD))[clado_to_colnum_TF]
	
	end1 = (clado_to_colnum-1)
	start2 = (clado_to_colnum)
	end2 = ncol(clado_events_table_OLD)
	if ( (clado_to_colnum+1) <= ncol(clado_events_table_OLD))
		{
		clado_events_table = cbind(clado_events_table_OLD[1:end1], clado_dispersal_from, clado_events_table_OLD[start2:end2])
		} else {
		clado_events_table = cbind(clado_events_table_OLD[1:end1], clado_dispersal_from)
		} # END if ( (clado_to_colnum+1) <= ncol(clado_events_table_OLD))
	
	return(clado_events_table)
	} # END simulate_source_area_clado <- function(clado_events_tables, ana_events_tables, areanames, dmat=NULL)






# Assumes equal weights of ancestral areas, when ancestor
# occupies more than one area
count_clado_dispersal_events <- function(clado_events_table, areanames, actual_names=areanames)
	{
	defaults='
	areanames = c("A", "B", "C", "D", "E", "F", "G")
	actual_names = c("SAm", "CAm", "Car", "NAm", "AF", "EU", "OZ")
	clado_events_table = DECJ_clado_events_tables[[1]]
	'

	# Error check
	if ("clado_dispersal_from" %in% names(clado_events_table) == FALSE)
		{
		txt = paste0("STOP ERROR in count_clado_dispersal_events(): input 'clado_events_table' must have a column 'clado_dispersal_from'. To add this column, run simulate_source_area_clado() or simulate_source_areas_ana_clado().")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		} # END error check

	
	TF = nchar(areanames) > 1
	if (any(TF))
		{
		txt = paste0("STOP ERROR: all 'areanames' must be only 1 character long for the counting of events in stochastic mapping to work. You have violations for these 'areanames': ", paste(areanames[TF], sep=", "), ".")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		} # END if (any(TF))


	
	# Built up the transition matrix
	events_df2 = NULL
	for (i in 1:length(areanames))
		{
		areaname = areanames[i]
		TF = clado_events_table$clado_dispersal_to == areaname
		TF[is.na(TF)] = FALSE
		tmp_table = clado_events_table[TF, ]
		
		if ((length(TF) > 0) && (sum(TF) > 0))
			{
			# Split them on "->"
			events_txt = tmp_table$clado_event_txt
	
			events_df = as.data.frame(matrix(unlist(sapply(X=events_txt, FUN=strsplit, split="->")), ncol=2, byrow=TRUE), stringsAsFactors=FALSE)
			names(events_df) = c("ancrange", "decrange")
			orig_rownum = 1:nrow(events_df)
			
			time_bp = tmp_table$time_bp
			ancrange = events_df$ancrange
			decrange = events_df$decrange
			clado_dispersal_from = tmp_table$clado_dispersal_from
			clado_dispersal_to = tmp_table$clado_dispersal_to
			
			events_df = cbind(clado_events_table$time_bp[TF], events_df, clado_dispersal_from, clado_events_table$clado_dispersal_to[TF])
			names(events_df) = c("time_bp", "ancrange", "decrange", "clado_dispersal_from", "clado_dispersal_to")
			
			events_df2 = rbind(events_df2, events_df)
			names(events_df2) = c("time_bp", "ancrange", "decrange", "clado_dispersal_from", "clado_dispersal_to")
			} # END if (sum(TF) > 0)
		} # # END for (i in 1:areanames)

	

	
	if (is.null(events_df2) == TRUE)
		{
		counts_df = NULL
		return(counts_df)
		} # END if (is.null(events_df2) == TRUE)

	counts_matrix = matrix(0, nrow=length(areanames), ncol=length(areanames))
	for (i in 1:length(areanames))
		{
		for (j in 1:length(areanames))
			{
			TF1 = events_df2$clado_dispersal_from == areanames[i]
			TF2 = events_df2$clado_dispersal_to == areanames[j]
			TF = (TF1 + TF2) == 2
			if (sum(TF) > 0)
				{
				counts_matrix[i,j] = sum(TF)
				} else {
				counts_matrix[i,j] = 0
				} # END if (sum(TF) > 0)
			} # END for (j in 1:length(areanames))
		} # END for (i in 1:length(areanames))

	counts_df = adf2(counts_matrix)
	names(counts_df) = actual_names
	row.names(counts_df) = actual_names
	
	return(counts_df)
	} # END count_clado_dispersal_events <- function(clado_events_table, areanames, actual_names=areanames)








# No longer assumes equal weights of ancestral areas, when ancestor
# occupies more than one area; assumes you have used dmat and 
# simulate_source_areas_ana_clado to add source area.
count_ana_dispersal_events <- function(ana_events_table, areanames, actual_names=areanames)
	{
	defaults='
	areanames = c("A", "B", "C", "D", "E", "F", "G")
	actual_names = c("SAm", "CAm", "Car", "NAm", "AF", "EU", "OZ")
	ana_events_table = DECJ_ana_events_tables[[1]]
	'

	if ("ana_dispersal_from" %in% names(ana_events_table) == FALSE)
		{
		txt = paste0("STOP ERROR in count_ana_dispersal_events(): input 'ana_events_table' must have a column 'ana_dispersal_from'. To add this column, run simulate_source_area_ana() or simulate_source_areas_ana_clado().")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		} # END error check
	
	TF = nchar(areanames) > 1
	if (any(TF))
		{
		txt = paste0("STOP ERROR in count_ana_dispersal_events(): all 'areanames' must be only 1 character long for the counting of events in stochastic mapping to work. You have violations for these 'areanames': ", paste(areanames[TF], sep=", "), ".")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		} # END if (any(TF))

	# Built up the transition matrix
	events_df2 = NULL
	extirp_df2 = NULL
	for (i in 1:length(areanames))
		{
		areaname = areanames[i]
		
		# Count the dispersal events
		TF = ana_events_table$dispersal_to == areaname
		TF[is.na(TF)] = FALSE
		tmp_table = ana_events_table[TF, ]
	
		# Get the ancestor range
		#nchars = nchar(tmp_table$current_rangetxt)
		time_bp = tmp_table$abs_event_time
		event_type = tmp_table$event_type
		ancrange = tmp_table$current_rangetxt
		decrange = tmp_table$new_rangetxt
		ana_dispersal_from = tmp_table$ana_dispersal_from
		ana_dispersal_to = tmp_table$dispersal_to
		events_df = adf2(cbind(time_bp, event_type, ancrange, decrange, ana_dispersal_from, ana_dispersal_to))
		events_df
		names(events_df) = c("time_bp", "event_type", "ancrange", "decrange", "ana_dispersal_from", "ana_dispersal_to")
		
		events_df2 = rbind(events_df2, events_df)
		
		# Count the extirpation events also!
		TF = ana_events_table$extirpation_from == areaname
		TF[is.na(TF)] = FALSE
		tmp_table = ana_events_table[TF, ]

		time_bp = tmp_table$abs_event_time
		event_type = tmp_table$event_type
		ancrange = tmp_table$current_rangetxt
		decrange = tmp_table$new_rangetxt
		extirpation_from = tmp_table$extirpation_from
		extirp_df = adf2(cbind(time_bp, event_type, ancrange, decrange, extirpation_from))
		names(extirp_df) = c("time_bp", "event_type", "ancrange", "decrange", "extirpation_from")
		extirp_df2 = rbind(extirp_df2, extirp_df)
		} # for (i in 1:areanames)
	
	if ( (is.null(events_df2) == TRUE) && (is.null(extirp_df2) == TRUE) )
		{
		ana_disp_counts = NULL
		ana_disp_counts$e_counts_df = NULL
		ana_disp_counts$a_counts_df = NULL
		ana_disp_counts$d_counts_df = NULL
		ana_disp_counts$counts_df = NULL
		return(ana_disp_counts)
		} # END if ( (is.null(events_df2) == TRUE) && (is.null(extirp_df2) == TRUE) )
	
	names(events_df2) = c("time_bp", "event_type", "ancrange", "decrange", "ana_dispersal_from", "ana_dispersal_to")
	events_df2
	names(extirp_df2) = c("time_bp", "event_type", "ancrange", "decrange", "extirpation_from")
	extirp_df2

	a_counts_matrix = matrix(0, nrow=length(areanames), ncol=length(areanames))
	d_counts_matrix = matrix(0, nrow=length(areanames), ncol=length(areanames))
	counts_matrix = matrix(0, nrow=length(areanames), ncol=length(areanames))
	e_counts_list = matrix(0, nrow=1, ncol=length(areanames))
	for (i in 1:length(areanames))
		{
		# Dispersal matrices
		for (j in 1:length(areanames))
			{
			TF1 = events_df2$ana_dispersal_from == areanames[i]
			TF2 = events_df2$ana_dispersal_to == areanames[j]
			TF = (TF1 + TF2) == 2
			if (sum(TF) > 0)
				{
				aTF = events_df2$event_type[TF] == "a"
				dTF = events_df2$event_type[TF] == "d"
				
				a_counts_matrix[i,j] = sum(aTF)
				d_counts_matrix[i,j] = sum(dTF)
				counts_matrix[i,j] = sum(TF)
				} else {
				a_counts_matrix[i,j] = 0
				d_counts_matrix[i,j] = 0
				counts_matrix[i,j] = 0
				} # END if (sum(TF) > 0)
			} # END for (j in 1:length(areanames))

		# Extinction/extirpation vector
		TF = events_df2$extirpation_from == areanames[i]
		if (sum(TF) > 0)
			{
			e_counts_list[,i] = sum(TF)
			} else {
			e_counts_list[,i] = 0
			}# END if (sum(TF) > 0)
		} # END for (i in 1:length(areanames))

	counts_df = adf2(counts_matrix)
	names(counts_df) = actual_names
	row.names(counts_df) = actual_names

	a_counts_df = adf2(a_counts_matrix)
	names(a_counts_df) = actual_names
	row.names(a_counts_df) = actual_names

	d_counts_df = adf2(d_counts_matrix)
	names(d_counts_df) = actual_names
	row.names(d_counts_df) = actual_names

	e_counts_df = adf2(e_counts_list)
	names(e_counts_list) = actual_names

	
	ana_disp_counts = NULL
	ana_disp_counts$e_counts_df = e_counts_df
	ana_disp_counts$a_counts_df = a_counts_df
	ana_disp_counts$d_counts_df = d_counts_df
	ana_disp_counts$counts_df = counts_df

	extract = '
	e_counts_df = ana_disp_counts$e_counts_df
	a_counts_df = ana_disp_counts$a_counts_df
	d_counts_df = ana_disp_counts$d_counts_df
	counts_df = ana_disp_counts$counts_df
	'
	
	return(ana_disp_counts)
	} # END count_ana_dispersal_events <- function(ana_events_table, areanames, actual_names=areanames)


# Convert e.g. AB->B,A to AB->A,B, since these are identical
uniquify_clado_events <- function(clado_events_table)
	{
	# Fix AB -> B,A

	# Eliminate the blank internal nodes
	TF = clado_events_table$clado_event_type != ""
	TF[is.na(clado_events_table$clado_event_type)] = FALSE
	clado_events_table = clado_events_table[TF,]
	clado_event_txt = clado_events_table$clado_event_txt
	dim(clado_events_table)

	# Lose something here?
	x = sapply(X=clado_event_txt, FUN=strsplit, split="->")
	# lx = unlist(lapply(X=x, FUN=length))
	# zero_TF = (lx == 0)
	# clado_events_table[zero_TF,]
	# Eliminate events with 0 hits first


	events_df = as.data.frame(matrix(unlist(x), ncol=2, byrow=TRUE), stringsAsFactors=FALSE)
	names(events_df) = c("ancrange", "decrange")
	dim(events_df)

	LRdf = as.data.frame(matrix(unlist(sapply(X=events_df$decrange, FUN=strsplit, split=",")), ncol=2, byrow=TRUE), stringsAsFactors=FALSE)
	names(LRdf) = c("L", "R")
	dim(LRdf)

	# Switch if left is bigger
	nchar_LRdf = as.data.frame(sapply(X=LRdf, FUN=nchar))
	names(nchar_LRdf) = c("L", "R")
	class(nchar_LRdf$L) = "numeric"
	class(nchar_LRdf$R) = "numeric"
	L_bigger_TF = nchar_LRdf$L > nchar_LRdf$R
	new_L = LRdf$R[L_bigger_TF]
	new_R = LRdf$L[L_bigger_TF]
	LRdf$L[L_bigger_TF] = new_L
	LRdf$R[L_bigger_TF] = new_R
	LRdf

	# In cases where left and right are equal, sort
	nchar_LRdf = as.data.frame(sapply(X=LRdf, FUN=nchar))
	names(nchar_LRdf) = c("L", "R")
	class(nchar_LRdf$L) = "numeric"
	class(nchar_LRdf$R) = "numeric"
	LR_equal_TF = nchar_LRdf$L == nchar_LRdf$R
	order_df = as.data.frame(matrix(apply(X=LRdf, MARGIN=1, FUN=order), ncol=2, byrow=TRUE))
	names(order_df) = c("L", "R")
	L_bigger_TF = order_df$L > order_df$R
	switch_TF = (L_bigger_TF + LR_equal_TF) == 2
	new_L = LRdf$R[switch_TF]
	new_R = LRdf$L[switch_TF]
	LRdf$L[switch_TF] = new_L
	LRdf$R[switch_TF] = new_R
	LRdf

	# Merge back
	LRtxt = apply(X=LRdf, MARGIN=1, paste, sep="", collapse=",")
	event_txt = apply(X=cbind(events_df$ancrange, LRtxt), MARGIN=1, paste, sep="", collapse="->")
	event_txt
	clado_events_table$clado_event_txt = event_txt
	
	return(clado_events_table)
	}





#######################################################
# Get huge events tables from saved stochastic maps
# Rdata files
#######################################################
get_huge_events_tables_from_BSM_Rdata <- function(BSM_tables_dir, BSM_fn_base, model_name="", suffix="")
	{
	defaults='
	# Working directory
	wd = "/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/_basic_example/"
	setwd(wd)

	# Biogeographical stochastic mappings (BSMs) saved
	BSM_tables_dir = "BSM_tables_M0_v1"
	BSM_fn_base = "BSM_M0_"
	model_name = "DEC"
	suffix = ".Rdata"
	' # END defaults

	# Necessary setting to avoid getting numbers etc. in the stochastic mapping output
	options(stringsAsFactors=FALSE)
	
	
	# Get the list of files
	searchstring = paste(BSM_fn_base, model_name, "_", sep="")#*.Rdata", sep="")
	tmpfns = list.files(path=BSM_tables_dir, pattern=searchstring)
	TF = grepl(pattern=suffix, x=tmpfns)
	BSM_fns = tmpfns[TF]
	BSM_fns
	num_BSM_fns = length(BSM_fns)
	
	# Huge tables, and indexes each BSM
	huge_cladogenetic_events_table = NULL
	huge_anagenetic_events_table = NULL
	ic = NULL
	ia = NULL


	# Get the output
	for (i in 1:num_BSM_fns)
		{
		BSM_fn = slashslash(paste(addslash(BSM_tables_dir), BSM_fns[i], sep=""))
	
		cat("\nExtracting file #", i, "/", num_BSM_fns, ": ", BSM_fn, "...", sep="")
	
		# Loads to BSM_results
		load(BSM_fn)
	
		# Is it stratified?
		if (class(BSM_results) == "data.frame")
			{
			# Non-stratified. Extract cladogenetic events table, 
			# calculate anagenetic events table
			BSM_strat_TF = FALSE
			stochastic_mapping_results = BSM_results
			master_table_cladogenetic_events = stochastic_mapping_results

			# Calculate anagenetic events table
			events_table = events_txt_list_into_events_table(events_txt_list=master_table_cladogenetic_events$anagenetic_events_txt_below_node, trtable=master_table_cladogenetic_events, recalc_abs_ages=TRUE)
			events_table
			table_w_anagenetic_events = events_table
		
			} else {
			# Stratified. Extract cladogenetic and anagenetic events tables
			BSM_strat_TF = TRUE
			stochastic_mapping_results = BSM_results
			master_table_cladogenetic_events = stochastic_mapping_results$master_table_cladogenetic_events
			table_w_anagenetic_events = stochastic_mapping_results$table_w_anagenetic_events
			} # END if (class(BSM_results) == "data.frame")
	
		# Add to huge tables
		huge_cladogenetic_events_table = rbind(huge_cladogenetic_events_table, master_table_cladogenetic_events)
		huge_anagenetic_events_table = rbind(huge_anagenetic_events_table, table_w_anagenetic_events)
		
		# Add to the ic, ia indexes
		if (is.null(master_table_cladogenetic_events) == FALSE)
			{
			tmp_ic = rep(i, nrow(master_table_cladogenetic_events))
			ic = c(ic, tmp_ic)
			} # END if (!is.null(master_table_cladogenetic_events))

		if (is.null(table_w_anagenetic_events) == FALSE)
			{
			tmp_ia = rep(i, nrow(table_w_anagenetic_events))
			ia = c(ia, tmp_ia)
			} # END if (!is.null(table_w_anagenetic_events))
		
		} # END
	
	# Add the indexes
	i = ic
	huge_cladogenetic_events_table = cbind(i, huge_cladogenetic_events_table)
	i = ia
	huge_anagenetic_events_table = cbind(i, huge_anagenetic_events_table)
	i = 0
	
	# Store as list of 2 tables
	huge_tables = NULL
	huge_tables$huge_cladogenetic_events_table = huge_cladogenetic_events_table
	huge_tables$huge_anagenetic_events_table = huge_anagenetic_events_table
	
	# To extract
	#huge_cladogenetic_events_table = huge_tables$huge_cladogenetic_events_table
	#huge_anagenetic_events_table = huge_tables$huge_anagenetic_events_table
	
	return(huge_tables)
	} # END get_huge_events_tables_from_BSM_Rdata






count_events_huge_tables <- function(huge_tables, timeperiod=NULL, area_abbr_to=NULL, area_abbr_from=NULL, BSM_i=NULL)
	{
	# Extract
	huge_cladogenetic_events_table = huge_tables$huge_cladogenetic_events_table
	huge_anagenetic_events_table = huge_tables$huge_anagenetic_events_table
	
	# If desired, count only within a timeperiod / stratum / time bin
	if (is.null(timeperiod) == FALSE)
		{
		# Error check
		if (length(timeperiod) != 2)
			{
			errortxt = paste("\n\nERROR in count_events_huge_tables(): timeperiod must either be NULL or be 2 numbers (min and max of time bin, in time before present).\n\n", sep="")
			cat(errortxt)
			stop(errortxt)
			}
		
		minT = min(timeperiod)
		maxT = max(timeperiod)
		
		# Subset cladogenetic events
		times_GT_min_TF = huge_cladogenetic_events_table$time_bp >= minT
		times_LT_max_TF = huge_cladogenetic_events_table$time_bp < maxT
		TF = (times_GT_min_TF + times_LT_max_TF) == 2
		huge_cladogenetic_events_table = huge_cladogenetic_events_table[TF,]
		
		# Subset anagenetic events
		times_GT_min_TF = huge_anagenetic_events_table$abs_event_time >= minT
		times_LT_max_TF = huge_anagenetic_events_table$abs_event_time < maxT
		TF = (times_GT_min_TF + times_LT_max_TF) == 2
		huge_anagenetic_events_table = huge_anagenetic_events_table[TF,]
		} # END if (is.null(timeperiod) == FALSE)


	# If desired, count events only dispersing to the area in area_abbr_to
	# (These are just d, j events)
	if (is.null(area_abbr_to) == FALSE)
		{
		# Subset cladogenetic events
		TF = huge_cladogenetic_events_table$clado_dispersal_to %in% area_abbr_to
		huge_cladogenetic_events_table = huge_cladogenetic_events_table[TF,]

		# Subset anagenetic events
		TF = huge_anagenetic_events_table$dispersal_to %in% area_abbr_to
		huge_anagenetic_events_table = huge_anagenetic_events_table[TF,]
		
		} # END if (is.null(area_abbr_to) == FALSE)

	# If desired, count events only dispersing to the area in area_abbr_from
	# (These are just d, j events)
	if (is.null(area_abbr_from) == FALSE)
		{
		# Subset cladogenetic events by the ancestor state before cladogenesis
		# (first, get clado_event_from events)
		tmpsplit1 <- function(x)
			{
			clado_event_from = strsplit(x=x, split="->")[[1]][1]
			return(clado_event_from)
			}
		
		clado_event_from = sapply(X=huge_cladogenetic_events_table$clado_event_txt, FUN=tmpsplit1) 
		names(clado_event_from)=NULL
		clado_event_from[is.na(clado_event_from)] = ""
		clado_event_from
		
		TF = clado_event_from %in% area_abbr_from
		huge_cladogenetic_events_table = huge_cladogenetic_events_table[TF,]

		# Subset anagenetic events
		ana_event_from = sapply(X=huge_anagenetic_events_table$event_txt, FUN=tmpsplit1) 
		names(ana_event_from)=NULL
		ana_event_from[is.na(ana_event_from)] = ""
		ana_event_from
		
		TF = ana_event_from %in% area_abbr_from
		huge_anagenetic_events_table = huge_anagenetic_events_table[TF,]
		
		} # END if (is.null(area_abbr_from) == FALSE)

	# If desired, count events only for a certain BSM realization index #i
	if (is.null(BSM_i) == FALSE)
		{
		# Subset cladogenetic events
		TF = huge_cladogenetic_events_table$i %in% BSM_i
		huge_cladogenetic_events_table = huge_cladogenetic_events_table[TF,]

		# Subset anagenetic events
		TF = huge_anagenetic_events_table$i %in% BSM_i
		huge_anagenetic_events_table = huge_anagenetic_events_table[TF,]
		}
	
	
	
	# Cladogenetic events
	# y events
	TF = grepl(pattern="\\(y\\)", x=huge_cladogenetic_events_table$clado_event_type)
	ycount = sum(TF, na.rm=TRUE)

	# s events
	TF = grepl(pattern="\\(s\\)", x=huge_cladogenetic_events_table$clado_event_type)
	scount = sum(TF, na.rm=TRUE)

	# v events
	TF = grepl(pattern="\\(v\\)", x=huge_cladogenetic_events_table$clado_event_type)
	vcount = sum(TF, na.rm=TRUE)

	# j events
	TF = grepl(pattern="\\(j\\)", x=huge_cladogenetic_events_table$clado_event_type)
	jcount = sum(TF, na.rm=TRUE)
	
	# Anagenetic events
	# d events
	TF = huge_anagenetic_events_table$event_type == "d"
	dcount = sum(TF, na.rm=TRUE)

	# e events
	TF = huge_anagenetic_events_table$event_type == "e"
	ecount = sum(TF, na.rm=TRUE)

	# a events
	TF = huge_anagenetic_events_table$event_type == "a"
	acount = sum(TF, na.rm=TRUE)
	
	events_count = matrix(data=c(dcount, ecount, acount, ycount, scount, vcount, jcount), nrow=1)
	events_count = as.data.frame(events_count)
	names(events_count) = c("d", "e", "a", "y", "s", "v", "j")
	events_count
		
	return(events_count)
	} # END count_events_huge_tables <- function(huge_tables, timeperiod=NULL, area_abbr_to=NULL, area_abbr_from=NULL, BSM_i=NULL)



# Copy of count_ana_clado_events (for back-compatibility; it DOES do jumps and anagenesis)
count_clado_events_nonjump <- function(clado_events_tables, ana_events_tables, areanames, actual_names)
	{
	count_ana_clado_events(clado_events_tables, ana_events_tables, areanames, actual_names)
	} # END count_clado_events_nonjump <- function(clado_events_tables, ana_events_tables, areanames, actual_names)






#######################################################
# Count anagenetic and cladogenetic events across
# a collection of BSMs
#######################################################

count_ana_clado_events <- function(clado_events_tables, ana_events_tables, areanames, actual_names)
	{
	defaults='
	areanames = c("A", "B", "C", "D", "E", "F", "G")
	actual_names = c("SAm", "CAm", "Car", "NAm", "AF", "EU", "OZ")
	clado_events_table = clado_events_tables[[1]]
	'
	
	TF = nchar(areanames) > 1
	if (any(TF))
		{
		txt = paste0("STOP ERROR: all 'areanames' must be only 1 character long for the counting of events in stochastic mapping to work. You have violations for these 'areanames': ", paste(areanames[TF], sep=", "), ".")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		} # END if (any(TF))


	# Error check
	if (length(areanames) != length(actual_names))
		{
		txt = paste0("STOP ERROR in count_ana_clado_events(): areanames and actual_names must have the same length. Instead, length(areanames)=", length(areanames), ", length(actual_names)=", length(actual_names), ".")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		} # END if (length(areanames) != length(actual_names))
	
	# Initialize
	dimdata = c(length(areanames), length(areanames))
	dims = c(dimdata, length(clado_events_tables))
	founder_counts_cube = array(NA, dim=dims)
	founder_totals_list = rep(NA, length(clado_events_tables))

	dimdata = c(length(areanames), length(areanames))
	dims = c(dimdata, length(ana_events_tables))
	ana_counts_cube = array(NA, dim=dims)
	ana_totals_list = rep(NA, length(ana_events_tables))
	a_counts_cube = array(NA, dim=dims)
	a_totals_list = rep(NA, length(ana_events_tables))
	d_counts_cube = array(NA, dim=dims)
	d_totals_list = rep(NA, length(ana_events_tables))

	e_counts_rectangle = matrix(NA, nrow=length(ana_events_tables), ncol=length(areanames))
	e_totals_list = rep(NA, length(ana_events_tables))

	
	vicariance_clado_events_table_bigTable = NULL
	subsetSymp_clado_events_table_bigTable = NULL
	sympatry_clado_events_table_bigTable = NULL

	vicariance_totals_list = rep(NA, length(clado_events_tables))
	subsetSymp_totals_list = rep(NA, length(clado_events_tables))
	sympatry_totals_list = rep(NA, length(clado_events_tables))


	# Make a list of each non-jump cladogenesis event
	for (i in 1:length(clado_events_tables))
		{
		# Cladogenesis events for this i
		clado_events_table = clado_events_tables[[i]]

		# Convert e.g. AB->B,A to AB->A,B, since these are identical
		clado_events_table = uniquify_clado_events(clado_events_table)


		# Count vicariance events
		vicariance_clado_events_table = clado_events_table[clado_events_table$clado_event_type == "vicariance (v)",]
		vicariance_clado_events_table
		vic_table = table(vicariance_clado_events_table$clado_event_txt)
		names_vic_table = names(vic_table)
		vic_table
	
		if (length(vic_table) > 0)
			{
			vic_table2 = cbind(i, names(vic_table), vic_table)
			} else {
			vic_table2 = NULL
			} # END if (length(vic_table) > 0)
	
		vicariance_clado_events_table_bigTable = rbind(vicariance_clado_events_table_bigTable, vic_table2)
	
		# Count subset sympatry events
		subsetSymp_clado_events_table = clado_events_table[clado_events_table$clado_event_type == "subset (s)",]
		subsetSymp_clado_events_table
		sub_table = table(subsetSymp_clado_events_table$clado_event_txt)
		names_sub_table = names(sub_table)
		sub_table

		if (length(sub_table) > 0)
			{
			sub_table2 = cbind(i, names(sub_table), sub_table)
			} else {
			sub_table2 = NULL
			} # END if (length(sub_table) > 0)

		subsetSymp_clado_events_table_bigTable = rbind(subsetSymp_clado_events_table_bigTable, sub_table2)

		# Count sympatry events
		sympatry_clado_events_table = clado_events_table[clado_events_table$clado_event_type == "sympatry (y)",]
		sympatry_clado_events_table
		sym_table = table(sympatry_clado_events_table$clado_event_txt)
		names_sym_table = names(sym_table)
		sym_table
	
		if (length(sym_table) > 0)
			{
			sym_table2 = cbind(i, names(sym_table), sym_table)
			} else {
			sym_table2 = NULL
			} # END if (length(sym_table2) > 0)
	
		sympatry_clado_events_table_bigTable = rbind(sympatry_clado_events_table_bigTable, sym_table2)

		} # END for (i in 1:length(clado_events_tables))

	
	# Defaults
	unique_vic_events = NULL
	unique_sub_events = NULL
	unique_sym_events = NULL

	# Populate if non-null
	if (is.null(vicariance_clado_events_table_bigTable) == FALSE)
		{
		unique_vic_events = sort(unique(vicariance_clado_events_table_bigTable[,2]))
		unique_vic_events
		} # END if (is.null(vicariance_clado_events_table_bigTable) == FALSE)
	if (is.null(subsetSymp_clado_events_table_bigTable) == FALSE)
		{
		unique_sub_events = sort(unique(subsetSymp_clado_events_table_bigTable[,2]))
		unique_sub_events
		} # END if (is.null(vicariance_clado_events_table_bigTable) == FALSE)
	if (is.null(sympatry_clado_events_table_bigTable) == FALSE)
		{
		unique_sym_events = sort(unique(sympatry_clado_events_table_bigTable[,2]))
		unique_sym_events
		} # END if (is.null(vicariance_clado_events_table_bigTable) == FALSE)

	unique_vic_counts = matrix(NA, nrow=length(clado_events_tables), ncol=length(unique_vic_events))
	unique_sub_counts = matrix(NA, nrow=length(clado_events_tables), ncol=length(unique_sub_events))
	unique_sym_counts = matrix(NA, nrow=length(clado_events_tables), ncol=length(unique_sym_events))

	cat("\ncount_ana_clado_events() is counting events over Biogeographical Stochastic Map (BSM) #:\n")
	for (i in 1:length(clado_events_tables))
	#for (i in 1:26)
		{
		cat(i, " ", sep="")
		
		# Extract the cladogenetic and anagenetic tables for this BSM
		clado_events_table = clado_events_tables[[i]]
		ana_events_table = ana_events_tables[[i]]
	
		# Count founder-events
		founder_counts_df = count_clado_dispersal_events(clado_events_table, areanames=areanames, actual_names=actual_names)
		if (is.null(founder_counts_df) == TRUE)
			{
			founder_counts_df = matrix(data=0, nrow=length(areanames), ncol=length(areanames))
			} # END if (is.null(founder_counts_df) == TRUE)

		# Fill in the jump-dispersal counts
		founder_counts_cube[,,i] = as.matrix(founder_counts_df)
		founder_totals_list[i] = sum(founder_counts_df)
		
		# Count range-expansion events
		if ( (length(ana_events_table) > 1) || (is.na(ana_events_table)==FALSE) )
			{
			# Parse out the different sorts of anagenetic dispersal
			#ana_counts_df = 
			ana_disp_counts = count_ana_dispersal_events(ana_events_table, areanames=areanames, actual_names=actual_names)
			
			e_counts_df = ana_disp_counts$e_counts_df
			a_counts_df = ana_disp_counts$a_counts_df
			d_counts_df = ana_disp_counts$d_counts_df
			ana_counts_df = ana_disp_counts$counts_df
			
			#round(ana_counts_df, 1)
			} # END if ( (length(ana_events_table) > 1) || (is.na(ana_events_table)==FALSE) )
	
		if ( (length(ana_events_table) > 1) || (is.na(ana_events_table)==FALSE) )
			{
			if (!is.null(ana_counts_df))
				{
				ana_counts_cube[,,i] = as.matrix(ana_counts_df)
				ana_totals_list[i] = sum(ana_counts_df)
				} else {
				ana_counts_cube[,,i] = 0
				ana_totals_list[i] = 0
				} # END if (!is.null(ana_counts_df))

			if (!is.null(a_counts_df))
				{
				a_counts_cube[,,i] = as.matrix(a_counts_df)
				a_totals_list[i] = sum(a_counts_df)
				} else {
				a_counts_cube[,,i] = 0
				a_totals_list[i] = 0
				} # END if (!is.null(a_counts_df))

			if (!is.null(d_counts_df))
				{
				d_counts_cube[,,i] = as.matrix(d_counts_df)
				d_totals_list[i] = sum(d_counts_df)
				} else {
				d_counts_cube[,,i] = 0
				d_totals_list[i] = 0
				} # END if (!is.null(d_counts_df))

			if (!is.null(e_counts_df))
				{
				e_counts_rectangle[i,] = as.matrix(e_counts_df)
				e_totals_list[i] = sum(e_counts_df)
				} else {
				e_counts_rectangle[i,] = 0
				e_totals_list[i] = 0
				} # END if (!is.null(e_counts_df))

			} # END if ( (length(ana_events_table) > 1) || (is.na(ana_events_table)==FALSE) )
	
		################################################################
		# Convert e.g. AB->B,A to AB->A,B, since these are identical
		################################################################
		clado_events_table = uniquify_clado_events(clado_events_table)

		# Count vicariance events
		vicariance_clado_events_table = clado_events_table[clado_events_table$clado_event_type == "vicariance (v)",]
		vicariance_clado_events_table
		if (nrow(vicariance_clado_events_table) > 0)
			{
			vic_table = table(vicariance_clado_events_table$clado_event_txt)
			names_vic_table = names(vic_table)
			#vic_table

			#ordernums = match(x=names_vic_table, table=unique_vic_events)
			#ordernums = ordernums[is.na(ordernums) == FALSE]
			event_in_uniqs_TF = unique_vic_events %in% names_vic_table
			unique_vic_counts[i, event_in_uniqs_TF] = vic_table
			vicariance_totals_list[i] = sum(unique_vic_counts[i,event_in_uniqs_TF])
			} else {
			vic_table = NA
			vicariance_totals_list[i] = 0
			}
		
		# Count subset sympatry events
		subsetSymp_clado_events_table = clado_events_table[clado_events_table$clado_event_type == "subset (s)",]
		subsetSymp_clado_events_table
		if (nrow(subsetSymp_clado_events_table) > 0)
			{
			sub_table = table(subsetSymp_clado_events_table$clado_event_txt)
			names_sub_table = names(sub_table)
			#ordernums = match(x=names_sub_table, table=unique_sub_events)
			#ordernums = ordernums[is.na(ordernums) == FALSE]
			event_in_uniqs_TF = unique_sub_events %in% names_sub_table
			unique_sub_counts[i, event_in_uniqs_TF] = sub_table
			subsetSymp_totals_list[i] = sum(unique_sub_counts[i,event_in_uniqs_TF])
			} else {
			sub_table = NA
			subsetSymp_totals_list[i] = 0
			}
		sub_table
	
		# Count sympatry events
		sympatry_clado_events_table = clado_events_table[clado_events_table$clado_event_type == "sympatry (y)",]
		sympatry_clado_events_table
		if (nrow(sympatry_clado_events_table) > 0)
			{
			sym_table = table(sympatry_clado_events_table$clado_event_txt)
			names_sym_table = names(sym_table)
			#sym_table
	
			#ordernums = match(x=names_sym_table, table=unique_sym_events)
			#ordernums = ordernums[is.na(ordernums) == FALSE]
			event_in_uniqs_TF = unique_sym_events %in% names_sym_table
			unique_sym_counts[i, event_in_uniqs_TF] = sym_table
			sympatry_totals_list[i] = sum(unique_sym_counts[i,event_in_uniqs_TF])
			} else {
			sym_table = NA
			sympatry_totals_list[i] = 0
			}
		} # END for (i in 1:length(clado_events_tables))


	unique_vic_counts = adf2(unique_vic_counts)
	names(unique_vic_counts) = unique_vic_events
	if (ncol(unique_vic_counts) > 0)
		{
		unique_vic_counts[is.na(unique_vic_counts)] = 0
		} # END if (nrow(unique_vic_events) > 0)

	unique_sub_counts = adf2(unique_sub_counts)
	names(unique_sub_counts) = unique_sub_events
	

	if (ncol(unique_sub_counts) > 0)
		{
		unique_sub_counts[is.na(unique_sub_counts)] = 0
		} # END if (nrow(unique_sub_counts) > 0)

	unique_sym_counts = adf2(unique_sym_counts)
	names(unique_sym_counts) = unique_sym_events
	if (ncol(unique_sym_counts) > 0)
		{
		unique_sym_counts[is.na(unique_sym_counts)] = 0
		} # END if (nrow(unique_sym_counts) > 0)

	head(unique_vic_counts)
	head(unique_sub_counts)
	head(unique_sym_counts)

	round(apply(X=unique_vic_counts, MARGIN=2, FUN=mean), 1)
	round(apply(X=unique_vic_counts, MARGIN=2, FUN=sd), 1)

	round(apply(X=unique_sub_counts, MARGIN=2, FUN=mean), 1)
	round(apply(X=unique_sub_counts, MARGIN=2, FUN=sd), 1)

	round(apply(X=unique_sym_counts, MARGIN=2, FUN=mean), 1)
	round(apply(X=unique_sym_counts, MARGIN=2, FUN=sd), 1)

	# Save the cladogenesis events
	vic_counts_means = apply(X=unique_vic_counts, MARGIN=2, FUN=mean)
	vic_counts_sds = apply(X=unique_vic_counts, MARGIN=2, FUN=sd)
	vic_counts_sums = apply(X=unique_vic_counts, MARGIN=2, FUN=sum)

	sub_counts_means = apply(X=unique_sub_counts, MARGIN=2, FUN=mean)
	sub_counts_sds = apply(X=unique_sub_counts, MARGIN=2, FUN=sd)
	sub_counts_sums = apply(X=unique_sub_counts, MARGIN=2, FUN=sum)

	sym_counts_means = apply(X=unique_sym_counts, MARGIN=2, FUN=mean)
	sym_counts_sds = apply(X=unique_sym_counts, MARGIN=2, FUN=sd)
	sym_counts_sums = apply(X=unique_sym_counts, MARGIN=2, FUN=sum)
	
	
	sum(sym_counts_means)
	sum(sub_counts_means)
	sum(vic_counts_means)


	
	# Save results
	counts_list = NULL
	
	# Fix NAs
	founder_counts_cube[is.na(founder_counts_cube)] = 0
	ana_counts_cube[is.na(ana_counts_cube)] = 0
	a_counts_cube[is.na(a_counts_cube)] = 0
	d_counts_cube[is.na(d_counts_cube)] = 0
	e_counts_rectangle[is.na(e_counts_rectangle)] = 0
	unique_sub_counts[is.na(unique_sub_counts)] = 0
	unique_vic_counts[is.na(unique_vic_counts)] = 0
	unique_sym_counts[is.na(unique_sym_counts)] = 0
	
	founder_totals_list[is.na(founder_totals_list)] = 0
	ana_totals_list[is.na(ana_totals_list)] = 0
	a_totals_list[is.na(a_totals_list)] = 0
	d_totals_list[is.na(d_totals_list)] = 0
	e_totals_list[is.na(e_totals_list)] = 0
	subsetSymp_totals_list[is.na(subsetSymp_totals_list)] = 0
	vicariance_totals_list[is.na(vicariance_totals_list)] = 0
	sympatry_totals_list[is.na(sympatry_totals_list)] = 0

	# Make an ALL dispersal events cube
	all_dispersals_counts_cube = founder_counts_cube
	all_dispersals_totals_list = founder_totals_list
	
	# Make an ALL anagenetic dispersal events cube
	anagenetic_dispersals_counts_cube = d_counts_cube
	anagenetic_dispersals_totals_list = d_totals_list
	
	for (i in 1:length(all_dispersals_totals_list))
		{
		all_dispersals_counts_cube[,,i] = all_dispersals_counts_cube[,,i] + a_counts_cube[,,i] + d_counts_cube[,,i]
		all_dispersals_totals_list[i] = all_dispersals_totals_list[i] + a_totals_list[i] + d_totals_list[i]
		anagenetic_dispersals_counts_cube[,,i] = anagenetic_dispersals_counts_cube[,,i] + a_counts_cube[,,i]
		anagenetic_dispersals_totals_list[i] = anagenetic_dispersals_totals_list[i] + a_totals_list[i]
		}
	
	# Events by category per BSM
	counts_list$all_dispersals_counts_cube = all_dispersals_counts_cube
	counts_list$anagenetic_dispersals_counts_cube = anagenetic_dispersals_counts_cube
	counts_list$founder_counts_cube = founder_counts_cube
	counts_list$ana_counts_cube = ana_counts_cube
	counts_list$a_counts_cube = a_counts_cube
	counts_list$d_counts_cube = d_counts_cube
	counts_list$e_counts_rectangle = e_counts_rectangle
	counts_list$unique_sub_counts = unique_sub_counts
	counts_list$unique_vic_counts = unique_vic_counts
	counts_list$unique_sym_counts = unique_sym_counts
	
	# Totals by category, per BSM
	counts_list$all_dispersals_totals_list = all_dispersals_totals_list
	counts_list$anagenetic_dispersals_totals_list = anagenetic_dispersals_totals_list
	counts_list$founder_totals_list = founder_totals_list
	counts_list$ana_totals_list = ana_totals_list
	counts_list$a_totals_list = a_totals_list
	counts_list$d_totals_list = d_totals_list
	counts_list$e_totals_list = e_totals_list
	counts_list$subsetSymp_totals_list = subsetSymp_totals_list
	counts_list$vicariance_totals_list = vicariance_totals_list
	counts_list$sympatry_totals_list = sympatry_totals_list

	clado_totals_list = founder_totals_list + subsetSymp_totals_list + vicariance_totals_list + sympatry_totals_list
	
	all_totals_list = founder_totals_list + subsetSymp_totals_list + vicariance_totals_list + sympatry_totals_list + ana_totals_list
	counts_list$clado_totals_list = clado_totals_list
	counts_list$all_totals_list = all_totals_list
	clado_totals_list
	all_totals_list
	# Means, SDs, sums across all BSMs
	
	means = c(mean(founder_totals_list), mean(a_totals_list), mean(d_totals_list), mean(e_totals_list), mean(subsetSymp_totals_list), mean(vicariance_totals_list), mean(sympatry_totals_list), mean(all_dispersals_totals_list), mean(anagenetic_dispersals_totals_list), mean(ana_totals_list), mean(clado_totals_list), mean(all_totals_list))
	sds = c(sd(founder_totals_list), sd(a_totals_list), sd(d_totals_list), sd(e_totals_list), sd(subsetSymp_totals_list), sd(vicariance_totals_list), sd(sympatry_totals_list), sd(all_dispersals_totals_list), sd(anagenetic_dispersals_totals_list), sd(ana_totals_list), sd(clado_totals_list), sd(all_totals_list))
	sums = c(sum(founder_totals_list), sum(a_totals_list), sum(d_totals_list), sum(e_totals_list), sum(subsetSymp_totals_list), sum(vicariance_totals_list), sum(sympatry_totals_list), sum(all_dispersals_totals_list), sum(anagenetic_dispersals_totals_list), sum(ana_totals_list), sum(clado_totals_list), sum(all_totals_list))
	
	# Summarize directionality of dispersal
	all_dispersals_counts_fromto_means = adf2(apply(X=all_dispersals_counts_cube, MARGIN=c(1,2), FUN=mean))
	ana_dispersals_counts_fromto_means = adf2(apply(X=anagenetic_dispersals_counts_cube, MARGIN=c(1,2), FUN=mean))
	founder_counts_fromto_means = adf2(apply(X=founder_counts_cube, MARGIN=c(1,2), FUN=mean))
	ana_events_counts_fromto_means = adf2(apply(X=ana_counts_cube, MARGIN=c(1,2), FUN=mean))
	a_counts_fromto_means = adf2(apply(X=a_counts_cube, MARGIN=c(1,2), FUN=mean))
	d_counts_fromto_means = adf2(apply(X=d_counts_cube, MARGIN=c(1,2), FUN=mean))

	all_dispersals_counts_fromto_sds = adf2(apply(X=all_dispersals_counts_cube, MARGIN=c(1,2), FUN=sd))
	ana_dispersals_counts_fromto_sds = adf2(apply(X=anagenetic_dispersals_counts_cube, MARGIN=c(1,2), FUN=sd))
	founder_counts_fromto_sds = adf2(apply(X=founder_counts_cube, MARGIN=c(1,2), FUN=sd))
	ana_counts_fromto_sds = adf2(apply(X=ana_counts_cube, MARGIN=c(1,2), FUN=sd))
	a_counts_fromto_sds = adf2(apply(X=a_counts_cube, MARGIN=c(1,2), FUN=sd))
	d_counts_fromto_sds = adf2(apply(X=d_counts_cube, MARGIN=c(1,2), FUN=sd))
	
	names(all_dispersals_counts_fromto_means) = actual_names
	row.names(all_dispersals_counts_fromto_means) = actual_names
	names(ana_dispersals_counts_fromto_means) = actual_names
	row.names(ana_dispersals_counts_fromto_means) = actual_names
	names(founder_counts_fromto_means) = actual_names
	row.names(founder_counts_fromto_means) = actual_names
	names(a_counts_fromto_means) = actual_names
	row.names(a_counts_fromto_means) = actual_names
	names(d_counts_fromto_means) = actual_names
	row.names(d_counts_fromto_means) = actual_names

	names(all_dispersals_counts_fromto_sds) = actual_names
	row.names(all_dispersals_counts_fromto_sds) = actual_names
	names(ana_dispersals_counts_fromto_sds) = actual_names
	row.names(ana_dispersals_counts_fromto_sds) = actual_names
	names(founder_counts_fromto_sds) = actual_names
	row.names(founder_counts_fromto_sds) = actual_names
	names(a_counts_fromto_sds) = actual_names
	row.names(a_counts_fromto_sds) = actual_names
	names(d_counts_fromto_sds) = actual_names
	row.names(d_counts_fromto_sds) = actual_names
	
	counts_list$all_dispersals_counts_fromto_means = all_dispersals_counts_fromto_means
	counts_list$ana_dispersals_counts_fromto_means = ana_dispersals_counts_fromto_means
	counts_list$founder_counts_fromto_means = founder_counts_fromto_means
	counts_list$a_counts_fromto_means = a_counts_fromto_means
	counts_list$d_counts_fromto_means = d_counts_fromto_means

	counts_list$all_dispersals_counts_fromto_sds = all_dispersals_counts_fromto_sds
	counts_list$ana_dispersals_counts_fromto_sds = ana_dispersals_counts_fromto_sds
	counts_list$founder_counts_fromto_sds = founder_counts_fromto_sds
	counts_list$a_counts_fromto_sds = a_counts_fromto_sds
	counts_list$d_counts_fromto_sds = d_counts_fromto_sds
	
	cat("\n\nRange-switching dispersal (all observed 'a' dispersals):\n")
	cat("\nmeans:\n")
	print(conditional_format_table(counts_list$a_counts_fromto_means))
	cat("\nstandard deviations:\n")
	print(conditional_format_table(counts_list$a_counts_fromto_sds))

	cat("\n\nRange-expansion dispersal (all observed 'd' dispersals):\n")
	cat("\nmeans:\n")
	print(conditional_format_table(counts_list$d_counts_fromto_means))
	cat("\nstandard deviations:\n")
	print(conditional_format_table(counts_list$d_counts_fromto_sds))

	cat("\n\nAnagenetic dispersal (mean of all observed anagenetic 'a' or 'd' dispersals):\n")
	cat("\nmeans:\n")
	print(conditional_format_table(counts_list$ana_dispersals_counts_fromto_means))
	cat("\nstandard deviations:\n")
	print(conditional_format_table(counts_list$ana_dispersals_counts_fromto_sds))

	cat("\n\nCladogenetic dispersal (mean of all observed jump 'j' dispersals):\n")
	cat("\nmeans:\n")
	print(conditional_format_table(counts_list$founder_counts_fromto_means))
	cat("\nstandard deviations:\n")
	print(conditional_format_table(counts_list$founder_counts_fromto_sds))

	cat("\n\nALL dispersal (mean of all observed anagenetic 'a', 'd' dispersals, PLUS cladogenetic founder/jump dispersal):\n")
	cat("\nmeans:\n")
	print(conditional_format_table(counts_list$all_dispersals_counts_fromto_means))
	cat("\nstandard deviations:\n")
	print(conditional_format_table(counts_list$all_dispersals_counts_fromto_sds))

	
	
	summary_counts_BSMs = rbind(means, sds, sums)
	summary_counts_BSMs = adf(summary_counts_BSMs)
	names(summary_counts_BSMs) = c("founder", "a", "d", "e", "subset", "vicariance", "sympatry", "ALL_disp", "ana_disp", "all_ana", "all_clado", "total_events")
	row.names(summary_counts_BSMs) = c("means", "stdevs", "sums")
	
	counts_list$summary_counts_BSMs = summary_counts_BSMs
	
	cat("\n\nSummary of event counts from ", length(clado_events_tables), " BSMs:\n", sep="")
	print(conditional_format_table(summary_counts_BSMs))
	
	# Code to extract output tables
	extract='
	all_dispersals_counts_cube = counts_list$all_dispersals_counts_cube
	anagenetic_dispersals_counts_cube = counts_list$anagenetic_dispersals_counts_cube
	founder_counts_cube = counts_list$founder_counts_cube
	ana_counts_cube = counts_list$ana_counts_cube
	a_counts_cube = counts_list$a_counts_cube
	d_counts_cube = counts_list$d_counts_cube
	e_counts_rectangle = counts_list$e_counts_rectangle
	unique_sub_counts = counts_list$unique_sub_counts
	unique_vic_counts = counts_list$unique_vic_counts
	unique_sym_counts = counts_list$unique_sym_counts

	all_dispersals_totals_list = counts_list$all_dispersals_totals_list
	anagenetic_dispersals_totals_list = counts_list$anagenetic_dispersals_totals_list
	founder_totals_list = counts_list$founder_totals_list
	ana_totals_list = counts_list$ana_totals_list
	a_totals_list = counts_list$a_totals_list
	d_totals_list = counts_list$d_totals_list
	e_totals_list = counts_list$e_totals_list
	subsetSymp_totals_list = counts_list$subsetSymp_totals_list
	vicariance_totals_list = counts_list$vicariance_totals_list
	sympatry_totals_list = counts_list$sympatry_totals_list
	

	summary_counts_BSMs = counts_list$summary_counts_BSMs
	'
	
	return(counts_list)
	} # END count_ana_clado_events <- function(clado_events_tables, ana_events_tables, areanames, actual_names)


calc_BSM_mean_node_states <- function(clado_events_tables, tr, numstates)
	{
	numBSMs = length(clado_events_tables)
	numtips = length(tr$tip.label)
	numnodes = numtips + tr$Nnode
	dims = c(numnodes, numstates, numBSMs)
	state_counts_1based = array(data=0, dim=dims)

	strat_TF = "SUBnode.type" %in% names(clado_events_tables[[1]])
	for (i in 1:numBSMs)
		{
		clado_events_table = clado_events_tables[[i]]
	
		if (strat_TF == TRUE)
			{
			# Stratified
			# Eliminate internal tips
			TF1 = clado_events_table$node.type == "internal"
			TF2 = clado_events_table$SUBnode.type == "tip"
			internal_tips_TF = (TF1 + TF2) == 2
			clado_events_table = clado_events_table[internal_tips_TF==FALSE,]
		
			# Order by node
			clado_events_table = clado_events_table[order(clado_events_table$node), ]

			nodestates = clado_events_table$sampled_states_AT_nodes
			# Store the states
			state_counts_temp = matrix(0, nrow=numnodes, ncol=numstates)
			for (j in 1:length(clado_events_table$node))
				{
				current_nodenum = clado_events_table$node[j]
				state_counts_temp[current_nodenum, nodestates[j]] = 1
				}
			} else {
			# NON-stratified
			nodestates = clado_events_table$sampled_states_AT_nodes
			# Store the states
			state_counts_temp = matrix(0, nrow=numnodes, ncol=numstates)
			for (j in 1:numnodes)
				{
				state_counts_temp[j, nodestates[j]] = 1
				}
			} # END if (strat_TF == TRUE)
		
		state_counts_1based[,,i] = state_counts_temp
		} # END for (i in 1:numBSMs)
	
	meanBSMprobs = apply(X=state_counts_1based, MARGIN=c(1,2), FUN=mean)
	sumBSMcounts = apply(X=state_counts_1based, MARGIN=c(1,2), FUN=sum)
	
	BSMstates_summary = NULL
	BSMstates_summary$state_counts_1based = state_counts_1based
	BSMstates_summary$meanBSMprobs = meanBSMprobs
	BSMstates_summary$sumBSMcounts = sumBSMcounts

	head(meanBSMprobs)
	dim(meanBSMprobs)
	colSums(meanBSMprobs)
	sum(colSums(meanBSMprobs))
	rowSums(meanBSMprobs)
	
	extract='
	state_counts_1based = BSMstates_summary$state_counts_1based
	meanBSMprobs = BSMstates_summary$meanBSMprobs
	sumBSMcounts = BSMstates_summary$sumBSMcounts
	'
	return(BSMstates_summary)
	} # END calc_BSM_mean_node_states <- function(clado_events_tables, tr, numstates)













#######################################################
# Histograms of event counts
#######################################################
hist_event_counts <- function(counts_list, pdffn="hist_event_counts.pdf", col="grey70")
	{
	title_cex = 1
	
	counts_tmp = c(counts_list$founder_totals_list, counts_list$ana_totals_list, counts_list$subsetSymp_totals_list, counts_list$vicariance_totals_list, counts_list$sympatry_totals_list)
	xmax = max(counts_tmp)
	xlims = c(0, max(pretty(counts_tmp)))
	count_names = seq(xlims[1], xlims[2])
	numBSMs = length(counts_list$all_totals_list)
	
	
	

	#pdffn = "hist_event_counts.pdf"
	pdffn = pdffn
	pdf(file=pdffn, width=4, height=8)

	par(mfrow=c(5,1))
	# c(bottom, left, top, right) 
	par(mar=c(2,4,0.5,1))
	
	ylabel = paste0("Freq. in ", numBSMs, " BSMs")
	#par(oma=c(0,0,0,0))
	#hist(counts_list$ana_totals_list, xlab=NULL, ylab=NULL, xlim=xlims, main="")
	mp = barplot(table2(counts_list$ana_totals_list, xlims=xlims), xlab=NULL, ylab=NULL, main="", col=col)
	axis(side=1, at=mp, labels=count_names, tick=FALSE, line=-0.75, cex.axis=0.8)
	#title("Anagenetic dispersal events", cex.main=title_cex, font.main=1)
	mtext(text="# anagenetic dispersal", side=2, cex=0.7, line=2.5)

	mp = barplot(table2(counts_list$sympatry_totals_list, xlims=xlims), xlab=NULL, ylab=NULL, main="", col=col)
	axis(side=1, at=mp, labels=count_names, tick=FALSE, line=-0.75, cex.axis=0.8)
	#title("Narrow sympatry", cex.main=title_cex, font.main=1)
	mtext(text="# narrow sympatry", side=2, cex=0.7, line=2.5)

	mp = barplot(table2(counts_list$subsetSymp_totals_list, xlims=xlims), xlab=NULL, ylab=NULL, main="", col=col)
	axis(side=1, at=mp, labels=count_names, tick=FALSE, line=-0.75, cex.axis=0.8)
	#title("Subset sympatry", cex.main=title_cex, font.main=1)
	mtext(text="# subset sympatry", side=2, cex=0.7, line=2.5)

	mp = barplot(table2(counts_list$vicariance_totals_list, xlims=xlims), xlab=NULL, ylab=NULL, main="", col=col)
	axis(side=1, at=mp, labels=count_names, tick=FALSE, line=-0.75, cex.axis=0.8)
	#title("Vicariance", cex.main=title_cex, font.main=1)
	mtext(text="# vicariance", side=2, cex=0.7, line=2.5)

	par(mar=c(3.5,4,0,1))
	mp = barplot(table2(counts_list$founder_totals_list, xlims=xlims), xlab=NULL, ylab=NULL, main="", col=col)
	axis(side=1, at=mp, labels=count_names, tick=FALSE, line=-0.75, cex.axis=0.8)
	#title("Founder-events", cex.main=title_cex, font.main=1)
	mtext(text="# founder events", side=2, cex=0.7, line=2.5)
	mtext(text=paste0("Event counts in each of ", numBSMs, " BSMs"), side=1, cex=0.7, line=2)
	
	dev.off()
	cmdstr = paste0("open ", pdffn)
	system(cmdstr)
	} # END hist_event_counts <- function(counts_list, pdffn="hist_event_counts.pdf")

# Max a table, with 0s for the integers that are unobserved
table2 <- function(event_counts, xlims)
	{
	count_names = seq(xlims[1], xlims[2])
	counts_for_barplot = matrix(0, nrow=1, ncol=length(count_names))
	names(counts_for_barplot) = count_names
	counts_for_barplot
	
	tmptable = table(as.integer(event_counts))
	indices = match(x=names(tmptable), table=count_names)
	counts_for_barplot[indices] = tmptable
	
	return(counts_for_barplot)
	}






#######################################################
# Check that ML ancestral state/range probabilities and
# the mean of the BSMs approximately line up
#######################################################
check_ML_vs_BSM <- function(res, clado_events_tables, model_name, tr=NULL, plot_each_node=FALSE, linreg_plot=TRUE, MultinomialCI=TRUE)
	{
	if (is.null(tr))
		{
		#tr = read.tree(res$inputs$trfn)
		tr = check_trfn(trfn=res$inputs$trfn)
		} # END if (is.null(tr))

	pdffn = paste0(model_name, "_ML_vs_BSM.pdf")
	pdf(file=pdffn, height=6, width=6)

	x_all = NULL
	y_all = NULL
	numtips = length(tr$tip.label)
	intnodenums = (numtips+1):(numtips+tr$Nnode)

	MLstateprobs = res$ML_marginal_prob_each_state_at_branch_top_AT_node
	numstates = ncol(MLstateprobs)
	BSMstates_summary = calc_BSM_mean_node_states(clado_events_tables, tr, numstates)
	state_counts_1based = BSMstates_summary$state_counts_1based
	meanBSMprobs = BSMstates_summary$meanBSMprobs
	sumBSMcounts = BSMstates_summary$sumBSMcounts

	node_history_samples_allNodes = NULL
	cat("\nCalculating BSM means for node #:", sep="")
	
	for (nodenum in intnodenums)
		{
		cat(nodenum, " ", sep="")
		
		node_history_samples = NULL
		for (i in 1:length(clado_events_tables))
			{
			clado_events_table = clado_events_tables[[i]]
			if (stratified == TRUE)
				{
				TF1 = clado_events_table$SUBnode.type != "tip"
				TF2 = clado_events_table$node.type != "tip"
				TF = (TF1+TF2)==2
				} else {
				TF = clado_events_table$node.type != "tip"
				}
			clado_events_table = clado_events_table[TF,]
			TF = clado_events_table$node == nodenum
			node_history_sample = clado_events_table[TF,]
			node_history_samples = rbind(node_history_samples, node_history_sample)
			} # END for (i in 1:length(clado_events_tables))
		head(node_history_samples)
		class(node_history_samples)
		dim(node_history_samples)

		# Save
		node_history_samples_allNodes = rbind(node_history_samples_allNodes, node_history_samples)
	
		node_history_samples$clado_event_txt
		node_history_samples$sampled_states_AT_nodes


		cbind(node_history_samples$sampled_states_AT_nodes, node_history_samples$sampled_states_AT_brbots)
		table(node_history_samples$sampled_states_AT_nodes)
		table(node_history_samples$sampled_states_AT_brbots)

		x=round(res$ML_marginal_prob_each_state_at_branch_top_AT_node[nodenum,],3)

		num_root_state = table(node_history_samples$sampled_states_AT_nodes)
		fract_root_state = num_root_state / sum(num_root_state)
		obs_BSM_probs_at_root = rep(0, length(x))
		obs_BSM_probs_at_root[as.numeric(names(fract_root_state))] = fract_root_state
	
		x=round(res$ML_marginal_prob_each_state_at_branch_top_AT_node[nodenum,],3)
		y=obs_BSM_probs_at_root
		x_all = c(x_all, x)
		y_all = c(y_all, y)
	
		if (plot_each_node == TRUE)
			{
			plot(x,y, xlim=c(0,1), ylim=c(0,1), xlab="ML marginal probs", ylab="BSM probs")
			segments(0,0,1,1)
			title(nodenum)
			} # END if (plot_each_node == TRUE)
		} # END for (nodenum in 20:37)
	cat("...done.\n\n")

	# Plot allnodes
	if (linreg_plot == FALSE)
		{
		plot(x_all, y_all, xlim=c(0,1), ylim=c(0,1), xlab="ML marginal probs", ylab="BSM probs")
		segments(0,0,1,1)
		title("all nodes")
		} # END if (linreg_plot == FALSE)

	if (linreg_plot == TRUE)
		{
		linear_regression_plot(x=x_all, y=y_all, tmppch=1, xlab="State probabilities under ML model", ylab="State probabilities as mean of BSMs", xlim=c(0,1), ylim=c(0,1))
		# Multiple R-squared:  0.9655,	Adjusted R-squared:  0.9655 
		title(paste0(model_name, ":\nML state probs vs. mean of BSMs"))
		} # END if (linreg_plot == TRUE)

	if (MultinomialCI == TRUE)
		{
		# Try doing confidence intervals
		require(MultinomialCI)
		sumBSMcounts2 = sumBSMcounts[-(1:numtips),]

		# Confident intervals for a particular node
		# CI95 = multinomialCI(x=sumBSMcounts2[1,], alpha=0.05, verbose=TRUE)
		# CI95

		# Add multinomial CI95s to plot
		cat("\nAdding CI95 segments for: ", sep="")
		for (i in 1:tr$Nnode)
			{
			cat(i+numtips, " ", sep="")
			tmpcounts = sumBSMcounts2[i,]
			CI95 = multinomialCI(x=tmpcounts, alpha=0.05, verbose=FALSE)
			#xvals = tmpcounts / sum(tmpcounts)
			xvals = MLstateprobs[-(1:numtips),][i,]
			segments(x0=xvals, x1=xvals, y0=CI95[,1], y1=CI95[,2])
			} # for (1 in 1:numnodes)
		} # END if (MultinomialCI == TRUE)
	cat("...done.\n\n")

	# Close PDF and open for viewer
	dev.off()
	cmdstr = paste("open ", pdffn)
	system(cmdstr)
	} # END check_ML_vs_BSM <- function(res, clado_events_tables, tr=NULL)


minmax_pretty <- function(x)
	{
	minmax = c( min(pretty(x)), max(pretty(x)) )
	return(minmax)
	}


#######################################################
# Linear regression plot, with stats
#######################################################
# If slope1=TRUE, subtract 1:1 slope from the line and test for differences

linear_regression_plot <- function(x, y, xlabel="x", ylabel="y", tmppch=".", pointscol="black", tmplinecol="black", tmplty=1, tmplwd=1, plottext=TRUE, legend_title="", textcol="black", legend_x="topleft", legend_y=NULL, xlim=minmax_pretty(x), ylim=minmax_pretty(y), increment_fraction=NULL, legend_cex=1, axis_cex=1, slope1=FALSE, intercept_in_legend=TRUE, add_to_plot=FALSE, printall=TRUE, ...)
	{
	# Make a linear regression plot
	model1 = lm(y~x)
	slope = model1$coefficients[2]
	intercept = model1$coefficients[1]
	
	if (printall)
		{
		print("Summary of model 1 (standard regression, y predicted by x)")
		print(summary(model1))
		}
	



	# Set up Legend
	if (is.null(increment_fraction))
		{
		increment_fraction = 0.05
		}

	# Legend positions
	if (legend_x == "topleft")
		{
		legend_x = min(xlim)
		legend_y = 1 * max(ylim)
		}
	
	increment = increment_fraction * (max(ylim) - min(ylim))
	
	
	# Plot the full plot, or just add points
	if (add_to_plot == FALSE)
		{
		plot(x, y, xlab=xlabel, ylab=ylabel, pch=tmppch, xlim=xlim, ylim=ylim, col=pointscol, cex.axis=axis_cex, ...)
		} else {
		points(x, y, pch=tmppch, xlim=xlim, ylim=ylim, col=pointscol, ...)
		}
	tmpx1 = min(x, na.rm=TRUE)
	tmpx2 = max(x, na.rm=TRUE)
	tmpy1 = slope*tmpx1 + intercept
	tmpy2 = slope*tmpx2 + intercept
	
	segments(x0=tmpx1, y0=tmpy1, x1=tmpx2, y1=tmpy2, col=tmplinecol, lty=tmplty, lwd=tmplwd)
	
	
	if (plottext)
		{
		# Subtract 1:1 slope if desired
		# (for Patrick Shih cyanobacteria analysis)
		if (slope1 == TRUE)
			{
			# Subtract 1:1 from the y values
			ynew = y-x
			model2 = lm(ynew~x)
			pval = summary(model2)$coefficients[2,4]

			print("Summary of model 2 (standard regression, (y-x) predicted by x; this means we've substracted out the 1:1 expected relationship.)")
			print(summary(model2))
			} else {
			pval = summary(model1)$coefficients[2,4]
			} # END if (slope1 == TRUE)
		
		# Start ypos for legend text at 0
		ypos = 0
		# Plot a legend title if desired
		if (legend_title != "")
			{
			text(x=legend_x, y=legend_y - ypos*increment, pos=4, labels=legend_title, cex=legend_cex, col=textcol)
			} # END if (legend_title != "")
	
		
		if (intercept_in_legend==TRUE)
			{		
			R2 = summary(model1)$r.squared
			slope = summary(model1)$coefficients[2,1]
			slopeSE = summary(model1)$coefficients[2,2]
			slope95 = 1.96*slopeSE

			intercept = summary(model1)$coefficients[1,1]
			interceptSE = summary(model1)$coefficients[1,2]
			intercept95 = 1.96*interceptSE
		
			R2txt = bquote(italic(R)^2 == .(format(R2, digits=3)))
			text(x=legend_x, y=legend_y - (ypos=ypos+1)*increment, pos=4, labels=R2txt, cex=legend_cex, col=textcol)
			slopetxt = bquote(italic(m) == .(format(slope, digits=3)) %+-% .(format(slope95, digits=3)))
			text(x=legend_x, y=legend_y - (ypos=ypos+1)*increment, pos=4, labels=slopetxt, cex=legend_cex, col=textcol)
			#intercepttxt = paste("intercept=", format(intercept, digits=3), " +/- ", format(intercept95, digits=3), sep="")
			intercepttxt = bquote(italic(b) == .(format(intercept, digits=3)) %+-% .(format(intercept95, digits=3)))
			text(x=legend_x, y=legend_y - (ypos=ypos+1)*increment, pos=4, labels=intercepttxt, cex=legend_cex, col=textcol)
		
			if (slope1 == TRUE)
				{
				#pvaltxt = paste("p = ", format(pval, digits=3), " (null: slope is 1:1)", sep="")
				pvaltxt = bquote(italic(p) == .(format(pval, digits=3))~" (null: slope is 1:1)")
				text(x=legend_x, y=legend_y - (ypos=ypos+1)*increment, pos=4, labels=pvaltxt, cex=legend_cex, col=textcol)
				} else {
				#pvaltxt = paste("p = ", format(pval, digits=3), sep="")
				pvaltxt = bquote(italic(p) == .(format(pval, digits=3)))
				text(x=legend_x, y=legend_y - (ypos=ypos+1)*increment, pos=4, labels=pvaltxt, cex=legend_cex, col=textcol)
				} # END if (slope1 == TRUE)

			} # END if (intercept_in_legend==TRUE)
		
# 		if (intercept_in_legend==TRUE)
# 			{
# 			txt_to_plot = paste(R2txt, slopetxt, intercepttxt, pvaltxt, sep="\n")
# 			} else {
# 			txt_to_plot = paste(R2txt, slopetxt, pvaltxt, sep="\n")
# 			}
		
		# http://stackoverflow.com/questions/3761410/how-can-i-plot-my-r-squared-value-on-my-scatterplot-using-r
		# bty suppresses box
		# print(legend_x)
		# print(legend_y)
		# Discounted
		#legend(x=legend_x, y=legend_y, bty="n", legend=txt_to_plot, cex=0.9)

		}
	
	if (printall == TRUE)	
		{
		print(model1)
		print(summary(model1))
		}
	
	
	
	return(model1)
	}




linear_regression_plot_OLD <- function(x, y, xlabel="x", ylabel="y", tmppch=".", printall=TRUE, tmplinecol="black", tmplty=1, tmplwd=1, plottext=TRUE, legend_x="topleft", legend_y=NULL, xlim=minmax_pretty(x), ylim=minmax_pretty(y), slope1=FALSE, ...)
	{
	# Make a linear regression plot
	model1 = lm(y~x)
	print("Summary of model 1 (standard regression, y predicted by x)")
	print(summary(model1))
	
	slope = model1$coefficients[2]
	intercept = model1$coefficients[1]
	
	plot(x, y, xlab=xlabel, ylab=ylabel, pch=tmppch, xlim=xlim, ylim=ylim, ...)
	tmpx1 = min(x, na.rm=TRUE)
	tmpx2 = max(x, na.rm=TRUE)
	tmpy1 = slope*tmpx1 + intercept
	tmpy2 = slope*tmpx2 + intercept
	
	segments(x0=tmpx1, y0=tmpy1, x1=tmpx2, y1=tmpy2, col=tmplinecol, lty=tmplty, lwd=tmplwd)
	
	
	if (plottext)
		{
		# Subtract 1:1 slope if desired
		# (for Patrick Shih cyanobacteria analysis)
		if (slope1 == TRUE)
			{
			# Subtract 1:1 from the y values
			ynew = y-x
			model2 = lm(ynew~x)
			pval = summary(model2)$coefficients[2,4]

			print("Summary of model 2 (standard regression, (y-x) predicted by x; this means we've substracted out the 1:1 expected relationship.)")
			print(summary(model2))
			} else {
			pval = summary(model1)$coefficients[2,4]
			} # END if (slope1 == TRUE)
		
		
		R2 = summary(model1)$r.squared
		slope = summary(model1)$coefficients[2,1]
		slopeSE = summary(model1)$coefficients[2,2]
		slope95 = 1.96*slopeSE
		
		R2txt = paste("R2 = ", format(R2, digits=3), sep="")
		slopetxt = paste("m=", format(slope, digits=3), " +/- ", format(slope95, digits=3), sep="")
		
		if (slope1 == TRUE)
			{
			pvaltxt = paste("p = ", format(pval, digits=3), " (null: slope is 1:1)", sep="")
			} else {
			pvaltxt = paste("p = ", format(pval, digits=3), sep="")
			} # END if (slope1 == TRUE)

		txt_to_plot = paste(R2txt, slopetxt, pvaltxt, sep="\n")
		
		# http://stackoverflow.com/questions/3761410/how-can-i-plot-my-r-squared-value-on-my-scatterplot-using-r
		# bty suppresses box
		# print(legend_x)
		# print(legend_y)
		legend(x=legend_x, y=legend_y, bty="n", legend=txt_to_plot, cex=0.9)

		}
	
	if (printall == TRUE)	
		{
		print(model1)
		print(summary(model1))
		}
	
	return(model1)
	} # END linear_regression_plot

