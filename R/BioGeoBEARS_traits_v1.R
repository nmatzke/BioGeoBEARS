
#######################################################
# Code to assist models where an evolving trait can
# influence dispersal ability
#######################################################

add_jts_to_BioGeoBEARS_run_object <- function(BioGeoBEARS_run_object, numtrait_states)
	{
	jts_txt_matrix = matrix(data="", nrow=numtrait_states, ncol=numtrait_states)
	for (jts_i in 1:numtrait_states)
		{
		for (jts_j in 1:numtrait_states)
			{
			newtxt = paste0("jt", jts_i, jts_j)
			jts_txt_matrix[jts_i,jts_j] = newtxt
			
			if (jts_i != jts_j)
				{
				tmprow = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d",]
				row.names(tmprow) = newtxt
				tmprow$type = "fixed"
				tmprow$init = 0.0
				tmprow$min = 0.0
				tmprow$max = 1.0
				tmprow$est = 0.0
				tmprow$desc = paste0("prop. j events where trait changes state ", jts_i, "->", jts_j)
				BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table = rbind(BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table, tmprow)
				}
			} # END for (jts_j in 1:numtrait_states)
		} # END for (jts_i in 1:numtrait_states)
	
	BioGeoBEARS_run_object$jts_txt_matrix = jts_txt_matrix
	
	return(BioGeoBEARS_run_object)
	}


# Add a trait (and the relevant parameters etc.

#' BioGeoBEARS_run_object A BioGeoBEARS_run_object.
#' traits_fn A traits filename. The traits should be in the same 
#' kind of format that is used for geography data. I.e., a 2-state 
#' trait will be formatted like a 2-area geography dataset would be.
add_trait_to_BioGeoBEARS_run_object <- function(BioGeoBEARS_run_object, trait_fn, block_allQs=TRUE, block_allQs_trait=TRUE)
	{
	defaults='
	# Now, load a trait dataset (trait "Fly"/"Non", for flightlessness)
	trait_fn = "flightlessness.txt"
	' # END defaults

	#######################################################
	# Load the trait
	#######################################################

	# Look at your geographic range data:
	trait = getranges_from_LagrangePHYLIP(lgdata_fn=trait_fn, block_allQs=block_allQs_trait)
	trait
	
	###################################
	# Error check on trait
	###################################
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=BioGeoBEARS_run_object$geogfn, block_allQs=block_allQs)
	species_from_geog = sort(row.names(tipranges@df))

	species_from_traits = sort(row.names(trait@df))
	matchTF1 = all(species_from_geog %in% species_from_traits)
	matchTF2 = all(species_from_traits %in% species_from_geog)
	
	if ((matchTF1 + matchTF2) != 2)
		{
		txt = "STOP ERROR in add_trait_to_BioGeoBEARS_run_object(). The species names in your traits file do not match the species names in your geography file. These must match exactly."
		
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		cat("Printing species from geography file:")
		cat("\n\n")
		species_from_geog
		cat("Printing species from traits file:")
		cat("\n\n")
		species_from_traits
		cat("\n\n")
		stop(txt)
		} # END if ((matchTF1 + matchTF2) != 2)
	###################################
	# END Error check on trait
	###################################
	
	
	# Add trait to BioGeoBEARS_run_object
	BioGeoBEARS_run_object$trait = trait
	

	#######################################################
	# Modify the BioGeoBEARS_model_object@params_table
	#######################################################
	ntrait_states = ncol(trait@df)
	ntrait_states

	#######################################################
	# Add Pmat for trait (transition rates between trait states)
	# (parameters t12, t21)
	#######################################################
	trait_Pmat_txt = matrix(data="0", nrow=ntrait_states, ncol=ntrait_states)
	param_row = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table[1,]
	BGB_trait_Pmat_params_table = NULL
	for (i in 1:ntrait_states)
		{
		for (j in 1:ntrait_states)
			{
			if (i==j)
				{
				trait_Pmat_txt[i,j] = "0"
				next()
				} else {
				diag_name = paste0("t", i, j)
				trait_Pmat_txt[i,j] = diag_name
				param_row$type = "fixed"
				param_row$init = 0.001
				param_row$est = 1
				param_row$min = 0.00000001
				param_row$max = 50
				param_row$desc = "trait transition rate"
				row.names(param_row) = diag_name
				} # END if (i==j)
			BGB_trait_Pmat_params_table = rbind(BGB_trait_Pmat_params_table, param_row)
			} # END for (j in 1:ntrait_states)
		} # END for (i in 1:ntrait_states)

	# Merge t12 and t21 into params table
	row_names_being_added = row.names(BGB_trait_Pmat_params_table)
	
	# Remove those rows, if already present
	current_row_names = row.names(BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table)
	keepTF = (current_row_names %in% row_names_being_added) == FALSE
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table[keepTF,]
	
	# Add the t12, t21 etc. rows
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table = rbind(BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table, BGB_trait_Pmat_params_table)

	# Add the Pmat for the trait to the BioGeoBEARS run object
	BioGeoBEARS_run_object$trait_Pmat_txt = trait_Pmat_txt


	#######################################################
	# Model for dispersal multiplier (m) when in different trait states
	# (parameter m1, m2)
	#######################################################
	mnames = paste("m", 1:ntrait_states, sep="")
	# Start with dummy rows of the right size
	BGB_trait_model_params_table = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table[1:ntrait_states,]
	rownames(BGB_trait_model_params_table) = mnames
	# Start all the multiplier parameters at 1
	BGB_trait_model_params_table$type = "fixed"
	BGB_trait_model_params_table$init = 1
	BGB_trait_model_params_table$est = 1

	for (i in 1:ntrait_states)
		{
		BGB_trait_model_params_table$desc[i] = paste0("trait-based dispersal rate multiplier when trait=", i-1)
		}

	# Merge m1 and m2 into params table
	row_names_being_added = row.names(BGB_trait_model_params_table)
	
	# Remove those rows, if already present
	current_row_names = row.names(BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table)
	keepTF = (current_row_names %in% row_names_being_added) == FALSE
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table[keepTF,]

  # Now add the row names
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table = rbind(BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table, BGB_trait_model_params_table)

	return(BioGeoBEARS_run_object)
	} # END add_trait_to_BioGeoBEARS_run_object <- function(BioGeoBEARS_run_object, trait_fn)



#######################################################
# Modify the base Qmat by the trait transition matrix
#######################################################
modify_Qmat_with_trait <- function(Qmat=NULL, BioGeoBEARS_run_object, numstates_geogtrait, areas_list, states_list, dispersal_multipliers_matrix, elist, force_sparse)
	{
	# Update BioGeoBEARS_model_object, just to be sure
	BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
	
	# Pmat for trait transition matrix
	trait_Pmat_txt = BioGeoBEARS_run_object$trait_Pmat_txt
	num_trait_states = ncol(trait_Pmat_txt)
	num_geog_states = numstates_geogtrait / num_trait_states
	
	# Trait multipliers
	trait_multiplier_rows_TF = grepl(x=BioGeoBEARS_model_object@params_table$desc, pattern="trait-based dispersal rate multiplier")
	trait_multiplier_rows_indices = (1:length(trait_multiplier_rows_TF))[trait_multiplier_rows_TF]
	ntrait_states = length(trait_multiplier_rows_indices)
	BGB_trait_model_params_table = BioGeoBEARS_model_object@params_table[trait_multiplier_rows_indices,]
	
	# Set the m values (dispersal multipliers, dependent on trait)
	m = BGB_trait_model_params_table$est
	#print("m:")
	#print(m)
	
	# Matrix for transitions in the trait
	trait_transition_rows_TF = grepl(x=BioGeoBEARS_model_object@params_table$desc, pattern="trait transition rate")
	trait_transition_rows_indices = (1:length(trait_transition_rows_TF))[trait_transition_rows_TF]
	BGB_trait_transition_params_table = BioGeoBEARS_model_object@params_table[trait_transition_rows_indices,]
	
	d = BioGeoBEARS_model_object@params_table["d","est"]
	a = BioGeoBEARS_model_object@params_table["a","est"]

				
	# Go through the 2x2 trait transition matrix, and
	# construct the big Rmat (which is then converted to Qmat)
	traitstate_i = 1
	traitstate_j = 1
	for (traitstate_i in 1:ntrait_states)
		{
		# For a particular trait state, apply the dispersal multipliers
		# trait_multiplier = m, basically
		# (for NJM: search on m = m, m=m, *m, *m[1])
		# This is equivalent to multiplying rates by *m[traitstate_i]
		trait_multiplier = BGB_trait_model_params_table$est[traitstate_i]
		#print("trait_multiplier:")
		#print(trait_multiplier)
	
		
		# Generate the Qmat, unless given!
		if (is.null(Qmat))
			{
			# trait_multiplier is in dmat and amat, so it's in the output matrices also
			dmat_times_d = trait_multiplier * dispersal_multipliers_matrix * matrix(d, nrow=length(areas_list), ncol=length(areas_list))
			amat = trait_multiplier * dispersal_multipliers_matrix * matrix(a, nrow=length(areas_list), ncol=length(areas_list))
			
			# By default, the COO version comes out BACKWARDS; this is accounted for elsewhere in regular 
			# sparse matrices, but here we have to reverse it!!
			Qmat_tmp = rcpp_states_list_to_DEmat(areas_list=areas_list, states_list=states_list, dmat=dmat_times_d, elist=elist, amat=amat, include_null_range=BioGeoBEARS_run_object$include_null_range, normalize_TF=TRUE, makeCOO_TF=force_sparse)


# HERE1			
# 			ja_tmp = Qmat_tmp$ja
# 			Qmat_tmp$ja = Qmat_tmp$ia
# 			Qmat_tmp$ia = ja_tmp
			} else {
			Qmat_tmp = Qmat
			warning("WARNING in modify_Qmat_with_trait - a user-specified Qmat is being used, this will not be scaled by m / the trait_multiplier")
			} # END if (is.null(Qmat))

		if (force_sparse == FALSE)
			{		
			if (ncol(Qmat_tmp) != num_geog_states)
				{
				txt = cat("STOP ERROR in modify_Qmat_with_trait: the number of states in num_geog_states (=", num_geog_states, "; derived from ncol(tip_condlikes_of_data_on_each_state)/num_trait_states) and ncol(Qmat_tmp) (=", ncol(Qmat_tmp), ") don't match.", sep="")
				cat("\n\n")
				cat(txt)
				cat("\n\n")
				stop(txt)
				} # END if (ncol(Qmat_tmp) != num_geog_states)
			} # END if (force_sparse == FALSE)


		
		# Print dmat
		#print("dmat_times_d:")
		#print(dmat_times_d)
	
		#print("Qmat_tmp:")
		#print(Qmat_tmp)
		
		
		if (force_sparse == FALSE)
			{
			# DENSE MATRIX METHOD
			# Convert the Qmat to Rmat
			# No base frequencies, so put 1s
			Rmat_tmp = Qmat_to_Rmat(Qmat_tmp, basefreqs=rep(1,ncol(Qmat_tmp)), return_scalefactor=FALSE)
			Rmat_tmp[is.nan(Rmat_tmp)] = 0.0

			# Scale by total amount of change
			# Since the first row of the Rmat (null range)
			# may have a rate of 0 leaving (absorbing state)
			# The divide-by-sum operation can cause a NaN.
			# Fix this
			# (the scalefactor is a vector of off-diagonal rowSums)
			# (this is all kind of unnecessary unless we ever do basefreqs or tree scaling)
			scalefactor = Qmat_to_Rmat(Qmat_tmp, basefreqs=rep(1,ncol(Qmat_tmp)), return_scalefactor=TRUE)
			scalefactor[is.nan(scalefactor)] = 0.0
	
			# Multiply by scalefactor to get scaled 
			# rate of change
			# (Basically, the results is the Qmat without the negatives on the diagonal, if basefreqs are 1)
			Rmat_tmp = Rmat_tmp * scalefactor
			#print("Rmat_tmp:")
			#print(Rmat_tmp)


			# Fill in the blocks of the big Rmat
			for (traitstate_j in 1:ntrait_states)
				{
				# If it's the first time, make the big matrix
				if ((traitstate_i == 1) && (traitstate_j == 1))
					{
					Rmat = matrix(data=0, ncol=ncol(Qmat_tmp)*ntrait_states, nrow=nrow(Qmat_tmp)*ntrait_states)
					num_ranges = ncol(Qmat_tmp)
					} # END if ((traitstate_i == 1) && (traitstate_j == 1))


				# Rows; (current state)
				startrow = ((traitstate_i-1) * num_ranges) + 1
				endrow = (traitstate_i * num_ranges)
			
				# Columns (ending state)
				startcol = ((traitstate_j-1) * num_ranges) + 1
				endcol = (traitstate_j * num_ranges)
			
				if (traitstate_i == traitstate_j)
					{
					# If on the diagonal of the big 2x2 matrix,
					# geography can change but the 
					# trait stays the same
					Rmat_tmp[Rmat_tmp < 0] = 0
	# 						print(startrow)
	# 						print(endrow)
	# 						print(startcol)
	# 						print(endcol)
	# 						print(Rmat)
	# 						print(Rmat_tmp)
				
					# If there is only 1 geographic area, need a slightly different function
					if ((endrow-startrow) > 0)
						{
						Rmat[startrow:endrow,][,startcol:endcol] = Rmat_tmp
						}
					if ((endrow-startrow) == 0)
						{
						Rmat[startrow:endrow,][endcol] = Rmat_tmp
						}
										
				
					} else {
					# Off-diagonal of the big 2x2 matrix:
					# geography CANNOT change, but the 
					# trait CAN CHANGE
					#
					# This means JUST DIAGONALS
					# in cells 1,2 and 2,1
				
	#						print(ntrait_states)
	#						print(dim(trait_Pmat_txt))
	#						print(traitstate)
	#						print(traitstate_j)

	# 					print(traitstate_i)
	# 					print("traitstate_i")
	# 					print(traitstate_j)
	# 					print("traitstate_j")
					txt_for_trait_transition_rate = trait_Pmat_txt[traitstate_i, traitstate_j]
	# 						print("txt_for_trait_transition_rate:")
	# 						print(txt_for_trait_transition_rate)
					params_matrix_match_TF = row.names(BioGeoBEARS_model_object@params_table) == txt_for_trait_transition_rate
					trait_transition_rate = BioGeoBEARS_model_object@params_table$est[params_matrix_match_TF]

	# 				print("trait_transition_rate:")
	# 				print(trait_transition_rate)
				
					zeros_tmp = matrix(data=0, nrow=nrow(Qmat_tmp), ncol=ncol(Qmat_tmp))
					diag(zeros_tmp) = trait_transition_rate

					# If there is only 1 geographic area, need a slightly different function
					#print(zeros_tmp)
					if ((endrow-startrow) > 0)
						{
						Rmat[startrow:endrow,][,startcol:endcol] = zeros_tmp
						}
					if ((endrow-startrow) == 0)
						{
						Rmat[startrow:endrow,][endcol] = zeros_tmp
						}
					} # END if (traitstate == traitstate_j)
			
				# Populate Rmat
	#					print(startrow)
	#					print(endrow)
	#					print(startcol)
	#					print(endcol)
	#					print(Rmat)
	#					print(Rmat_tmp)
	#					print(trait_transition_rate)

					} # END for (traitstate_j in 1:ntrait_states)

				} else {
				# SPARSE MATRIX BUILDING OF LARGE RATE TRAIT-GEOGRAPHY MATRIX
		
				# Get the Rmat (zeros on diagonal, positive values on off-diagonal
				# COO-formatted

				# Convert the Qmat to Rmat, and add to the big trait-geog matrix
				# (so, the input Qmat_tmp has only num_geog_states)
				# No base frequencies, so put 1s
				cooRmat_df_nDiags = Qmat_to_Rmat_COO(cooQmat=Qmat_tmp, n=num_geog_states, basefreqs=rep(1,num_geog_states), return_scalefactor=FALSE)

				# Scale by total amount of change
				# Since the first row of the Rmat (null range)
				# may have a rate of 0 leaving (absorbing state)
				# The divide-by-sum operation can cause a NaN (in sparse, we exclude the diagonals)
				cooQmat=Qmat_tmp
				scalefactor = Qmat_to_Rmat_COO(cooQmat=Qmat_tmp, n=num_geog_states, basefreqs=rep(1,num_geog_states), return_scalefactor=TRUE)
				
				# Multiply by scalefactors (sums of each row) to get scaled 
				# rate of change
				# i.e. Rmat, rate matrix (no diagonals)
				cooRmat_df_nDiags$a = cooRmat_df_nDiags$a * scalefactor$a[cooRmat_df_nDiags$ia]
				sum(cooRmat_df_nDiags$a)
				
				
				# Fill in the blocks of the big Rmat
				num_ranges = num_geog_states
				for (traitstate_j in 1:ntrait_states)
					{
					# If it's the first time, make the big matrix
					if ((traitstate_i == 1) && (traitstate_j == 1))
						{
						# For COO, we can just append rows
						cooRmat_geogtrait_df = NULL
						} # END if ((traitstate_i == 1) && (traitstate_j == 1))

					# Rows; (current state)
					startrow = ((traitstate_i-1) * num_ranges) + 1
					endrow = (traitstate_i * num_ranges)
					newrow <- function(geog_row, traitstate_i, num_ranges)
						{
						trait_geog_row = ((traitstate_i-1) * num_ranges) + geog_row
						return(trait_geog_row)
						}
		
					# Columns (ending state)
					startcol = ((traitstate_j-1) * num_ranges) + 1
					endcol = (traitstate_j * num_ranges)
					newcol <- function(geog_col, traitstate_j, num_ranges)
						{
						trait_geog_col = ((traitstate_j-1) * num_ranges) + geog_col
						return(trait_geog_col)
						}
			
					# When the trait-states stay the same, the geography changes
					if (traitstate_i == traitstate_j)
						{
						# $a stays the same, but $ia and $ja might change to the bigger coordinates
						tmp_cooRmat_df = cooRmat_df_nDiags
						tmp_cooRmat_df$ia = newrow(geog_row=tmp_cooRmat_df$ia, traitstate_i=traitstate_i, num_ranges=num_ranges)
						tmp_cooRmat_df$ja = newcol(geog_col=tmp_cooRmat_df$ja, traitstate_j=traitstate_j, num_ranges=num_ranges)
				
						# Add to the COO-formatted Rmat
						cooRmat_geogtrait_df = rbind(cooRmat_geogtrait_df, tmp_cooRmat_df)
						} else {
						#### traitstate_i != traitstate_j ####
						# When the trait-states are different, the geography stays the same
						# Off-diagonal of the big 2x2 matrix:
						# geography CANNOT change, but the 
						# trait CAN CHANGE
						#
						# This means JUST DIAGONALS
						# in trait cells 1,2 and 2,1
				
						# Get the trait transition rate
						
						# FIX 2018-04-30 -- sparse matrices are assembled TRANSPOSED from dense matrices
						#txt_for_trait_transition_rate = trait_Pmat_txt[traitstate_i, traitstate_j]
						# Also have to sum the COLUMNS not the ROWS for negative diagonal
						txt_for_trait_transition_rate = trait_Pmat_txt[traitstate_j, traitstate_i]
						params_matrix_match_TF = row.names(BioGeoBEARS_model_object@params_table) == txt_for_trait_transition_rate
						trait_transition_rate = BioGeoBEARS_model_object@params_table$est[params_matrix_match_TF]
						
						# Build the Rmat matrix (COO version) for this quadrant
						tmp_ia = newrow(geog_row=1:num_ranges, traitstate_i=traitstate_i, num_ranges=num_ranges)
						tmp_ja = newcol(geog_col=1:num_ranges, traitstate_j=traitstate_j, num_ranges=num_ranges)
						tmp_a = rep(trait_transition_rate, times=num_ranges)
						tmp_cooRmat_df = as.data.frame(cbind(tmp_ia, tmp_ja, tmp_a))
						names(tmp_cooRmat_df) = c("ia", "ja", "a")
						
						# Add to the COO-formatted Rmat
						cooRmat_geogtrait_df = rbind(cooRmat_geogtrait_df, tmp_cooRmat_df)				
						} # END if (traitstate_i == traitstate_j)
					} # END for (traitstate_j in 1:ntrait_states)
				} # END if (force_sparse == FALSE)
			} # END for (traitstate_i in 1:ntrait_states)
		#print(Rmat)
	
	# Error check
# 	txt = paste0(cooRmat_geogtrait_df$ia, cooRmat_geogtrait_df$ja, sep="_")
# 	length(txt)
# 	length(unique(txt))

	if (force_sparse == FALSE)
		{
		# Convert the Rmat to Qmat
		# Print the Rmat, if desired
		#print("round(Rmat, 4):")
		#print(round(Rmat, 4))
		Qmat = Rmat_to_Qmat_noScaling(Rmat, basefreqs=rep(1,ncol(Rmat)) )
		} else {
		tmp_cooQmat_df3 = Rmat_to_Qmat_noScaling_COO(cooRmat_df=cooRmat_geogtrait_df, n=numstates_geogtrait, basefreqs=rep(1,numstates_geogtrait), use_scaling=FALSE, use_basefreqs=FALSE, scale_by_ja_instead_of_ia=TRUE)
		
		# Check: each row sums to 0
# 		aggregate(a~ia, data=tmp_cooQmat_df3, FUN=sum)
# 		tmp1 = tmp_cooQmat_df3
# 		tmp2 = tmp_cooQmat_df3
		
		# Sum diagonals and off-diagonals by row
# 		tmp1$a[tmp1$a<0] = 0
# 		tmp2$a[tmp1$a>0] = 0
# 		aggregate(a~ia, data=tmp1, FUN=sum)
# 		aggregate(a~ia, data=tmp2, FUN=sum)
		
# 		sum((aggregate(a~ia, data=tmp1, FUN=sum))$a)
# 		sum((aggregate(a~ia, data=tmp2, FUN=sum))$a)
		
		# Transpose here
		tmp_cooQmat_df4 = tmp_cooQmat_df3
		tmp_cooQmat_df4$ia = tmp_cooQmat_df3$ja
		tmp_cooQmat_df4$ja = tmp_cooQmat_df3$ia
		
		Qmat = tmp_cooQmat_df3
		} # END if (force_sparse == FALSE)
	
	res = NULL
	res$Qmat = Qmat
	res$m = m
	return(res)
	} # END modify_Qmat_with_trait <- function(BioGeoBEARS_run_object)











# Collapse the traits+geog states to traits-only, and geog-only, ancestral states results
# (everything else is the same as the input res)
traits_results_to_geogOnly_traitsOnly <- function(res=NULL, Rdata_fn=NULL)
	{
	if ( (is.null(res) == TRUE) &&  (is.null(Rdata_fn) == TRUE) )
		{
		stoptxt = "STOP ERROR in traits_results_to_geogOnly_traitsOnly(): inputs 'res' and 'Rdata_fn' cannot both be NULL"
		cat("\n\n")
		cat(stoptxt)
		cat("\n\n")
		stop(stoptxt)
		} # END if ( (is.null(res) == TRUE) &&  (is.null(Rdata_fn) == TRUE) )


	if (is.null(res) == TRUE)
		{
		# Force load to res
		res_name <- load(Rdata_fn)
	
		cmdtxt = paste0("res = ", res_name)
		eval(parse(text=cmdtxt))
		} # END if (is.null(res) == TRUE)


	# Input file
	Rdata_fn = "DEC+t12+t21+m2_inf.Rdata"

	# Output files
	traits_Rdata = gsub(pattern="\\.Rdata", replacement="_TRAITS.Rdata", x=Rdata_fn)
	geog_Rdata = gsub(pattern="\\.Rdata", replacement="_GEOG.Rdata", x=Rdata_fn)


	# Load original file
	# Loads to res
	load(Rdata_fn)
	res


	# Error check
	if ( is.null(res$inputs$trait) == TRUE )
		{
		stoptxt = paste0("STOP ERROR IN traits_results_to_geogOnly_traitsOnly(): Your input results object has nothing in res$inputs$trait -- it was probably not using the trait-based dispersal model. This function requires a trait-based dispersal model result in 'res'.")
	
		cat("\n\n")
		cat(stoptxt)
		cat("\n\n")
		stop(stoptxt)
		}

	# Copy to new objects
	traits_res = res
	geog_res = res





	# Extract the number of states, trait states, and geog states
	numnodes = nrow(res$ML_marginal_prob_each_state_at_branch_top_AT_node)
	num_allstates = ncol(res$ML_marginal_prob_each_state_at_branch_top_AT_node)
	numtraits = ncol(res$inputs$trait@df)
	num_geog_states = num_allstates / numtraits

	# Extract the number of geographic states

	# Error check
	if ( (num_geog_states == round(num_geog_states) ) == FALSE)
		{
		stoptxt = paste0("STOP ERROR IN traits_results_to_geogOnly_traitsOnly(): Your input results object from the trait-based dispersal model has to have the total number of states (here, ", num_allstates, ") has to be evenly dividable by ncol(res$inputs$trait@df) (here, ", numtraits, "). However, this is not the case here (num_geog_states=", num_geog_states, "), so the function has failed. Something is badly screwed up in your results object.")
	
		cat("\n\n")
		cat(stoptxt)
		cat("\n\n")
		stop(stoptxt)
		}





	#################################################
	# Collapse the traits
	#################################################
	startcols = seq(from=1, to=num_allstates, by=num_geog_states)
	endcols = seq(from=num_geog_states, to=num_allstates, by=num_geog_states)
	startcols
	endcols

	# Make a blank trait matrix
	blank_trait_matrix = matrix(data=0, ncol=numtraits, nrow=numnodes)
	computed_likelihoods_at_each_node = blank_trait_matrix
	relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS = blank_trait_matrix
	condlikes_of_each_state = blank_trait_matrix
	relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS = blank_trait_matrix
	relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS = blank_trait_matrix
	relative_probs_of_each_state_at_branch_top_AT_node_UPPASS = blank_trait_matrix
	ML_marginal_prob_each_state_at_branch_bottom_below_node = blank_trait_matrix
	ML_marginal_prob_each_state_at_branch_top_AT_node = blank_trait_matrix
	relative_probs_of_each_state_at_bottom_of_root_branch = rep(0, times=numtraits)

	# Sum in each trait category
	for (i in 1:length(startcols))
		{
		startcol = startcols[i]
		endcol = endcols[i]
	
		relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[,i] = rowSums(res$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[,startcol:endcol])
	condlikes_of_each_state[,i] = rowSums(res$condlikes_of_each_state[,startcol:endcol])
	relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[,i] = rowSums(res$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[,startcol:endcol])
	relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS[,i] = rowSums(res$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS[,startcol:endcol])
	relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[,i] = rowSums(res$relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[,startcol:endcol])
	ML_marginal_prob_each_state_at_branch_bottom_below_node[,i] = rowSums(res$ML_marginal_prob_each_state_at_branch_bottom_below_node[,startcol:endcol])
	ML_marginal_prob_each_state_at_branch_top_AT_node[,i] = rowSums(res$ML_marginal_prob_each_state_at_branch_top_AT_node[,startcol:endcol])
	relative_probs_of_each_state_at_bottom_of_root_branch[i] = sum(res$relative_probs_of_each_state_at_bottom_of_root_branch[startcol:endcol])
		}

	# Store the results
	traits_res$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS
	traits_res$condlikes_of_each_state = condlikes_of_each_state
	traits_res$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS
	traits_res$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS = relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS
	traits_res$relative_probs_of_each_state_at_branch_top_AT_node_UPPASS = relative_probs_of_each_state_at_branch_top_AT_node_UPPASS
	traits_res$ML_marginal_prob_each_state_at_branch_bottom_below_node = ML_marginal_prob_each_state_at_branch_bottom_below_node
	traits_res$ML_marginal_prob_each_state_at_branch_top_AT_node = ML_marginal_prob_each_state_at_branch_top_AT_node
	traits_res$relative_probs_of_each_state_at_bottom_of_root_branch = relative_probs_of_each_state_at_bottom_of_root_branch

	traits_res$ML_marginal_prob_each_state_at_branch_top_AT_node
	dim(traits_res$ML_marginal_prob_each_state_at_branch_top_AT_node)
	rowSums(traits_res$ML_marginal_prob_each_state_at_branch_top_AT_node)


	# Loads to "traits_res"
	txt = paste0("\nSaving results from trait-based dispersal model, ancestral trait estimates, to: '", traits_Rdata, "'...")
	cat(txt)
	save(traits_res, file=traits_Rdata)




	#################################################
	# Collapse the geographic ranges
	#################################################
	startcols = seq(from=1, to=num_allstates, by=num_geog_states)
	endcols = seq(from=num_geog_states, to=num_allstates, by=num_geog_states)
	startcols
	endcols

	# Make a blank geog matrix
	blank_geog_matrix = matrix(data=0, ncol=num_geog_states, nrow=numnodes)
	computed_likelihoods_at_each_node = blank_geog_matrix
	relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS = blank_geog_matrix
	condlikes_of_each_state = blank_geog_matrix
	relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS = blank_geog_matrix
	relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS = blank_geog_matrix
	relative_probs_of_each_state_at_branch_top_AT_node_UPPASS = blank_geog_matrix
	ML_marginal_prob_each_state_at_branch_bottom_below_node = blank_geog_matrix
	ML_marginal_prob_each_state_at_branch_top_AT_node = blank_geog_matrix
	relative_probs_of_each_state_at_bottom_of_root_branch = rep(0, times=num_geog_states)

	# Sum in each geog category
	for (i in 1:length(startcols))
		{
		startcol = startcols[i]
		endcol = endcols[i]
	
		relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[1:numnodes,] = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[1:numnodes,] + res$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[,startcol:endcol]

		condlikes_of_each_state[1:numnodes,] = condlikes_of_each_state[1:numnodes,] + res$condlikes_of_each_state[,startcol:endcol]

		relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[1:numnodes,] = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[1:numnodes,] + res$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[,startcol:endcol]

		relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS[1:numnodes,] = relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS[1:numnodes,] + res$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS[,startcol:endcol]

		relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[1:numnodes,] = relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[1:numnodes,] + res$relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[,startcol:endcol]

		ML_marginal_prob_each_state_at_branch_bottom_below_node[1:numnodes,] = ML_marginal_prob_each_state_at_branch_bottom_below_node[1:numnodes,] + res$ML_marginal_prob_each_state_at_branch_bottom_below_node[,startcol:endcol]

		ML_marginal_prob_each_state_at_branch_top_AT_node[1:numnodes,] = ML_marginal_prob_each_state_at_branch_top_AT_node[1:numnodes,] + res$ML_marginal_prob_each_state_at_branch_top_AT_node[,startcol:endcol]

		relative_probs_of_each_state_at_bottom_of_root_branch[1:num_geog_states] = relative_probs_of_each_state_at_bottom_of_root_branch[1:num_geog_states] + res$relative_probs_of_each_state_at_bottom_of_root_branch[startcol:endcol]
		}

	# Store the results
	geog_res$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS
	geog_res$condlikes_of_each_state = condlikes_of_each_state
	geog_res$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS
	geog_res$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS = relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS
	geog_res$relative_probs_of_each_state_at_branch_top_AT_node_UPPASS = relative_probs_of_each_state_at_branch_top_AT_node_UPPASS
	geog_res$ML_marginal_prob_each_state_at_branch_bottom_below_node = ML_marginal_prob_each_state_at_branch_bottom_below_node
	geog_res$ML_marginal_prob_each_state_at_branch_top_AT_node = ML_marginal_prob_each_state_at_branch_top_AT_node
	geog_res$relative_probs_of_each_state_at_bottom_of_root_branch = relative_probs_of_each_state_at_bottom_of_root_branch

	geog_res$ML_marginal_prob_each_state_at_branch_top_AT_node
	dim(geog_res$ML_marginal_prob_each_state_at_branch_top_AT_node)
	rowSums(geog_res$ML_marginal_prob_each_state_at_branch_top_AT_node)


	# Loads to "geog_res"
	txt = paste0("\nSaving results from trait-based dispersal model, ancestral geographic range estimates, to: '", geog_Rdata, "'...")
	cat(txt)
	save(geog_res, file=geog_Rdata)


	mods = NULL
	mods$geog_res = geog_res
	mods$traits_res = traits_res
	
	return(mods)
	} # END traits_results_to_geogOnly_traitsOnly <- function(res=NULL, Rdata_fn=NULL)










#######################################################
# Functions for re-running optimization
# BioGeoBEARS_rerun_optimization_v1.R
#######################################################
# 
# Given a BioGeoBEARS_results_object, modify the parameters and re-run
# 

extract_param_inference <- function(res)
	{
	BioGeoBEARS_run_object = res$inputs
	
	# Find the free parameters
	free_TF = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$type == "free"
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$est[free_TF]

	# Extract the original inference
	lnL = res$total_loglikelihood
	nparam = sum(free_TF)
	params = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$est[free_TF]
	tmp_rownames = row.names(BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table)[free_TF]
	orig_inf = c(lnL, nparam, params)
	names(orig_inf) = c("lnL", "nparam", tmp_rownames)
	orig_inf
	
	return(orig_inf)
	} # END extract_param_inference <- function(BioGeoBEARS_run_object)


# Re-run the optimization from low and high starting points
rerun_optimization_w_HiLow <- function(res=NULL, Rdata_fn=NULL, fraction_change=0.25, use_optimx=NULL, runslow=TRUE)
	{
	if ( (is.null(res)) && (is.null(Rdata_fn)) )
		{
		txt = paste0("STOP ERROR in rerun_optimization_w_HiLow(): inputs 'res' and 'Rdata_fn' cannot both be NULL.")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		} # END if ( (is.null(res)) && (is.null(Rdata_fn)) )
	
	if (is.null(Rdata_fn) == FALSE)
		{
		# Loads to something, probably "res"
		# Find out what it loaded to
		name_of_last_object <- load(Rdata_fn)

		print("name_of_last_object:")
		print(name_of_last_object)
		
		# Put it into "res"
		cmdtxt = paste0("res = ", name_of_last_object)
		eval(parse(text=cmdtxt))
		} # END if (is.null(Rdata_fn) == FALSE)
	
	# Save original res object
	orig_res = res
	
	
# 	print("orig_res$inputs$use_optimx:")
# 	print(orig_res$inputs$use_optimx)
	
	# Change the optimizer if desired
	if (is.null(use_optimx) == FALSE)
		{
		orig_res$inputs$use_optimx = use_optimx
		}
	

	# Extract the original inference
	orig_inf = extract_param_inference(orig_res)
	orig_inf

	# Get the run object from res
	BioGeoBEARS_run_object = orig_res$inputs

	# Load the tree
	tr = read.tree(BioGeoBEARS_run_object$trfn)
# 	phy = tr
	
# 	print(getwd())
# 	print(BioGeoBEARS_run_object$trfn)
# 	print(tr)
	
	# Find the free parameters
	free_TF = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$type == "free"
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$est[free_TF]


	# Re-run inference from previous ML parameters
	BioGeoBEARS_run_object = orig_res$inputs
	
	runslow = runslow
	resfn = gsub(pattern="\\.Rdata", replacement="_res_redo.Rdata", x=Rdata_fn)
	if (runslow)
			{
			res = bears_optim_run(BioGeoBEARS_run_object)
			res_redo = res

			save(res_redo, file=resfn)
			} else {
			# Loads to "res_redo"
			load(resfn)
			} # END if (runslow)

	# Extract parameters
	orig_redo = extract_param_inference(res=res_redo)



	# Given a res object, lower all parameters by fraction_change, filter for min, and re-estimate
	BioGeoBEARS_run_object = orig_res$inputs

	# Change the estimates by lowering by fraction_change
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$est[free_TF] = (1-fraction_change) * BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$est[free_TF]

	# Put the estimates into initial
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$init[free_TF] = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$est[free_TF]

	# Fix any parameters that are under the minimum
	BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
	check_BioGeoBEARS_run(BioGeoBEARS_run_object)

	# Update linked parameters
	BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
	BioGeoBEARS_run_object$BioGeoBEARS_model_object = calc_linked_params_BioGeoBEARS_model_object(BioGeoBEARS_model_object)

	# Re-run inference
	runslow = runslow
	resfn = gsub(pattern="\\.Rdata", replacement="_res_lowstart.Rdata", x=Rdata_fn)
	if (runslow)
			{
			res = bears_optim_run(BioGeoBEARS_run_object)
			res_lowstart = res

			save(res_lowstart, file=resfn)
			} else {
			# Loads to "res_lowstart"
			load(resfn)
			} # END if (runslow)

	# Extract parameters
	lowstart = extract_param_inference(res=res_lowstart)


	# Given a res object, multiply all parameters by (1+fraction_change), filter for max, and re-estimate
	# Get the run object from res
	BioGeoBEARS_run_object = orig_res$inputs

	# Find the free parameters
	free_TF = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$type == "free"
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$est[free_TF]

	# Change the estimates by multiplying by (1+fraction_change)
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$est[free_TF] = (1+fraction_change) * BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$est[free_TF]

	# Put the estimates into initial
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$init[free_TF] = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$est[free_TF]

	# Fix any parameters that are under the minimum
	BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
	check_BioGeoBEARS_run(BioGeoBEARS_run_object)

	# Update linked parameters
	BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
	BioGeoBEARS_run_object$BioGeoBEARS_model_object = calc_linked_params_BioGeoBEARS_model_object(BioGeoBEARS_model_object)

	# Re-run inference
	runslow = runslow
	resfn = gsub(pattern="\\.Rdata", replacement="_res_histart.Rdata", x=Rdata_fn)
	if (runslow)
			{
			res = bears_optim_run(BioGeoBEARS_run_object)
			res_histart = res

			save(res_histart, file=resfn)
			} else {
			# Loads to "res_histart"
			load(resfn)
			} # END if (runslow)

	# Extract parameters
	histart = extract_param_inference(res=res_histart)


	# Make table comparing inferences
	#rerun_optim_table_fn = "rerun_optim_table_v1.txt"
	rerun_optim_table_fn = gsub(pattern="\\.Rdata", replacement="_rerun_optim_table_v1.txt", x=Rdata_fn)
	rerun_optim_table = as.data.frame(rbind(orig_inf, orig_redo, lowstart, histart), stringsAsFactors=FALSE)
	rerun_optim_table = dfnums_to_numeric(rerun_optim_table)
	write.table(rerun_optim_table, file=rerun_optim_table_fn, sep="\t", quote=FALSE)

	cat("\n\nPrinting '", rerun_optim_table_fn, "':\n", sep="")
	print(rerun_optim_table)

	return(rerun_optim_table)
	} # END rerun_optimization_w_HiLow <- function(res=NULL, Rdata_fn=NULL)


