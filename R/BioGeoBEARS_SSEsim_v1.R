
#######################################################
# Set up the default inputs to SSEsim
#######################################################

#######################################################
# If trait_fn is NULL, skip traits (default)
# If trait_fn is a number, make a fake trait_fn with that many states
# Otherwise, load the given trait_fn
# (This is a lazy way to set up the traits model for simulation...but who cares?)
#######################################################


SSEsim_setup_inputs <- function(SSEmodel=NULL, BioGeoBEARS_run_object=NULL, trait_fn=NULL)
	{
	# This object holds the SSE sim inputs
	SSEsim_inputs = NULL
	
	# Set up the Birth Rate, Death Rate, etc.
	# These are the parameters for the tree, and for 
	# dependence of speciation/extinction on states
	if (is.null(SSEmodel))
		{
		SSEmodel = NULL
		
		# Birth rate (lambda)
		# This is the ML estimate under Yule on the Psychotria tree
		SSEmodel$brate = 0.3289132
		
		# Death rate (omega)
		SSEmodel$drate = SSEmodel$brate/3
		
		# Exponent on rangesize multiplier
		# Positive exponent means larger ranges have
		# increased speciation rates
		SSEmodel$rangesize_b_exponent = 1	

		# Exponent on rangesize multiplier
		# More negative exponents mean larger ranges have
		# decreases extinction rates
		SSEmodel$rangesize_d_exponent = -1

		# To get the values out:
		# 		brate = SSEmodel$brate
		# 		drate = SSEmodel$drate
		# 		rangesize_b_exponent = SSEmodel$rangesize_b_exponent
		# 		rangesize_d_exponent = SSEmodel$rangesize_d_exponent
		} # END if (is.null(SSEmodel))
	
	
	# Set up the BioGeoBEARS model (the model of range evolution)
	if (is.null(BioGeoBEARS_run_object))
		{
		# Use a default BioGeoBEARS_run_object, but specify some different parameters
		BioGeoBEARS_run_object = define_BioGeoBEARS_run()
		BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
		include_null_range = BioGeoBEARS_run_object$include_null_range

		#######################################################
		# Define the areas and states
		#######################################################
		# Get geographic ranges at tips
		tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=np(BioGeoBEARS_run_object$geogfn))

		# Get the list of geographic areas
		areas = getareas_from_tipranges_object(tipranges)
		areas_list = seq(0, length(areas)-1, 1)		# 0-base indexes
		areanames = areas
		areas
		areas_list

		max_range_size = length(areas)	# if Psychotria, this is 4
		BioGeoBEARS_run_object$max_range_size = max_range_size
		
		states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
		states_list

		state_indices_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=BioGeoBEARS_run_object$max_range_size, include_null_range=TRUE)
		state_indices_0based

		# Get the ranges
		ranges_list = areas_list_to_states_list_new(areas=areas, maxareas=length(areas),
		include_null_range=include_null_range, split_ABC=FALSE)
		ranges_list
		ranges = unlist(ranges_list)
		ranges
		

		###############################################
		# Default params for SSEsim
		###############################################
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","type"] = "free"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = 0.03
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = 0.03

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","type"] = "free"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = 0.03
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = 0.03

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.1
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.1

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["a","type"] = "fixed"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["a","init"] = 0
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["a","est"] = 0

		} else {
		# Read the input files, if any
		BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
		include_null_range = BioGeoBEARS_run_object$include_null_range

		#######################################################
		# Define the areas and states
		#######################################################
		# Get geographic ranges at tips
		tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=np(BioGeoBEARS_run_object$geogfn))

		# Get the list of geographic areas
		areas = getareas_from_tipranges_object(tipranges)
		areas_list = seq(0, length(areas)-1, 1)		# 0-base indexes
		areanames = areas
		areas
		areas_list

		max_range_size = BioGeoBEARS_run_object$max_range_size
		
		states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)
		states_list

		state_indices_0based = states_list

		# Get the ranges
		ranges_list = areas_list_to_states_list_new(areas=areas, maxareas=length(areas),
		include_null_range=include_null_range, split_ABC=FALSE)
		ranges_list
		ranges = unlist(ranges_list)
		ranges

		} # end is.null(BioGeoBEARS_run_object)
		




	#######################################################
	# If trait_fn is NULL, skip traits (default)
	# If trait_fn is a number, make a fake trait_fn with that many states
	# Otherwise, load the given trait_fn
	# (This is a lazy way to set up the traits model for simulation...but who cares?)
	#######################################################
	number_of_trait_states = NULL
	# There must be more than one trait state!
	if (!is.null(trait_fn) == TRUE)
		{
		if (is.numeric(trait_fn) == TRUE && (trait_fn > 1))
			{
			# Make a fake trait_fn:
			number_of_trait_states = trait_fn
			
			# Load tree and geography
			trfn = BioGeoBEARS_run_object$trfn
			tr = read.tree(trfn)
			geogfn = BioGeoBEARS_run_object$geogfn
			tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
			trait_values = tipranges
			traits_df = trait_values@df
			
			if (ncol(traits_df) <= number_of_trait_states)
				{
				traits_df = traits_df[,1:number_of_trait_states]
				} else {
				cols = rep(1, times=number_of_trait_states)
				traits_df = traits_df[, cols]
				
				# Edit for a 1-column problem
				if (length(cols) == 1)
					{
					traits_mat = as.matrix(x=traits_df, ncol=1)
					traits_df = as.data.frame(traits_mat, stringsAsFactors=FALSE)
					row.names(traits_df) = row.names(tipranges@df)
					}
				}
			
			# Randomly generate trait states at tips (doesn't matter, just no blanks)
			for (i in 1:nrow(traits_df))
				{
				statenum = sample(1:number_of_trait_states, size=1)
				traits_df[i,] = rep(0, times=number_of_trait_states)
				traits_df[i,][statenum] = 1
				}
			
			names(traits_df) = LETTERS[1:number_of_trait_states]
			trait_values@df = traits_df
			
			# Write to a traits file
			trait_fn = "random_traits_file.data"
			save_tipranges_to_LagrangePHYLIP(tipranges_object=trait_values, lgdata_fn=trait_fn)
			} # END if (is.numeric(trait_fn) == TRUE)
		
		# Load the traits file (given, or random)
		if (is.character(trait_fn) == TRUE)
			{
			# Load a given traits file
			trait_values = getranges_from_LagrangePHYLIP(lgdata_fn=trait_fn)

			# Add the traits data and model
			BioGeoBEARS_run_object = add_trait_to_BioGeoBEARS_run_object(BioGeoBEARS_run_object, trait_fn=trait_fn)
			
			if (is.null(SSEmodel$dej_params))
				{
				txt = paste0("STOP ERROR in SSEsim_setup_inputs(): For a traits analysis, SSEmodel must have SSEmodel$dej_params.")
				cat("\n\n")
				cat(txt)
				cat("\n\n")
				} # END if (is.null(SSEmodel$dej_params))
			
			# Input the parameter values for the traits model
			number_of_trait_states = ncol(trait_values@df)
			dej_params = SSEmodel$dej_params
			for (i in 1:number_of_trait_states)
				{
				# m parameters
				# e.g. m1 = SSEmodel$m1
				txt = paste0("BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['m", i, "', 'init'] = dej_params$m", i, "[param_iter]")
				eval(parse(text=txt))
				txt = paste0("BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['m", i, "', 'est'] = dej_params$m", i, "[param_iter]")
				eval(parse(text=txt))
			
				for (j in 1:number_of_trait_states)
					{
					# trait rate parameters			
					if (i != j)
						{
						# e.g. t12 = SSEmodel$t12
						txt = paste0("BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['t", i, j, "', 'init'] = dej_params$t", i, j, "[param_iter]")
						eval(parse(text=txt))
						txt = paste0("BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['t", i, j, "', 'est'] = dej_params$t", i, j, "[param_iter]")
						eval(parse(text=txt))
						}
					} # END for (j in 1:number_of_trait_states)
				} # END for (i in 1:number_of_trait_states)

			
			} # END if (is.character(trait_fn) == TRUE)
		} # END if (!is.null(trait_fn) == TRUE)
	# Store number of trait states
	SSEsim_inputs$number_of_trait_states = number_of_trait_states
		
	# Update linked parameters
	BioGeoBEARS_run_object$BioGeoBEARS_model_object = calc_linked_params_BioGeoBEARS_model_object(BioGeoBEARS_run_object$BioGeoBEARS_model_object)
	BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object

	
	#######################################################
	# Read in the SSE params
	#######################################################
	brate = SSEmodel$brate
	drate = SSEmodel$drate
	rangesize_b_exponent = SSEmodel$rangesize_b_exponent
	rangesize_d_exponent = SSEmodel$rangesize_d_exponent
	



	#######################################################
	# Use the anagenetic and cladogenetic parameters to make a Qmat and COOmat
	#######################################################	
	# Set the dispersal and extinction rate
	d = BioGeoBEARS_model_object@params_table["d","est"]
	e = BioGeoBEARS_model_object@params_table["e","est"]
	a = BioGeoBEARS_model_object@params_table["a","est"]
	
	# More inputs
	force_sparse = BioGeoBEARS_run_object$force_sparse
	areas = areas_list
	
	# Calculate the dispersal_multipliers_matrix
	dispersal_multipliers_matrix = dispersal_multipliers_matrix_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object)
	
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
	u = BioGeoBEARS_model_object@params_table["u","est"]
	extinction_modifier_list = area_of_areas ^ (1 * u)
	
	# Apply to extinction rate
	elist = extinction_modifier_list * rep(e, length(areas))
	
	# Set up the instantaneous rate matrix (Q matrix)
	Qmat = rcpp_states_list_to_DEmat(areas_list=areas_list, states_list=states_list, dmat=dmat_times_d, elist=elist, amat=amat, include_null_range=BioGeoBEARS_run_object$include_null_range, normalize_TF=TRUE, makeCOO_TF=force_sparse)



	# Is this a traits-based analysis?
	traitTF = is.null(BioGeoBEARS_run_object$trait) == FALSE

	# Analysis with a trait modifying dispersal rate
	if (traitTF == TRUE)
		{
		res = modify_Qmat_with_trait(Qmat=NULL, BioGeoBEARS_run_object=BioGeoBEARS_run_object, areas_list=areas_list, states_list=states_list, dispersal_multipliers_matrix=dispersal_multipliers_matrix, elist=elist, force_sparse=force_sparse)
		Qmat = res$Qmat
		m = res$m

		# If the trait can change during jump events
		jts_matrix = NULL
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
		m = NULL
		jts_matrix = NULL
		} # END if (if (traitTF == TRUE))




	#######################################################
	# Get the cladogenesis matrix
	#######################################################
	spPmat_inputs = spPmat_inputs_from_BioGeoBEARS_model_object(BioGeoBEARS_run_object, states_list, dispersal_multipliers_matrix)
	
	
	cppSpMethod=3
	printmat=FALSE
	m=m
	include_null_range=BioGeoBEARS_run_object$include_null_range
	jts_matrix=NULL

	# numstates_in_cladogenesis_matrix equals length(states_list), times length m, minus null range
	if (traitTF == FALSE)
		{
		if (include_null_range == TRUE)
			{
			numstates_in_cladogenesis_matrix = length(states_list) - 1
			}
		if (include_null_range == FALSE)
			{
			numstates_in_cladogenesis_matrix = length(states_list) - 0
			}
		} else {
		if (include_null_range == TRUE)
			{
			numstates_in_cladogenesis_matrix = (length(m) * length(states_list)) - 1
			}
		if (include_null_range == FALSE)
			{
			numstates_in_cladogenesis_matrix = (length(m) * length(states_list)) - 0
			}		
		} # END if (traitTF == FALSE)
	numstates_in_cladogenesis_matrix

	cppSpMethod=3
	numstates_in_cladogenesis_matrix=numstates_in_cladogenesis_matrix
	printmat=FALSE
	m=m
	include_null_range=BioGeoBEARS_run_object$include_null_range
	jts_matrix=NULL
	
	#print("m")
	#print(m)
	
	COOmat_Rsp_rowsums = spPmat_inputs_to_COO_weights_columnar(spPmat_inputs=spPmat_inputs, cppSpMethod=cppSpMethod, numstates_in_cladogenesis_matrix=numstates_in_cladogenesis_matrix, printmat=printmat, m=m, include_null_range=include_null_range, jts_matrix=jts_matrix)
	COO_weights_columnar = COOmat_Rsp_rowsums$COO_weights_columnar
	Rsp_rowsums = COOmat_Rsp_rowsums$Rsp_rowsums
	
	# Store the SSEsim inputs
	SSEsim_inputs$BioGeoBEARS_run_object = BioGeoBEARS_run_object
	SSEsim_inputs$Qmat = Qmat
	SSEsim_inputs$COO_weights_columnar = COO_weights_columnar
	SSEsim_inputs$Rsp_rowsums = Rsp_rowsums
	SSEsim_inputs$state_indices_0based = state_indices_0based
	SSEsim_inputs$ranges = ranges
	SSEsim_inputs$areanames = areanames
	SSEsim_inputs$SSEmodel = SSEmodel

	# To get them back out
# 	SSEmodel = SSEsim_inputs$SSEmodel
# 	Qmat = SSEsim_inputs$Qmat
# 	COO_weights_columnar = SSEsim_inputs$COO_weights_columnar
# 	Rsp_rowsums = SSEsim_inputs$Rsp_rowsums
#	state_indices_0based = SSEsim_inputs$state_indices_0based
#	ranges = SSEsim_inputs$ranges

	return(SSEsim_inputs)
	}



SSEsim_run <- function(SSEsim_inputs, rootstate=2, time_stop=1000, taxa_stop=1000, seedval=as.numeric(Sys.time()), printlevel=2, testwd="~")
	{
	defaults='
	rootstate = 2
	time_stop=100
	taxa_stop=50
	seed=1
	printlevel=4
	testwd=testwd
	'

	# Is this a traits-based analysis?
	BioGeoBEARS_run_object = SSEsim_inputs$BioGeoBEARS_run_object
	BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
	include_null_range = BioGeoBEARS_run_object$include_null_range
	traitTF = is.null(BioGeoBEARS_run_object$trait) == FALSE

	
	# Get the model parameters
	Qmat = SSEsim_inputs$Qmat
	COO_weights_columnar = SSEsim_inputs$COO_weights_columnar
	Rsp_rowsums = SSEsim_inputs$Rsp_rowsums
	
	#print("COO_weights_columnar")
	#print(COO_weights_columnar)
	#print("Rsp_rowsums")
	#print(Rsp_rowsums)
	
	state_indices_0based = SSEsim_inputs$state_indices_0based
	ranges = SSEsim_inputs$ranges
	SSEmodel = SSEsim_inputs$SSEmodel
	
	#numstates = length(state_indices_0based)
	numstates = dim(Qmat)[1]  # This is (the number of geog states) times (number of trait states), because dim(Qmat)

	#######################################################
	# Set up the parameters of a full forward simulation
	#######################################################
	brate = SSEmodel$brate
	drate = SSEmodel$drate
	rangesize_b_exponent = SSEmodel$rangesize_b_exponent
	rangesize_d_exponent = SSEmodel$rangesize_d_exponent

	# SSE rates for states
	# Set the birthrate to be a function of the number of areas in each range
	range_sizes = unlist(lapply(X=state_indices_0based, FUN=length))
	brates = brate * (range_sizes ^ rangesize_b_exponent)
	brates

	# Set the deathrate to be a function of the number of areas in each range
	drates = drate * (range_sizes ^ rangesize_d_exponent)
	drates
	
	# Modify brates and drates by the number of trait states
	if (traitTF == TRUE)
		{
		# Trait multipliers
		trait_multiplier_rows_TF = grepl(x=BioGeoBEARS_model_object@params_table$desc, pattern="trait-based dispersal rate multiplier")
		trait_multiplier_rows_indices = (1:length(trait_multiplier_rows_TF))[trait_multiplier_rows_TF]
		ntrait_states = length(trait_multiplier_rows_indices)
		BGB_trait_model_params_table = BioGeoBEARS_model_object@params_table[trait_multiplier_rows_indices,]
		# Set the m values (dispersal multipliers, dependent on trait)
		m = BGB_trait_model_params_table$est
		
		brates = rep(brates, times=length(m))
		drates = rep(drates, times=length(m))
		} else {
		m = NULL
		}
	
	
	# Get the rates from the Qmat
	# Source: https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/R/simulation.R?root=proteinevoutk20
	rates_de = -diag(Qmat) #exponential rate of waiting, diagonal of Qmat
	rates_de


	#######################################################
	# Run the forward simulation!
	#######################################################
	# Set the seed
	if (seedval > 2147483647)
		{
		seedval = seedval %% 2147483647
		}
	seed_input = set.seed(seedval)
	
	
	if (time_stop == 0 & taxa_stop == 0)
		{
		stop("Must have stopping criterion\n")
		}

	# Number of tries
	trynum = 0
	

	#######################################################
	# START 'WHILE' LOOP
	#######################################################
	# Run this while loop until you hit a successful simulation
	while (1)
		{
		trynum = trynum + 1
		
		if (printlevel >= 1)
			{
		cat("\n\n=================================\nTree simulation attempt #", trynum, "\n\n=================================\n\n", sep="")
			}

		###########################################################
		# Set up edge.ranges and begin tree		
		###########################################################
		# "edge" is a 2x2 matrix to start
		# rows are lineages, col#1 is parent node, col#2 is daughter node
		edge <- rbind(c(1, 2), c(1, 3))
		edge_length <- rep(NA, 2)
	
		# have another matrix giving zone 1 (tropics) or zone 2 (temperate)
		# start in the tropics for now
		# edge "state" is the state of the parent node
		edge.zone <- c(1, 1)
	
		# Simulate the first split
		# edge "range" contains the geographical range
		initial_state_1based = rootstate		# Start simulation in Kauai
	
		COO_probs_columnar = COO_weights_columnar
		daughter_states = given_a_starting_state_simulate_split(index_Qmat_0based_of_starting_state=(initial_state_1based-1), COO_probs_columnar=COO_probs_columnar, numstates=numstates)
	
		# Convert to 1-based
		daughter_states = daughter_states + 1
		
		edge_ranges <- rbind(daughter_states[1], daughter_states[2])
	
		# Create array to hold branch lengths
		stem_depth = numeric(2)
		alive_TF = rep(TRUE, 2)
		t = 0.0
		next_node = 4
		
		
		#######################################################
		# You need to store the root state and its descendants
		#######################################################
		# Set up tables to store events
		tmp_edgelength = 0
		repeatnum = 0
		de_event_num = 0
		cl_event_num = 0

		# Record anagenetic and cladogenetic events
		table_of_range_change_events = NULL
		table_of_cladogenetic_events = NULL
		
		# The root node has no edgenum, really, because there is no edge/branch below it
		edgenum_w_event = NA
		parent_range = initial_state_1based
		tmprow = c(edgenum_w_event, parent_range, daughter_states[1], daughter_states[2])

		# Add to the cladogenesis event table
		cmdstr = paste("cl", cl_event_num, " = tmprow", sep="")
		eval(parse(text=cmdstr))

		cmdstr = paste("table_of_cladogenetic_events = rbind(table_of_cladogenetic_events, cl", cl_event_num, ")", sep="")
		eval(parse(text=cmdstr))
		# Done storing event at root node
		#######################################################


		
		# A flag to stop if you hit a taxa_stop
		stopped_due_to_taxon_count = FALSE
	
		# Repeat until one of the break statements is reached
		repeat
			{
			repeatnum = repeatnum+1
			popsize = sum(alive_TF)

			if (printlevel >= 2)
				{
				cat("\n")
				cat("Treebuilding step #", repeatnum, ", # species living: ", popsize, "\n", sep="")
				#print(" ")
				cat(paste("Time: ", format(t, digits=3), " my, #alive = ", sum(alive_TF), "\n", sep=""))
				}
		
			# Density-dependence...
			#d = original_d + ((popsize/k) * (original_b-original_d))
			#print(paste("birth rate=b=", b, ", death rate=d=", d))
			#d = (events$top[14]-events$top[7])
			#b = 1-(d)
			#cat("birth events total: ", b, ", death rate=d=", d, "\n", sep="")

			# change things with 0 range to dead
	# 		for (i in 1:nrow(edge.ranges))
	# 			{
	# 			tmp_ranges = edge.ranges[i, ]
	# 			if (sum(tmp_ranges) == 0)
	# 				{
	# 				alive[i] = FALSE
	# 				}
	# 			}
		
		
			# Stop the simulation if everything is dead
			if (sum(alive_TF) == 0)
				{
				break
				}


			# Get waiting time to the next event (scaling factor * number alive); 
			# depends on number of lineages, perhaps times area
		
			# The rates depend on the states at each tip
			# Go through each tip
			current_ntips = length(edge_ranges[alive_TF])
			rates_of_events_per_tip = matrix(data=NA, nrow=current_ntips, ncol=3)
			for (j in 1:current_ntips)
				{
				rates_of_events_per_tip[j,1] = rates_de[edge_ranges[alive_TF][j]]
				rates_of_events_per_tip[j,2] = brates[edge_ranges[alive_TF][j]]
				rates_of_events_per_tip[j,3] = drates[edge_ranges[alive_TF][j]]
				}
		
		
			#rate = sf * get_rate(alive, edge_ranges, rate_calc="per_lineage")
			#print("rates_of_events_per_tip")
			#print(rates_of_events_per_tip)
			rate = sum(rates_of_events_per_tip)
			rate
		
			# dt is the amount of time until the next event on the tree
			# (the branches with non-events will be extended)
			# Higher rates result in shorter waiting times
			dt <- rexp(n=1, rate = rate)
		
			# This is the total time since start
			t <- t + dt
		
			# Stop if out of time before the next event
			time_stop_hit = FALSE
			if (time_stop)
				{
				#print(dt)
				#print(edge.zone[alive_TF]==1)
				#print(edge.zone[alive_TF]==2)
				if (t >= time_stop)
					{
					t <- time_stop
					time_stop_hit = TRUE
					break
					}
				} # END if (time_stop)

			# If you didn't hit the time barrier after the last speciation event,
			# stop the simulation if you've reached the # taxa_stop
			# You don't want to stop instantly, because then you will have
			# 2 zero-length branches
			# Instead, stop after (dt/1), i.e. when the next event happens
			# (without actually doing that next event; just extend the branches;
			# do this at the end of the while loop, for both time_stop and taxa_stop)
			if (taxa_stop)
				{
				#print(paste("Stop? #alive=", sum(alive_TF), "taxa_stop=", taxa_stop, sep=""))
				if (sum(alive_TF) >= taxa_stop)
					{
					break					
					}
				}


		
			# Choose which event happened
			event_probs = c(rates_of_events_per_tip) / rate
			event_probs

			# Which event happened?
			eventnum = sample(x=1:length(event_probs), size=1, replace=FALSE, prob=event_probs)
			eventnum
		
			eventnums = matrix(data=1:length(event_probs), nrow=current_ntips, ncol=3, byrow=FALSE)
			event_TF = eventnums == eventnum
			eventnums
			event_TF
		
			event_rownums = matrix(data=1:current_ntips, nrow=current_ntips, ncol=3, byrow=FALSE)
			event_colnums = matrix(data=1:3, nrow=current_ntips, ncol=3, byrow=TRUE)
			event_rownums
			event_colnums
		
			rownum = event_rownums[event_TF]
			colnum = event_colnums[event_TF]
		
			# Get the tip to change
			tip_to_change = rownum
		
			# Get event type
			event_type = colnum
		
			# List the edges to modify
			e <- matrix(edge[alive_TF, ], ncol = 2)
		
			edgenums_alive = (1:nrow(edge))[alive_TF]
			edgenums_alive
		
			# Which edge had the event?
			edgenum_w_event = edgenums_alive[tip_to_change]
		
			# Get the parent
			parent <- edge[edgenum_w_event, 2]
		
			#######################################################
			# Make the changes, based on the event
			#######################################################
		
			# Anagenetic range expansion/contraction event
			# TRAITS: OR, TRAIT CHARACTER CHANGE
			# (can result in extinction if ranges of size 1 drop to 0)
			if (event_type == 1)
				{
				if (printlevel >= 2)
					{
					cat("Event type: #1 (anagenetic range-change or trait-change)\n", sep="")
					}
			
				starting_range = edge_ranges[edgenum_w_event,]
			
				# Get the probabilities of new ranges
				# (zeroing out the diagonal, since we know the 
				# range doesn't stay the same)
				probs_of_new_ranges = Qmat[starting_range, ]
				probs_of_new_ranges[probs_of_new_ranges < 0] = 0
				probs_of_new_ranges = probs_of_new_ranges / sum(probs_of_new_ranges)
				
				#print("probs_of_new_ranges")
				#print(probs_of_new_ranges)
				new_rangenum_1based = sample(x=1:numstates, size=1, replace=FALSE, prob=probs_of_new_ranges)
				edge_ranges[edgenum_w_event, ] = new_rangenum_1based
			
				# Check if the range contraction led to extinction:
				if (include_null_range == FALSE)
					{
					range_contraction_to_zero = FALSE
					} else {
					if (traitTF == FALSE)
						{
						# There is a null range, and the range contracted to rangestate1, i.e. null range
						if (new_rangenum_1based == 1)
							{
							# Range contraction to 0 range
							alive_TF[edgenum_w_event] <- FALSE
							}
						} else {
						# Traits-including model, with null ranges
						num_trait_states = length(m)
						# The null geographic range will be every mth state
						mth_states = seq(from=1, to=numstates, by=(numstates/num_trait_states))
						if ((new_rangenum_1based %in% mth_states) == TRUE)
							{
							# Range contraction to 0 range
							alive_TF[edgenum_w_event] <- FALSE
							}
						} # END if (traitTF == FALSE)
					} # END if (include_null_range == FALSE)
			
				# Also keep track of this event
				de_event_num = de_event_num + 1
				tmprow = c(edgenum_w_event, t, dt, starting_range, new_rangenum_1based)
				cmdstr = paste("de", de_event_num, " = tmprow", sep="")
				eval(parse(text=cmdstr))
			
				if (printlevel >= 2)
					{
					print(tmprow)
					}
			
				# Add row to the table
				cmdstr = paste("table_of_range_change_events = rbind(table_of_range_change_events, de", de_event_num, ")", sep="")
				eval(parse(text=cmdstr))
				} # end anagenetic range-changing event
		
		
			# Speciation event
			# (requires sampling an cladogenetic range-changing event)
			if (event_type == 2)
				{
				if (printlevel >= 2)
					{
					cat("Event type: #2 (cladogenesis, including range-scenario)\n", sep="")
					}

				# Speciation event
				# take the lineage and temporarily "kill" it
				#parental_zone <- edge.zone[alive_TF][random_lineage]
				alive_TF[edgenum_w_event] <- FALSE

				# assign two new daughter nodes
				edge <- rbind(edge, c(parent, next_node), c(parent, next_node + 1))

				# move the "next node"
				next_node <- next_node + 2
			
				# Add two new living lineages
				alive_TF <- c(alive_TF, TRUE, TRUE)

				# Add the range of the parents to the daughters
				# (direct inheritance)
				#parent.range <- edge.ranges[edge_to_change,]
				# Copy the parent range to each daughter
				#edge.ranges <- rbind(edge.ranges, parent.range)
				#edge.ranges <- rbind(edge.ranges, parent.range)
			
				# Get the parent range
				parent_range = edge_ranges[edgenum_w_event,]
			
				# Sample a cladogenetic range-changing event
				#if (include_null_range == TRUE)
				#	{
					index_Qmat_0based_of_starting_state = (parent_range-1)
				#	} else {
#					index_Qmat_0based_of_starting_state = (parent_range-0)
				#	}
				
				
				daughter_states1 = given_a_starting_state_simulate_split(index_Qmat_0based_of_starting_state=index_Qmat_0based_of_starting_state, COO_probs_columnar=COO_probs_columnar, numstates=numstates)
				
				
				# Convert to 1-based
				daughter_states2 = daughter_states1 + 1
				
				
				daughter_states = daughter_states2
				
				# Add the two daughter ranges to edge_ranges:
				edge_ranges <- rbind(edge_ranges, daughter_states[1], daughter_states[2])
				
				# add to stem depth array
				stem_depth <- c(stem_depth, t, t)
			
				# add to edge length array
				#x <- which(edge[, 2] == parent)
				#edge_length[x] = t - stem_depth[x]
				edge_length = c(edge_length, NA, NA)


				# Write the event
				if (printlevel >= 2)
					{
					event_txt = paste(parent_range, " -> ", daughter_states[1], ", ", daughter_states[2], "\n", sep="")
					cat(event_txt)
					}
			
				# Add to the cladogenesis event table
				cl_event_num = cl_event_num + 1
				tmprow = c(edgenum_w_event, parent_range, daughter_states[1], daughter_states[2])
				cmdstr = paste("cl", cl_event_num, " = tmprow", sep="")
				eval(parse(text=cmdstr))
			
				cmdstr = paste("table_of_cladogenetic_events = rbind(table_of_cladogenetic_events, cl", cl_event_num, ")", sep="")
				eval(parse(text=cmdstr))
			
				} # end cladogenesis event

			# Extinction event (lineage-wide event; extinctions can also
			# occur through range-contraction)
			if (event_type == 3)
				{
				if (printlevel >= 2)
					{
					cat("Event type: #3 (lineage extinction)\n", sep="")
					}
				
				# Get the parent edge
				starting_range = edge_ranges[edgenum_w_event,]

				# Lineage extinction event
				# take the lineage and *permanently* "kill" it
				alive_TF[edgenum_w_event] <- FALSE
			
				# The range nevertheless gets copied up
				new_rangenum_1based = starting_range
				edge_ranges[edgenum_w_event, ] = new_rangenum_1based
			
				} # end extinction event
		
			# Update the edge lengths on everything that was alive during this step
			# update to edge length array

			# List the edges to modify
			#e <- matrix(edge[alive_TF, ], ncol = 2)
			x <- which(edge[, 2] == parent)
			edge_length[x] <- t - stem_depth[x]		
	
			} #end repeat

		# break unless you need to delete extinct lineages (?)
		#if (return.all.extinct == T | sum(alive_TF) > 1) 

		# Don't break if everything died!!
		if (sum(alive_TF) == 0)
			{
			pass_txt = "\nThis simulation failed i.e. died-out.\n"
			cat(pass_txt)
			} else {

			# Don't break if it was a time-stop hit!
			if (time_stop_hit == TRUE)
				{
				pass_txt = paste("\n\nSimulation trynum#", trynum, " failed as it didn't hit taxa_stop=", taxa_stop, " within ", time_stop, " my; trying again.\n", sep="")
				cat(pass_txt)
				} else {
				# Only break if you got a success!!
				# After you get a successful simulation, 
				# add the remaining amount of time until the next event
				edge_length[alive_TF] <- t - stem_depth[alive_TF]
			
				# You got a successful simulation, so exit
				success = TRUE
				break
				}
			}

		# If failure, print the output
		num_alive = sum(alive_TF)
		d = SSEsim_inputs$BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"]
		e = SSEsim_inputs$BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"]
		j = SSEsim_inputs$BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"]
		
		txtnums = adf2(matrix(data=c(trynum, num_alive, t, d, e, j, SSEmodel), nrow=1))
		names(txtnums) = c("try", "num_alive", "t", "d", "e", "j", "brate", "drate", "b_exp", "d_exp")
		print(txtnums)
		
		# Brian_crazy
 		#if (trynum > 100)
 		if (trynum > 5000)
 			{
 			simtr = NA
 			# (This would cause problems later on, though, let's skip for now)
 			success = FALSE
 			break
 			}
		} # end while loop
	#######################################################
	# END 'WHILE' LOOP
	#######################################################

	# raw edge matrix from simulation
	edges_from_sim = edge


	
	# Error check, i.e. if tree failed
	if (success == FALSE)
		{
		cat("\n\nSimulation of tree failed after many tries\n\n", sep="")
		# Convert the simulation into a real tree
		# Get the branch lengths
		edge_length[alive_TF] <- t - stem_depth[alive_TF]

		# Store the original edges
		original_edges = edge


		# change internal node numbers to (-1, -2, ... , -n)
		n = -1
		for (i in 1:max(edge))
			{
			if (any(edge[, 1] == i))
				{
				edge[which(edge[, 1] == i), 1] = n
				edge[which(edge[, 2] == i), 2] = n
				n = n - 1
				}
			}

		# re-number tips (everything that's left)
		edge[edge > 0] <- 1:sum(edge > 0)

		# Make an edge conversion table
		edge_conversion_table = cbind(original_edges, edge)
		edge_conversion_table = as.data.frame(edge_conversion_table, stringsAsFactors=FALSE)
		names(edge_conversion_table) = c("orig_ancnode", "orig_decnode", "ancnode", "decnode")
		edge_conversion_table

		# Create numbers for tip labels
		tip.label <- 1:sum(edge > 0)

		# change arrays "edge" and "tip.label" from numeric to character
		mode(edge) <- "character"
		mode(tip.label) <- "character"

		# Create the "phylo" object
		obj <- list(edge = edge, edge.length = edge_length, tip.label = tip.label)

		# call it a class phylo object
		class(obj) <- "phylo"

		# replace the old object with a valid APE object
		# This is the "raw" simulation, and the node numbers
		# are still screwed up compared to default APE
		obj2 <- old2new.phylo(obj)

		
		SSEsim_results = NULL
		SSEsim_results$states_list = SSEsim_inputs$state_indices_0based
		SSEsim_results$ranges_list = SSEsim_inputs$ranges
		SSEsim_results$areanames = SSEsim_inputs$areanames
		SSEsim_results$rootstate = rootstate
		SSEsim_results$rootnode = length(obj2$tip.label) + 1
		SSEsim_results$trynum = trynum
		SSEsim_results$simtr = obj2
		SSEsim_results$edge_conversion_table = NA
		SSEsim_results$edge_ranges = NA
		SSEsim_results$table_of_cladogenetic_events_translated = NA
		SSEsim_results$table_of_range_change_events_translated = NA
		SSEsim_results$success = FALSE
	
		return(SSEsim_results)
		} # end error check




	# Convert the table of range-changing (anagenetic) events to 
	# a data.frame
	# Error check, i.e. if no events
	if ((length(table_of_range_change_events) == 0) || (length(table_of_range_change_events) == 1))
		{
		# Fill with NAs
		table_of_range_change_events = adf2(matrix(data=NA, ncol=5))
		names(table_of_range_change_events) = c("edgenum_w_event", "t", "dt", "starting_range", "new_range")
		table_of_range_change_events_PROBLEM = TRUE
		} else {
		table_of_range_change_events = adf2(matrix(data=table_of_range_change_events, ncol=5))
		names(table_of_range_change_events) = c("edgenum_w_event", "t", "dt", "starting_range", "new_range")
		table_of_range_change_events		
		table_of_range_change_events_PROBLEM = FALSE
		}

	if ((length(table_of_cladogenetic_events) == 0) || (length(table_of_cladogenetic_events) == 1))
		{
		# Fill with NAs
		table_of_cladogenetic_events = NA
		table_of_cladogenetic_events_PROBLEM = TRUE
		} else {
		table_of_cladogenetic_events = data.frame(table_of_cladogenetic_events)
		names(table_of_cladogenetic_events) = c("edgenum_w_event", "parent_range", "Left_state", "Right_state")
		table_of_cladogenetic_events_PROBLEM = FALSE
		}









	# Convert the simulation into a real tree
	edge_ranges

	# Get the branch lengths
	edge_length[alive_TF] <- t - stem_depth[alive_TF]

	# Store the original edges
	original_edges = edge


	# change internal node numbers to (-1, -2, ... , -n)
	n = -1
	for (i in 1:max(edge))
		{
		if (any(edge[, 1] == i))
			{
			edge[which(edge[, 1] == i), 1] = n
			edge[which(edge[, 2] == i), 2] = n
			n = n - 1
			}
		}

	# re-number tips (everything that's left)
	edge[edge > 0] <- 1:sum(edge > 0)

	# Make an edge conversion table
	edge_conversion_table = cbind(original_edges, edge)
	edge_conversion_table = as.data.frame(edge_conversion_table)
	names(edge_conversion_table) = c("orig_ancnode", "orig_decnode", "ancnode", "decnode")
	edge_conversion_table

	# Create numbers for tip labels
	tip.label <- 1:sum(edge > 0)

	# change arrays "edge" and "tip.label" from numeric to character
	mode(edge) <- "character"
	mode(tip.label) <- "character"

	# Create the "phylo" object
	obj <- list(edge = edge, edge.length = edge_length, tip.label = tip.label)

	# call it a class phylo object
	class(obj) <- "phylo"

	# replace the old object with a valid APE object
	# This is the "raw" simulation, and the node numbers
	# are still screwed up compared to default APE
	obj2 <- old2new.phylo(obj)


	# Set a temporary directory for the newick file
	orig_wd = getwd()
	
	#testwd = "/drives/SkyDrive/_________thesis/_doc2/ch2_submission/2014-02-11_reviews/testsim/"
	setwd(testwd)

	# Convert the simulated tree node labels to default APE nodelabels
	# Write the tree out and read it back in
	trfn = "obj2.newick"
	write.tree(obj2, file=trfn)
	#obj4 = read.tree(file=trfn)
	obj4 = check_trfn(trfn=trfn)

	sim_tip_nodenums = as.numeric(obj4$tip.label)
	ape_tip_nodenums = 1:length(obj4$tip.label)
	tip_nodenums = cbind(sim_tip_nodenums, ape_tip_nodenums)
	tip_nodenums

	sim_int_nodenums = get_lagrange_nodenums(obj2)
	ape_int_nodenums = get_lagrange_nodenums(obj4)
	int_nodenums = cbind(sim_int_nodenums[,1], ape_int_nodenums[,1])
	int_nodenums
	int_nodenums = int_nodenums[order(int_nodenums[,2]), ]
	int_nodenums

	translation_table = rbind(tip_nodenums, int_nodenums)
	translation_table = adf2(translation_table)
	names(translation_table) = c("sim_nodenums", "ape_nodenums")
	translation_table



	# The edge conversion table has to convert negative tipnums
	# to obj2 nodenums, this is numtips-(negative nodenum)
	# (see old2new.phylo)
	ntips = length(obj2$tip.label)
	TF = edge_conversion_table < 0
	edge_conversion_table2 = edge_conversion_table
	edge_conversion_table2[TF] = ntips - as.numeric(edge_conversion_table[TF])

	edge_conversion_table = cbind(edge_conversion_table, edge_conversion_table2[,3:4], edge_conversion_table2[,3:4])
	names(edge_conversion_table) = c("orig_ancnode", "orig_decnode", "ancnode", "decnode", "sim_ancnode", "sim_decnode", "ape_ancnode", "ape_decnode")

	edge_conversion_table


	# Translate obj2 (sim nodenums) into obj3 (ape nodenums)
	obj3 = obj2
	obj4_tiplabels = as.numeric(obj4$tip.label)
	
	# Change the node numbers to match default APE node numbers
	for (i in 1:nrow(translation_table))
		{
		# Translate nodes in edges
		TF = obj2$edge == translation_table$sim_nodenums[i]
		obj3$edge[TF] = translation_table$ape_nodenums[i]

		# Translate tipnames
		TF = obj4_tiplabels == translation_table$sim_nodenums[i]
		TF
		obj3$tip.label[TF] = translation_table$ape_nodenums[i]
		
		# Translate edge_conversion_table
		TF = edge_conversion_table$sim_ancnode == translation_table$sim_nodenums[i]
		edge_conversion_table$ape_ancnode[TF] = translation_table$ape_nodenums[i]

		TF = edge_conversion_table$sim_decnode == translation_table$sim_nodenums[i]
		edge_conversion_table$ape_decnode[TF] = translation_table$ape_nodenums[i]
		
		# Translate raw sim node numbers into 
		#TF = edge_conversion_table$sim_ancnode == translation_table$sim_nodenums[i]
		#raw_simnode = unique(edge_conversion_table$orig_ancnode[TF])
		#simnode = unique(edge_conversion_table$sim_ancnode[TF])
		}



	# Plot the raw simtree, and the APE simtree
	if (printlevel >= 4)
		{
		plot(obj2, label.offset=0.2)
		axisPhylo()
		title("Simulation tree: raw node labels")
		tiplabels()
		nodelabels()

		plot(obj3, label.offset=0.2)
		axisPhylo()
		title("Simulation tree: APE node labels")
		tiplabels()
		nodelabels()

		cbind(obj2$edge, obj3$edge)
		}

	# Relabel the tips in obj3 so it says "sp1", "sp2", etc.
	obj3$tip.label = paste("sp", obj3$tip.label, sep="")
	

	
	
	# Translate the simulation edgenums_w_events into APE anc and decscendent nodes
	ape_ancnode = obj3$edge[table_of_cladogenetic_events$edgenum_w_event,1]
	ape_decnode = obj3$edge[table_of_cladogenetic_events$edgenum_w_event,2]
	table_of_cladogenetic_events_translated = cbind(ape_ancnode, ape_decnode, table_of_cladogenetic_events)
	names(table_of_cladogenetic_events_translated)[3] = "sim_edgenum_w_event"
	table_of_cladogenetic_events_translated
	
	# 2014-05-29_NJM fixing so that the root node is accurately recorded:
	# Fix the first row (representing the root, so no ancnode, but the decnode
	# number is the root number
	table_of_cladogenetic_events_translated$ape_decnode[1] = table_of_cladogenetic_events_translated$ape_ancnode[2]
	table_of_cladogenetic_events_translated$sim_edgenum_w_event[1] = 1
	
	if (table_of_range_change_events_PROBLEM == FALSE)
		{
		table_of_range_change_events_translated = table_of_range_change_events
		ape_ancnode = obj3$edge[table_of_range_change_events$edgenum_w_event,1]
		ape_decnode = obj3$edge[table_of_range_change_events$edgenum_w_event,2]
		table_of_range_change_events_translated = cbind(ape_ancnode, ape_decnode, table_of_range_change_events)
		names(table_of_range_change_events_translated)[3] = "sim_edgenum_w_event"
		table_of_range_change_events_translated
		} else {
		table_of_range_change_events_translated = table_of_range_change_events
		}

	
	#######################################################
	# Store the simulation results, tables, etc.
	#######################################################
	SSEsim_results = NULL
	
	# Store the inputs, because you always should
	SSEsim_results$SSEsim_inputs = SSEsim_inputs
	
	# We will need the states_list and ranges_list later
	SSEsim_results$states_list = SSEsim_inputs$state_indices_0based
	SSEsim_results$ranges_list = SSEsim_inputs$ranges
	SSEsim_results$areanames = SSEsim_inputs$areanames

	SSEsim_results$rootstate = rootstate
	SSEsim_results$rootnode = length(obj3$tip.label) + 1
	
	# Number of tries to get a successful simulation
	SSEsim_results$trynum = trynum

	# Simulated tree with APE node numbers
	# (should save and read the same, thankfully)
	SSEsim_results$simtr = obj3

	# The ranges at the ends of the branches/edges, from the 
	# original simulation.
	# The edges remain in the original order
	# See the edge conversion table; the nodes holding these edge
	# ranges are edge_conversion_table$sim_decnode
	SSEsim_results$edge_conversion_table = edge_conversion_table
	SSEsim_results$edge_ranges = edge_ranges
	SSEsim_results$edge_length = edge_length
	
	SSEsim_results$table_of_cladogenetic_events_translated = table_of_cladogenetic_events_translated
	SSEsim_results$table_of_range_change_events_translated = table_of_range_change_events_translated
	
	
	SSEsim_results$success = success
	
	# Double-check the output
	# cbind(obj3$edge, obj2$edge, obj3$edge.length, obj2$edge.length, edge_conversion_table, edge_ranges)
	
	# Return to the original working directory
	setwd(orig_wd)
	
	return(SSEsim_results)
	}


get_simtr_full_to_observed_node_translation <- function(simtr_complete, simtr_observed, simtr_observed_table)
	{
	# Build a table to translate between the full simulation tree nodes, and the 
	fulltr_internal_nodenums = as.numeric(gsub(pattern="fulltr_node", replacement="", x=simtr_observed$node.label))
	fulltr_internal_nodenums
	
	# Get the species names of the smaller observed tree
	fulltr_tip_labels = simtr_observed$tip.label
	fulltr_tip_labels
	
	match_positions_of_first_in_second = match(x=fulltr_tip_labels, table=simtr_complete$tip.label)
	match_positions_of_first_in_second
	
	fulltr_nodenums = c(match_positions_of_first_in_second, fulltr_internal_nodenums)
	obstr_nodenums = simtr_observed_table$node

	simtr_full_to_observed_node_translation = cbind(fulltr_nodenums, obstr_nodenums)
	simtr_full_to_observed_node_translation = adf2(simtr_full_to_observed_node_translation)
	names(simtr_full_to_observed_node_translation) = c("full", "obs")
	simtr_full_to_observed_node_translation
	
	return(simtr_full_to_observed_node_translation)
	} # END get_simtr_full_to_observed_node_translation





# Take the SSEsim results and convert them to:
# 1. Newick and geography files WITH extinct tips
# 2. Newick and geography files WITHOUT extinct tips
# 3. True state probabilities for each node in complete tree 
# 4. True state probabilities for each node in observed tree
# 5. DE event counts (true and observed)
# 6. Cladogenetic event counts
# 7. Cladogenetic event counts preserved in observed tree
SSEsim_to_files <- function(SSEsim_results, simdir, fossils_older_than=0.001, printlevel=2, write_files=TRUE)
	{
	defaults='
	simdir = "/drives/SkyDrive/_________thesis/_doc2/ch2_submission/2014-02-11_reviews/testsim/"
	fossils_older_than=0.001
	printlevel=4
	write_files=TRUE
	'

	# Working directories
	orig_wd = getwd()
	setwd(simdir)
	
	simtr = SSEsim_results$simtr
	states_list = SSEsim_results$states_list
	ranges_list = SSEsim_results$ranges_list
	areanames = SSEsim_results$areanames
	rootstate = SSEsim_results$rootstate
	rootnode = SSEsim_results$rootnode
	
	BioGeoBEARS_run_object = SSEsim_results$SSEsim_inputs$BioGeoBEARS_run_object
	BioGeoBEARS_model_object = SSEsim_results$SSEsim_inputs$BioGeoBEARS_run_object$BioGeoBEARS_model_object
	states_list = SSEsim_results$states_list

	
	# Is there a trait?
	trait_TF = !is.null(SSEsim_results$SSEsim_inputs$BioGeoBEARS_run_object$trait)
	if (trait_TF == TRUE)
		{
		# Trait multipliers
		trait_multiplier_rows_TF = grepl(x=BioGeoBEARS_model_object@params_table$desc, pattern="trait-based dispersal rate multiplier")
		trait_multiplier_rows_indices = (1:length(trait_multiplier_rows_TF))[trait_multiplier_rows_TF]
		ntrait_states = length(trait_multiplier_rows_indices)
		BGB_trait_model_params_table = BioGeoBEARS_model_object@params_table[trait_multiplier_rows_indices,]
	
		# Get the m values (dispersal multipliers, dependent on trait)
		m = BGB_trait_model_params_table$est
		num_trait_states = length(m)
		
		# Numbering of geog vs. trait states
		num_geog_states = length(SSEsim_results$states_list)
		
		trait_state_indices_list_1based = NULL
		max_index_1based = 0
		for (i in 1:length(m))
			{
			tmp_indices_list = (max_index_1based+1) : (max_index_1based + num_geog_states)
			max_index_1based = max(tmp_indices_list)
			trait_state_indices_list_1based[[i]] = tmp_indices_list
			}
		
		numstates = length(states_list) * length(m)
		
		} else {
		num_geog_states = length(SSEsim_results$states_list)
		m = NULL
		trait_state_indices_list_1based = NULL
		numstates = length(states_list)
		} # END if (trait_TF == TRUE)
	
	
	
	# Label the original nodes of the full tree, so we can link
	# them to the observed tree
	ntips = length(simtr$tip.label)
	tipnums = 1:ntips
	nodenums = (ntips+1):(ntips+simtr$Nnode)
	nodenums
	all_nodenums = c(tipnums, nodenums)
	numnodes = length(all_nodenums)
	
	# PUT IN NODE LABELS SO THAT YOU CAN IDENTIFY NODES ON THE OBSERVED TREE LATER
	simtr$node.label = paste("fulltr_node", nodenums, sep="")
	
	# Get the tips to drop
	simtr_table = prt(simtr, fossils_older_than=fossils_older_than, printflag=FALSE, get_tipnames=TRUE)
	tips_to_drop_TF = simtr_table$fossils[1:ntips]
	num_fossil_tips = sum(tips_to_drop_TF)

	if (num_fossil_tips > 0)
		{
		tips_to_drop = simtr$tip.label[tips_to_drop_TF]
		simtr_observed = ape::drop.tip(phy=simtr, tip=tips_to_drop, trim.internal=TRUE)
		#simtr_observed$node.label = (length(simtr_observed$tip.label)+1):(length(simtr_observed$tip.label)+simtr_observed$Nnode)
		simtr_observed_table = prt(simtr_observed, fossils_older_than=fossils_older_than, printflag=FALSE, get_tipnames=TRUE)
		
		# Build a table to translate between the full simulation tree nodes, and the 
		fulltr_internal_nodenums = as.numeric(gsub(pattern="fulltr_node", replacement="", x=simtr_observed$node.label))
		fulltr_internal_nodenums
		
		# Get the species names of the smaller observed tree
		fulltr_tip_labels = simtr_observed$tip.label
		fulltr_tip_labels
		
		match_positions_of_first_in_second = match(x=fulltr_tip_labels, table=simtr$tip.label)
		match_positions_of_first_in_second
		
		fulltr_nodenums = c(match_positions_of_first_in_second, fulltr_internal_nodenums)
		obstr_nodenums = simtr_observed_table$node

		simtr_full_to_observed_node_translation = cbind(fulltr_nodenums, obstr_nodenums)
		simtr_full_to_observed_node_translation = adf2(simtr_full_to_observed_node_translation)
		names(simtr_full_to_observed_node_translation) = c("full", "obs")
		simtr_full_to_observed_node_translation

		} else {
		# No tips to drop
		simtr_observed = simtr
		simtr_observed_table = simtr_table
		
		simtr_full_to_observed_node_translation = cbind(simtr_table$node, simtr_table$node)
		simtr_full_to_observed_node_translation = adf2(simtr_full_to_observed_node_translation)
		names(simtr_full_to_observed_node_translation) = c("full", "obs")
		}
	
	# Calculate the difference in root heights, to later fix t in the 
	# observed tree
	difference_in_rootheight_between_simulated_and_observed = get_max_height_tree(simtr) -  get_max_height_tree(simtr_observed)
	
	# Plot the tree
	if (printlevel >= 4)
		{
		plot(simtr_observed, label.offset=0.2)
		axisPhylo()
		title("Simulated tree (observed)")
		nodelabels()
		tiplabels()
		}

	# Get the nodenumbers of the original full tree, which correspond to the
	# ordered node numbers of the subset tree
	list_of_full_nodelabels_in_subset_tree = c(simtr_observed$tip.label, simtr_observed$node.label)
	list_of_full_nodenums_in_subset_tree = gsub(pattern="sp", replacement="", x=list_of_full_nodelabels_in_subset_tree)
	list_of_full_nodenums_in_subset_tree = gsub(pattern="fulltr_node", replacement="", x=list_of_full_nodenums_in_subset_tree)
	list_of_full_nodenums_in_subset_tree = as.numeric(list_of_full_nodenums_in_subset_tree)
	list_of_full_nodenums_in_subset_tree
	
	list_of_full_tipnums_in_subset_tree = list_of_full_nodenums_in_subset_tree[1:length(simtr_observed$tip.label)]
	list_of_full_tipnums_in_subset_tree
	
	# Label the cladogenetic events
	table_of_cladogenetic_events_translated = SSEsim_results$table_of_cladogenetic_events_translated
	
	track_dispersal_dest = TRUE
	cladogenetic_event_labels = label_table_of_cladogenetic_events(table_of_cladogenetic_events_translated, states_list, ranges_list, areanames=areanames, track_dispersal_dest=track_dispersal_dest, m=m)
	table_of_cladogenetic_events_translated = cbind(table_of_cladogenetic_events_translated, cladogenetic_event_labels)
	table_of_cladogenetic_events_translated
	
	# Label the anagenetic events
	table_of_range_change_events_translated = SSEsim_results$table_of_range_change_events_translated
	
	# Test for the range change events table being NA
	# If it's not, add the anagenetic event labels
	naTF1 = ( (length(table_of_range_change_events_translated) == 1) && (is.na(table_of_range_change_events_translated[1])) )
	# Just check everything for NA
	naTF2 = all(is.na(table_of_range_change_events_translated))
	naTF = ((naTF1 + naTF2) >= 1)
	if (naTF == FALSE)
		{
		anagenetic_event_labels = label_table_of_anagenetic_events(table_of_range_change_events_translated, states_list, ranges_list, areanames=areanames, m=m)
	
		table_of_range_change_events_translated = cbind(table_of_range_change_events_translated, anagenetic_event_labels)
		}

	#######################################################
	# Make tipranges object for complete, and observed trees
	#######################################################
	edge_ranges = SSEsim_results$edge_ranges

	# Edit edge ranges if traits analysis
	if (trait_TF == TRUE)
		{
		edge_traits = rep(NA, times=length(edge_ranges))
		for (i in 1:length(edge_ranges))
			{
			# Translate the edge range statenums
			TFs = rep(FALSE, times=num_trait_states)
			for (j in 1:num_trait_states)
				{
				TF = edge_ranges[i] %in% trait_state_indices_list_1based[[j]]
				TFs[j] = TF
				}
			edge_trait_statenum = (1:num_trait_states)[TFs]
			edge_geog_state_tmp = (edge_ranges[i] - ((edge_trait_statenum-1)*num_geog_states) )
			edge_ranges[i] = edge_geog_state_tmp
			edge_traits[i] = edge_trait_statenum
			} # END for (i in 1:length(edge_ranges))

		edges_traits_table = cbind(SSEsim_results$edge_conversion_table, edge_traits)
		} # END if (trait_TF == TRUE)
	
	
	edges_table = cbind(SSEsim_results$edge_conversion_table, edge_ranges)
	edges_table2 = edges_table[order(edges_table$ape_decnode), ]
	edges_table2
	
	tiplabels_for_tipranges = simtr$tip.label
	tipranges_simulated = edges_table_to_tipranges_object(edges_table2, states_list, ntips, tiplabels_for_tipranges, areanames, addval=0)	

	# Store the traits
	if (trait_TF == TRUE)
		{
		edges_traits_table2 = edges_traits_table[order(edges_table$ape_decnode), ]
		edges_traits_table2$edge_ranges = edges_traits_table2$edge_traits

		tiplabels_for_tipranges = simtr$tip.label
		traits_states_list = as.list( 0:(length(m)-1) )
		traitnames = 1:length(m)
		traits_simulated = edges_table_to_tipranges_object(edges_table2=edges_traits_table2, states_list=traits_states_list, ntips=ntips, tiplabels_for_tipranges=tiplabels_for_tipranges, areanames=traitnames, addval=0)	
		} # END if (trait_TF == TRUE)





	# Get the (TRUE!) probabilities of the ancestral states
	# on the complete tree
	ancstate_probs_simfull = matrix(data=0, nrow=numnodes, ncol=numstates)
	for (i in 1:numnodes)
		{
		tmpstate_index = edges_table2$edge_ranges[i]
		tmpnodenum_index = edges_table2$ape_decnode[i]
		ancstate_probs_simfull[tmpnodenum_index, tmpstate_index] = 1
		}
	# Put in the root state
	ancstate_probs_simfull[rootnode, rootstate] = 1
	ancstate_probs_simfull
	rowSums(ancstate_probs_simfull)
	

	
	# Reduce the event lists based on what went extinct
	internal_nodes_that_are_observed_TF = simtr$node.label %in% simtr_observed$node.label
	internal_nodes_that_are_observed_TF
	internal_nodes_that_went_extinct_TF = internal_nodes_that_are_observed_TF == FALSE
	if (sum(internal_nodes_that_went_extinct_TF) > 0)
		{
		internal_nodes_that_went_extinct_labels = simtr$node.label[internal_nodes_that_went_extinct_TF]
		internal_nodes_that_went_extinct_labels
		internal_nodenums_that_went_extinct = as.numeric(gsub(pattern="fulltr_node", replacement="", x=internal_nodes_that_went_extinct_labels))
		
		tipnums_that_went_extinct = as.numeric(gsub(pattern="sp", replacement="", x=tips_to_drop))
		
		internal_nodenums_that_went_extinct
		tipnums_that_went_extinct
		
		nodenums_extinct = sort(c(tipnums_that_went_extinct, internal_nodenums_that_went_extinct))
		
		# Now, subset the de_events table and the cladogenetic events table
		
		# cladogenetic
		TF = table_of_cladogenetic_events_translated$ape_decnode %in% nodenums_extinct == FALSE
		table_of_cladogenetic_events_observed = table_of_cladogenetic_events_translated[TF, ]
		table_of_cladogenetic_events_observed
		
		
		
		
			#######################################################
			# 
			# 2014-06-01_NJM:
			# Convert the missing cladogenetic events to anagenetic events
			# 
			#######################################################
			if (nrow(table_of_cladogenetic_events_observed) > 0)
			#if (sum(TF == FALSE) > 0)
				{
				table_of_cladogenetic_events_extinct = table_of_cladogenetic_events_translated[TF == FALSE,]
				table_of_cladogenetic_events_extinct
		
				# First, check for cladogenesis events where both daughters go extinct
				node_completely_extinct_TF = rep(NA, nrow(table_of_cladogenetic_events_extinct))
				for (m in 1:length(table_of_cladogenetic_events_extinct$ape_decnode))
					{
					tmp_nodenum = table_of_cladogenetic_events_extinct$ape_decnode[m]
					# Get the tips 
					tmp_tipnames = get_all_daughter_tips_of_a_node(nodenum=tmp_nodenum, t=simtr)
			
					# Check if any of the tips are in the observed tree
					tmp_tipnames_in_observed_tree_TF = tmp_tipnames %in% simtr_observed$tip.label
					if (sum(tmp_tipnames_in_observed_tree_TF) == 0)
						{
						# Then this node is has no sampled descendants
						node_completely_extinct_TF[m] = TRUE
						} else {
						node_completely_extinct_TF[m] = FALSE
						}
					} # END for (m in 1:length(table_of_cladogenetic_events_extinct$ape_decnode))
				table_of_cladogenetic_events_w_1_descendant = table_of_cladogenetic_events_extinct[node_completely_extinct_TF == FALSE,]
				dim(table_of_cladogenetic_events_translated)
				dim(table_of_cladogenetic_events_extinct)
				dim(table_of_cladogenetic_events_w_1_descendant)


				######################################################
				# Extract anagenetic events on observed tree
				######################################################
				# 2014-06-03_NJM:
				# Get the anagenetic events on the nodes corresponding to the observed tree,
				# AND on the branches below that which each had cladogenesis events
				# with 1 descendant
			
				# Just living/observed nodes
				TF = table_of_range_change_events_translated$ape_decnode %in% nodenums_extinct == FALSE
				table_of_range_change_events_on_branches_below_nodes_preserved_in_observed_tree = table_of_range_change_events_translated[TF, ]
			
			
				# Nodes in full tree that have ONE descendant in the observed tree
				fulltr_nodenums_w_1_desc = table_of_cladogenetic_events_w_1_descendant$ape_decnode
			
				# Those nodes in the fulltr anagenetic table
				TF = table_of_range_change_events_translated$ape_decnode %in% fulltr_nodenums_w_1_desc
			
			
				table_of_range_change_events_on_branches_below_nodes_w1desc_in_observed_tree = table_of_range_change_events_translated[TF,]
			
				table_of_range_change_events_observed = rbind(table_of_range_change_events_on_branches_below_nodes_preserved_in_observed_tree, table_of_range_change_events_on_branches_below_nodes_w1desc_in_observed_tree)
				# 2014-06-03_NJM: THAT SHOULD FIX THE ISSUE WHERE EARLY d EVENTS IN OBSERVED
				#                 TREE ARE NOT PLOTTED!!!!!!!!!!!!!!!!!
			
				} else {
				# No extinction events causing unsampled nodes
				# (not really -- currently the test is for any observed
				#  clado events, which will always be positive, so this will
				#  never be reached
				table_of_cladogenetic_events_extinct = NA
				clado_events_converted_to_anagenetic = NA

				# anagenetic
				# 2014-06-03_NJM: this is the old/bad way of doing it,
				# just captures the d/e events on the branches below
				# the nodes preserved in the full tree
				TF = table_of_range_change_events_translated$ape_decnode %in% nodenums_extinct == FALSE
				table_of_range_change_events_observed = table_of_range_change_events_translated[TF, ]
				table_of_range_change_events_observed

				} # END if (sum(TF == FALSE) > 0)
		

		
		
		
		
		
		# Correct "t" (the absolute time above the full simulated tree root)
		# right here
		table_of_range_change_events_observed$t = table_of_range_change_events_observed$t - difference_in_rootheight_between_simulated_and_observed
		
		# subset ancestral probabilities
		rownums = 1:nrow(ancstate_probs_simfull)
		keepTF = (rownums %in% nodenums_extinct) == FALSE
		ancstate_probs_simobs = ancstate_probs_simfull[list_of_full_nodenums_in_subset_tree, ]
		
		# subset tipranges
		tipranges_observed = tipranges_simulated
		tipranges_observed@df = tipranges_simulated@df[list_of_full_tipnums_in_subset_tree, ]
		tipranges_observed
	
		if (trait_TF == TRUE)
			{
			traits_observed = traits_simulated
			# Subset to observed tips:
			traits_observed@df = traits_simulated@df[list_of_full_tipnums_in_subset_tree, ]
			} # END if (trait_TF == TRUE)
		
		
		
		} else {
		# Nothing went extinct, just use these
		table_of_cladogenetic_events_extinct = NA
		table_of_cladogenetic_events_w_1_descendant = NA
		table_of_cladogenetic_events_observed = table_of_cladogenetic_events_translated
		table_of_range_change_events_observed = table_of_range_change_events_translated
		ancstate_probs_simobs = ancstate_probs_simfull
		tipranges_observed = tipranges_simulated

		if (trait_TF == TRUE)
			{
			traits_observed = traits_simulated
			# Subset to observed tips:
			traits_observed@df = traits_simulated@df[list_of_full_tipnums_in_subset_tree, ]
			} # END if (trait_TF == TRUE)

		} # END if (sum(internal_nodes_that_went_extinct_TF) > 0)


	
	#######################################################
	# Get some summary statistics
	#######################################################
	simstats = NULL
	
	num_fossil_tips = sum(tips_to_drop_TF)
	num_extinct_nodes = simtr$Nnode - simtr_observed$Nnode

	y_actual = sum(table_of_cladogenetic_events_translated$event_type == "sympatry (y)")
	s_actual = sum(table_of_cladogenetic_events_translated$event_type == "subset (s)")
	v_actual = sum(table_of_cladogenetic_events_translated$event_type == "vicariance (v)")
	j_actual = sum(table_of_cladogenetic_events_translated$event_type == "founder (j)")

	y_observed = sum(table_of_cladogenetic_events_observed$event_type == "sympatry (y)")
	s_observed = sum(table_of_cladogenetic_events_observed$event_type == "subset (s)")
	v_observed = sum(table_of_cladogenetic_events_observed$event_type == "vicariance (v)")
	j_observed = sum(table_of_cladogenetic_events_observed$event_type == "founder (j)")

	d_actual = sum(table_of_range_change_events_translated$event_type == "expansion (d)")
	e_actual = sum(table_of_range_change_events_translated$event_type == "contraction (e)")
	a_actual = sum(table_of_range_change_events_translated$event_type == "range-switching (a)")

	d_observed = sum(table_of_range_change_events_observed$event_type == "expansion (d)")
	e_observed = sum(table_of_range_change_events_observed$event_type == "contraction (e)")
	a_observed = sum(table_of_range_change_events_observed$event_type == "range-switching (a)")
	
	simstats = c(num_fossil_tips, num_extinct_nodes, y_actual, s_actual, v_actual, j_actual, y_observed, s_observed, v_observed, j_observed, d_actual, e_actual, a_actual, d_observed, e_observed, a_observed)

	simstats = adf2(matrix(data=simstats, nrow=1))
	names(simstats) = c("num_fossil_tips", "num_extinct_nodes", "y_actual", "s_actual", "v_actual", "j_actual", "y_observed", "s_observed", "v_observed", "j_observed", "d_actual", "e_actual", "a_actual", "d_observed", "e_observed", "a_observed")
	
	
	
	#######################################################
	# Output object
	#######################################################
	SSEsim_results_processed = NULL
	
	# Also save the original results, as there was some info there
	SSEsim_results_processed$SSEsim_results_raw = SSEsim_results
	
	SSEsim_results_processed$simtr = simtr
	SSEsim_results_processed$simtr_observed = simtr_observed

	SSEsim_results_processed$ancstate_probs_simfull = ancstate_probs_simfull
	SSEsim_results_processed$ancstate_probs_simobs = ancstate_probs_simobs

	SSEsim_results_processed$tipranges_simulated = tipranges_simulated
	SSEsim_results_processed$tipranges_observed = tipranges_observed

	if (trait_TF == TRUE)
		{
		SSEsim_results_processed$traits_simulated = traits_simulated
		SSEsim_results_processed$traits_observed = traits_observed
		} # END if (trait_TF == TRUE)


	SSEsim_results_processed$table_of_cladogenetic_events_translated = table_of_cladogenetic_events_translated
	SSEsim_results_processed$table_of_cladogenetic_events_observed = table_of_cladogenetic_events_observed

	SSEsim_results_processed$table_of_range_change_events_translated = table_of_range_change_events_translated
	SSEsim_results_processed$table_of_range_change_events_observed = table_of_range_change_events_observed

	
	# Saving nodes that are extinct (translate later)
	SSEsim_results_processed$table_of_cladogenetic_events_extinct = table_of_cladogenetic_events_extinct
	SSEsim_results_processed$table_of_cladogenetic_events_w_1_descendant = table_of_cladogenetic_events_w_1_descendant
	# Translating full simtree nodenums to observed simtree nodenums
	SSEsim_results_processed$simtr_full_to_observed_node_translation = simtr_full_to_observed_node_translation

	SSEsim_results_processed$simstats = simstats
	SSEsim_results_processed$simdir = simdir

	if (write_files == TRUE)
		{
		write.tree(phy=simtr, file="simtr_complete.newick")
		write.tree(phy=simtr_observed, file="simtr_observed.newick")
		write.tree(phy=simtr_observed, file="tree.newick")
		
		save_tipranges_to_LagrangePHYLIP(tipranges=tipranges_simulated, lgdata_fn="geog_sim_complete.txt")
		save_tipranges_to_LagrangePHYLIP(tipranges=tipranges_observed, lgdata_fn="geog_sim_observed.txt")
		save_tipranges_to_LagrangePHYLIP(tipranges=tipranges_observed, lgdata_fn="geog.data")
		
		if (trait_TF == TRUE)
			{
			areanames = LETTERS[1:num_trait_states]
			save_tipranges_to_LagrangePHYLIP(tipranges=traits_simulated, lgdata_fn="traits_sim_complete.txt", areanames=areanames)
			save_tipranges_to_LagrangePHYLIP(tipranges=traits_observed, lgdata_fn="traits_sim_observed.txt", areanames=areanames)
			save_tipranges_to_LagrangePHYLIP(tipranges=traits_observed, lgdata_fn="traits.data", areanames=areanames)

			# Also make geog files with just 1 area, always occupied
			# (needed input for a traits-only analysis)
			tipranges_simulated_1area = tipranges_simulated
			tipranges_observed_1area = tipranges_observed
		
			# Take just the first column
			tipranges_simulated_1area@df = dfnums_to_numeric(tipranges_simulated_1area@df)
			tipranges_observed_1area@df = dfnums_to_numeric(tipranges_observed_1area@df)

			tmpmat = tipranges_simulated_1area@df
			tmpmat[,1][tmpmat[,1] != 1] = 1
			tmpmat = as.data.frame(matrix(tmpmat[,1], ncol=1, byrow=TRUE), stringsAsFactors=FALSE)
			tipranges_simulated_1area@df = tmpmat
			names(tipranges_simulated_1area@df) = names(tipranges_simulated)[1]
			row.names(tipranges_simulated_1area@df) = row.names(tipranges_simulated@df)
	
			tmpmat = tipranges_observed_1area@df
			tmpmat[,1][tmpmat[,1] != 1] = 1
			tmpmat = as.data.frame(matrix(tmpmat[,1], ncol=1, byrow=TRUE), stringsAsFactors=FALSE)
			tipranges_observed_1area@df = tmpmat
			names(tipranges_observed_1area@df) = names(tipranges_observed)[1]
			row.names(tipranges_observed_1area@df) = row.names(tipranges_observed@df)

			save_tipranges_to_LagrangePHYLIP(tipranges=tipranges_simulated_1area, lgdata_fn="tipranges_simulated_1area.data", areanames=c("A"))
			save_tipranges_to_LagrangePHYLIP(tipranges=tipranges_observed_1area, lgdata_fn="tipranges_observed_1area.data", areanames=c("A"))
			} # END if (trait_TF == TRUE)
		} # END if (write_files == TRUE)
	
	# Return to original working directory
	setwd(orig_wd)

	return(SSEsim_results_processed)
	} # end SSEsim_to_files()



# Label the recorded anagenetic events as d, e...
label_table_of_anagenetic_events <- function(table_of_range_change_events_translated, states_list, ranges_list, areanames, track_dispersal_dest=TRUE, m=NULL)
	{
	# Check for all NAs
	if (all(is.na(table_of_range_change_events_translated)))
		{
		stop("\n\nError in label_table_of_anagenetic_events(): input 'table_of_range_change_events_translated' was all NAs!\n\n")
		}


	# Get the m values (dispersal multipliers, dependent on trait)
	# Trait states
	if (is.null(m) == FALSE)
		{
		num_trait_states = length(m)
	
		# Numbering of geog vs. trait states
		num_geog_states = length(states_list)
	
		trait_state_indices_list_1based = NULL
		max_index_1based = 0
		for (i in 1:length(m))
			{
			tmp_indices_list = (max_index_1based+1) : (max_index_1based + num_geog_states)
			max_index_1based = max(tmp_indices_list)
			trait_state_indices_list_1based[[i]] = tmp_indices_list
			}
	
		numstates = length(states_list) * length(m)
		} else {
		numstates = length(states_list)
		} # END if (is.null(m) == FALSE)

	
	starting_ranges = as.numeric(table_of_range_change_events_translated$starting_range)
	new_ranges = as.numeric(table_of_range_change_events_translated$new_range)



	numevents = nrow(table_of_range_change_events_translated)
	event_type = rep(NA, numevents)
	event_txt = rep(NA, numevents)
	
	if (is.null(m) == FALSE)
		{
		trait_events = rep(NA, numevents)
		}
	
	
	# List the destination of dispersal events
	if (track_dispersal_dest == TRUE)
		{
		dispersal_to = rep(NA, numevents)
		extirpation_from = rep(NA, numevents)
		} # END if (track_dispersal_dest == TRUE)
	
	for (i in 1:numevents)
		{
		# Text describing the event
		#print(i)
		#print(starting_ranges)
		#print(starting_ranges[i])
		#print(ranges_list[[starting_ranges[i]]])

		# Translate the starting ranges:
		if (is.null(m) == FALSE)
			{
			# Translate the starting range
			TFs = rep(FALSE, times=num_trait_states)
			for (j in 1:num_trait_states)
				{
				TF = starting_ranges[i] %in% trait_state_indices_list_1based[[j]]
				TFs[j] = TF
				}
			starting_trait_statenum = (1:num_trait_states)[TFs]
			starting_geog_state_tmp = (starting_ranges[i] - ((starting_trait_statenum-1)*num_geog_states) )
			starting_ranges[i] = starting_geog_state_tmp
			
			# Translate the ending range
			TFs = rep(FALSE, times=num_trait_states)
			for (j in 1:num_trait_states)
				{
				TF = new_ranges[i] %in% trait_state_indices_list_1based[[j]]
				TFs[j] = TF
				}
			ending_trait_statenum = (1:num_trait_states)[TFs]
			ending_geog_state_tmp = (new_ranges[i] - ((ending_trait_statenum-1)*num_geog_states) )
			new_ranges[i] = ending_geog_state_tmp
			
			trait_events[i] = paste0(starting_trait_statenum, "->", ending_trait_statenum)
			} # END if (is.null(m) == FALSE)


		
		event_txt[i] = paste(ranges_list[[starting_ranges[i]]], "->", ranges_list[[new_ranges[i]]], sep="")
		
		# "dispersal" (range expansion)
		if (length(states_list[[starting_ranges[i]]]) < length(states_list[[new_ranges[i]]]))
			{
			event_type[i] = "expansion (d)"
			if (track_dispersal_dest == TRUE)
				{
				states_list_in_new_range_0based = states_list[[new_ranges[i]]]
				states_list_in_old_range_0based = states_list[[starting_ranges[i]]]
				
				new_area_TF = (states_list_in_new_range_0based %in% states_list_in_old_range_0based) == FALSE
				states_list_in_new_range_0based_just_new_area = states_list[[new_ranges[i]]][new_area_TF]
				dispersal_to[i] = areanames[1+states_list_in_new_range_0based_just_new_area]
				extirpation_from[i] = NA
				} # END if (track_dispersal_dest == TRUE)
			next()
			}
		
		# "extinction" (range contraction / local extirpation)
		if (length(states_list[[starting_ranges[i]]]) > length(states_list[[new_ranges[i]]]))
			{
			event_type[i] = "contraction (e)"

			# Find and record the area that was lost
			if (track_dispersal_dest == TRUE)
				{
				states_list_in_new_range_0based = states_list[[new_ranges[i]]]
				states_list_in_old_range_0based = states_list[[starting_ranges[i]]]
				lostarea_TF = (states_list_in_old_range_0based %in% states_list_in_new_range_0based) == FALSE
				states_list_in_new_range_0based_just_area_lost = states_list[[starting_ranges[i]]][lostarea_TF]
				dispersal_to[i] = NA
				extirpation_from[i] = areanames[1+states_list_in_new_range_0based_just_area_lost]
				} # END if (track_dispersal_dest == TRUE)
			next()
			}

		# If ranges are the same length, this would be anagenesis
		# BUT, if there are traits, it could be a trait-switch instead
		if (length(states_list[[starting_ranges[i]]]) == length(states_list[[new_ranges[i]]]))
			{
			# Trait switch possible
			if (is.null(m) == FALSE)
				{
				if (starting_trait_statenum != ending_trait_statenum)
					{
					event_type[i] = "trait change (t)"
					next()
					}
				} # END if (is.null(m) == FALSE)
			
			
			# Contraction will result in the same length, check for this
			if (ranges_list[[new_ranges[i]]] == "_")
				{
				event_type[i] = "contraction (e)"
				
				# Find and record the area that was lost
				if (track_dispersal_dest == TRUE)
					{
					states_list_in_new_range_0based = states_list[[new_ranges[i]]]
					states_list_in_old_range_0based = states_list[[starting_ranges[i]]]
					lostarea_TF = (states_list_in_old_range_0based %in% states_list_in_new_range_0based) == FALSE
					states_list_in_new_range_0based_just_area_lost = states_list[[starting_ranges[i]]][lostarea_TF]
					dispersal_to[i] = NA
					extirpation_from[i] = areanames[1+states_list_in_new_range_0based_just_area_lost]
					} # END if (track_dispersal_dest == TRUE)

				next()				
				} else {
				event_type[i] = "range-switching (a)"

				# Find and record the area that was "lost", and "gained"
				if (track_dispersal_dest == TRUE)
					{
					states_list_in_new_range_0based = states_list[[new_ranges[i]]]
					states_list_in_old_range_0based = states_list[[starting_ranges[i]]]
					
					# Area that was gained
					new_area_TF = (states_list_in_new_range_0based %in% states_list_in_old_range_0based) == FALSE
					states_list_in_new_range_0based_just_new_area = states_list[[new_ranges[i]]][new_area_TF]
					dispersal_to[i] = areanames[1+states_list_in_new_range_0based_just_new_area]
					
					# Area that was lost
					lostarea_TF = (states_list_in_old_range_0based %in% states_list_in_new_range_0based) == FALSE
					states_list_in_new_range_0based_just_area_lost = states_list[[starting_ranges[i]]][lostarea_TF]
					extirpation_from[i] = areanames[1+states_list_in_new_range_0based_just_area_lost]
					} # END if (track_dispersal_dest == TRUE)

				next()
				}
			}
		} # end for-loop
	
	if (track_dispersal_dest == TRUE)
		{
		anagenetic_event_labels = as.data.frame(cbind(event_type, event_txt, dispersal_to, extirpation_from), stringsAsFactors=FALSE)
		anagenetic_event_labels
		} else {
		anagenetic_event_labels = as.data.frame(cbind(event_type, event_txt), stringsAsFactors=FALSE)
		anagenetic_event_labels
		} # END if (track_dispersal_dest == TRUE)

	# Add the trait-change events
	if (is.null(m) == FALSE)
		{
		anagenetic_event_labels = cbind(anagenetic_event_labels, trait_events)		
		} # END if (is.null(m) == FALSE)
	
	return(anagenetic_event_labels)
	} # end label_table_of_anagenetic_events()







# Label the recorded cladogenetic events as j, v, s, y...
label_table_of_cladogenetic_events <- function(table_of_cladogenetic_events_translated, states_list, ranges_list, areanames=areanames, track_dispersal_dest=TRUE, m=NULL)
	{

	# Get the m values (dispersal multipliers, dependent on trait)
	# Trait states
	if (is.null(m) == FALSE)
		{
		num_trait_states = length(m)
	
		# Numbering of geog vs. trait states
		num_geog_states = length(states_list)
	
		trait_state_indices_list_1based = NULL
		max_index_1based = 0
		for (i in 1:length(m))
			{
			tmp_indices_list = (max_index_1based+1) : (max_index_1based + num_geog_states)
			max_index_1based = max(tmp_indices_list)
			trait_state_indices_list_1based[[i]] = tmp_indices_list
			}
	
		numstates = length(states_list) * length(m)
		} else {
		numstates = length(states_list)
		} # END if (is.null(m) == FALSE)





	Parent_geog_states = table_of_cladogenetic_events_translated$parent_range
	Left_geog_states = table_of_cladogenetic_events_translated$Left_state
	Right_geog_states = table_of_cladogenetic_events_translated$Right_state

	numevents = nrow(table_of_cladogenetic_events_translated)
	event_type = rep(NA, numevents)
	event_txt = rep(NA, numevents)

	if (is.null(m) == FALSE)
		{
		trait_events = rep(NA, numevents)
		}



	# List the destination of dispersal events
	if (track_dispersal_dest == TRUE)
		{
		dispersal_to = rep("", numevents)
		# No extirpation at cladogenesis
		} # END if (track_dispersal_dest == TRUE)
	
	for (i in 1:numevents)
		{
		# Text describing the event
		
		# No trait
		if (is.null(m) == TRUE)
			{
			Parent_geog_state_tmp = Parent_geog_states[i]
			Left_geog_state_tmp = Left_geog_states[i]
			Right_geog_state_tmp = Right_geog_states[i]
			}
		
		# Yes trait -- this requires "classifying" the state number down
		# to the no-traits state number (i.e., the ranges)
		if (is.null(m) == FALSE)
			{
			# Parent state
			TFs = rep(FALSE, times=num_trait_states)
			for (j in 1:num_trait_states)
				{
				TF = Parent_geog_states[i] %in% trait_state_indices_list_1based[[j]]
				TFs[j] = TF
				}
			parent_trait_statenum = (1:num_trait_states)[TFs]
			Parent_geog_state_tmp = (Parent_geog_states[i] - ((parent_trait_statenum-1)*num_geog_states) )
# 			print(Parent_geog_states[i])
# 			print(parent_trait_statenum)
# 			print(parent_trait_statenum-1)
# 			print(num_geog_states)
# 			print(Parent_geog_state_tmp)
			
			
			# Left state
			TFs = rep(FALSE, times=num_trait_states)
			for (j in 1:num_trait_states)
				{
				TF = Left_geog_states[i] %in% trait_state_indices_list_1based[[j]]
				TFs[j] = TF
				}
			Left_trait_statenum = (1:num_trait_states)[TFs]
			Left_trait_statenum = (1:num_trait_states)[TFs]
			Left_geog_state_tmp = (Left_geog_states[i] - ((Left_trait_statenum-1)*num_geog_states) )

			# Right state
			TFs = rep(FALSE, times=num_trait_states)
			for (j in 1:num_trait_states)
				{
				TF = Right_geog_states[i] %in% trait_state_indices_list_1based[[j]]
				TFs[j] = TF
				}
			Right_trait_statenum = (1:num_trait_states)[TFs]
			Right_trait_statenum = (1:num_trait_states)[TFs]
			Right_geog_state_tmp = (Right_geog_states[i] - ((Right_trait_statenum-1)*num_geog_states) )
			
			trait_events_tmp = paste0(parent_trait_statenum, "->", Left_trait_statenum, ",", Right_trait_statenum)
			trait_events[i] = trait_events_tmp
			} # END if (is.null(m) == FALSE)
		
		
		
		event_txt[i] = paste(ranges_list[[Parent_geog_state_tmp]], "->", ranges_list[[Left_geog_state_tmp]], ",",  ranges_list[[Right_geog_state_tmp]], sep="")
		
		#print("Left_geog_states[i]:")
		#print(Left_geog_states[i])
		#print("Right_geog_states[i]:")
		#print(Right_geog_states[i])
		
		#print("table_of_cladogenetic_events_translated:")
		#print(table_of_cladogenetic_events_translated)
		
		# Sympatry
		if (Left_geog_state_tmp == Right_geog_state_tmp)
			{
			event_type[i] = "sympatry (y)"
			next()
			}
		
		parent_indices = sort(states_list[[Parent_geog_state_tmp]])
		Left_indices = states_list[[Left_geog_state_tmp]]
		Right_indices = states_list[[Right_geog_state_tmp]]
		desc_merged = sort(c(Left_indices, Right_indices))
		desc_indices = sort(unique(c(Left_indices, Right_indices)))
		
		# Founder-events (j, jump dispersal)
		if (length(desc_indices) > length(parent_indices))
			{
			event_type[i] = "founder (j)"

			if (track_dispersal_dest == TRUE)
				{
				TF = desc_indices %in% parent_indices
				new_range_TF = TF == FALSE
				new_area_indices_1based = 1+desc_indices[new_range_TF]
				areas_txt = areanames[unlist(new_area_indices_1based)]
				range_txt = paste(areas_txt, collapse="", sep="")
				dispersal_to[i] = range_txt
				#dispersal_to[i] = ranges_list[[new_indices]]
				
				} # END if (track_dispersal_dest == TRUE)

			next()			
			}
		
		# Vicariance
		if ((length(parent_indices) == length(desc_merged)) && (all(parent_indices == desc_merged) == TRUE))
			{
			event_type[i] = "vicariance (v)"
			next()						
			} else {
			event_type[i] = "subset (s)"
			next()			
			}
		
		} # end for-loop
		
		


	if (track_dispersal_dest == TRUE)
		{
		cladogenetic_event_labels = as.data.frame(cbind(event_type, event_txt, dispersal_to), stringsAsFactors=FALSE)
		cladogenetic_event_labels
		} else {
		cladogenetic_event_labels = as.data.frame(cbind(event_type, event_txt), stringsAsFactors=FALSE)
		cladogenetic_event_labels
		} # END if (track_dispersal_dest == TRUE)

	if (is.null(m) == FALSE)
		{
		cladogenetic_event_labels = cbind(cladogenetic_event_labels, trait_events)
		}

	
	
	return(cladogenetic_event_labels)
	} # end label_table_of_cladogenetic_events()












edges_table_to_tipranges_object <- function(edges_table2, states_list, ntips, tiplabels_for_tipranges, areanames, addval=0)
	{
	defaults = '
	tiplabels = paste("sp", 1:ntips, sep="")
	areanames = c("K", "O", "M", "H")
	addval=0
	'
	binary_table = matrix(data=0, nrow=ntips, ncol=length(areanames))

	for (i in 1:ntips)
		{
		# Get the index of the range, from the simulation (APE nodenums)
		tmp_range_index = edges_table2$edge_ranges[i]
		
		# Get the 1-based indexes of the areas
		area_indices_1based = states_list[[tmp_range_index]] + 1 + addval
		
		# Convert to binary array
		binary_table[i, area_indices_1based] = 1
		}

	binary_table

	# Convert to tipranges object
	tmpdf2 = adf2(data.matrix(binary_table))
	names(tmpdf2) = areanames
	rownames(tmpdf2) = tiplabels_for_tipranges
	tmpdf2
	
	tipranges_simulated = define_tipranges_object(tmpdf=tmpdf2)
	tipranges_simulated
	
	return(tipranges_simulated)
	}








#######################################################
# Prepare for plotting SSE simulations
#######################################################

# Put the results of an SSEsim into a BioGeoBEARS results_object (res)
# This assumes certain default filenames in the directory
SSEsim_results_processed_into_res <- function(res, SSEsim_results_processed, simtr_complete_or_observed="complete")
	{
	defaults='
	# Prelim
	SSEsim_results_fn = "SSEsim_results_processed.Rdata"
	# Loads to SSEsim_results_processed
	load(SSEsim_results_fn)
	SSEsim_results_processed
	
	# Get a res object
	DEC_inf_fn = "DEC_inf.Rdata"
	load(DEC_inf_fn)	# Loads to DEC_inf
	DEC_inf
	res = DEC_inf
	
	# This function
	simtr_complete_or_observed = "complete"
	
	' # end defaults
	
	# Make a modified results object
	resmod = res
	numstates = ncol(resmod$ML_marginal_prob_each_state_at_branch_top_AT_node)
	
	# Get basic tree info
	if (simtr_complete_or_observed == "complete")
		{
		tmptr = SSEsim_results_processed$simtr
		} else {
		tmptr = SSEsim_results_processed$simtr_observed
		}

	# Tree info
	ntips = length(tmptr$tip.label)
	num_internal_nodes = tmptr$Nnode
	tipnums = 1:ntips
	internal_nodenums = (ntips+1):(ntips+num_internal_nodes)
	all_nodenums = c(tipnums, internal_nodenums)
	numnodes_all = length(all_nodenums)

	
	# Null out a lot of the items
	resmod$computed_likelihoods_at_each_node = NULL
	resmod$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS = NULL
	resmod$condlikes_of_each_state = NULL 
	resmod$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS = NULL
	resmod$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS = NULL
	resmod$relative_probs_of_each_state_at_branch_top_AT_node_UPPASS = NULL
	resmod$relative_probs_of_each_state_at_bottom_of_root_branch = NULL

	# Get the simulated ancestral state "probabilities" (1 or 0, really, since it's simulated)
	if (simtr_complete_or_observed == "complete")
		{
		# States at all tips & internal nodes (as probabilities)
		ancstate_probs_sim = SSEsim_results_processed$ancstate_probs_simfull
	
		# Setup getting the ancstates at the bottom of branches
		# These are the split events at the top of branches (internal nodes only)
		ancstates = SSEsim_results_processed$table_of_cladogenetic_events_translated
	
		} else {
		# States at all tips & internal nodes (as probabilities)
		ancstate_probs_sim = SSEsim_results_processed$ancstate_probs_simobs

		# Setup getting the ancstates at the bottom of branches
		# These are the split events at the top of branches (internal nodes only)
		ancstates = SSEsim_results_processed$table_of_cladogenetic_events_observed
		
		internal_nodenums = (length(tmptr$tip.label)+tmptr$Nnode)
		
		matchnums = match(x=ancstates$ape_decnode, table=SSEsim_results_processed$simtr_full_to_observed_node_translation$full)
		
		ancstates$ape_decnode = SSEsim_results_processed$simtr_full_to_observed_node_translation$obs[matchnums]
		
		
		for (a in 1:length(ancstates$ape_decnode))
			{
			if (is.na(ancstates$ape_decnode[a]) == TRUE)
				{
				ancstates$ape_ancnode[a] = NA
				next()
				}
			
			nodenum_to_get_ancestors_of = ancstates$ape_decnode[a]
			edge_ending_in_nodenum_TF = tmptr$edge[,2] == nodenum_to_get_ancestors_of
			anc_nodenum = tmptr$edge[edge_ending_in_nodenum_TF, 1]
			
			if (length(anc_nodenum) > 0)
				{
				#print(ancstates$ape_ancnode[a])
				#print(anc_nodenum)
				ancstates$ape_ancnode[a] = anc_nodenum
				} else {
				ancstates$ape_ancnode[a] = NA
				}
			}
		ancstates
		}
	dim(ancstate_probs_sim)
	dim(ancstates)

	# Fill these in with the observed states (0)
	resmod$ML_marginal_prob_each_state_at_branch_bottom_below_node = matrix(data=0, nrow=numnodes_all, ncol=numstates)
	resmod$ML_marginal_prob_each_state_at_branch_top_AT_node = matrix(data=0, nrow=numnodes_all, ncol=numstates)

	for (nodenum in 1:numnodes_all)
		{
		stateprobs_for_this_node = ancstate_probs_sim
		resmod$ML_marginal_prob_each_state_at_branch_top_AT_node[nodenum,] = stateprobs_for_this_node[nodenum,]
	
		# Getting the states at branch bottoms below nodes is more complex
		# Is the current nodenum an internal node stored in ancstates?
		#print(ancstates$ape_decnode)
		#print(nodenum)
		internal_TF = ancstates$ape_decnode == nodenum
		internal_TF[is.na(internal_TF)] = FALSE	# 2014-06-04_NJM
		
		if (sum(internal_TF) > 1)
			{
			stop("Error: sum(internal_TF) > 1")
			} # if (sum(internal_TF) > 1)
	
		# Find the left and right decnodenums
		if (sum(internal_TF) == 1)
			{
			row_TF = ancstates$ape_decnode == nodenum
			row_TF[is.na(row_TF)] = FALSE	# 2014-06-04_NJM
			rownum = (1:length(ancstates$ape_decnode))[row_TF]
			tmp_nodenum_to_find_decs_of = ancstates$ape_decnode[rownum]
			
			# Get the edges that have nodenum as the ancestor
			edgerow_TF_decnodes = tmptr$edge[,1] == nodenum
			if (sum(edgerow_TF_decnodes, na.rm=TRUE) == 0)
				{
				next()
				}
			
			# Descendant nodenums
			dec_nodenums = tmptr$edge[edgerow_TF_decnodes,2]
			dec_nodenums
		
			# Hope to God these are always left, right...
			Left_dec_nodenum = dec_nodenums[1]
			Right_dec_nodenum = dec_nodenums[2]
		
			# Store the probabilities of each state at the bottom of Left descending branch
			state_1based_at_branch_bottom_of_Left_decnode = ancstates$Left_state[rownum]
			probs_at_branch_bottom_of_Left_decnode = rep(0, times=numstates)
			probs_at_branch_bottom_of_Left_decnode[state_1based_at_branch_bottom_of_Left_decnode] = 1

			# Store the probabilities of each state at the bottom of Right descending branch
			state_1based_at_branch_bottom_of_Right_decnode = ancstates$Right_state[rownum]
			probs_at_branch_bottom_of_Right_decnode = rep(0, times=numstates)
			probs_at_branch_bottom_of_Right_decnode[state_1based_at_branch_bottom_of_Right_decnode] = 1
		
			# Store in the results matrix for branch bottoms
			resmod$ML_marginal_prob_each_state_at_branch_bottom_below_node[Left_dec_nodenum,] = probs_at_branch_bottom_of_Left_decnode
			resmod$ML_marginal_prob_each_state_at_branch_bottom_below_node[Right_dec_nodenum,] = probs_at_branch_bottom_of_Right_decnode
			} # if (sum(internal_TF) == 1)
		} # for (nodenum in 1:numnodes_all)

	rowSums(resmod$ML_marginal_prob_each_state_at_branch_top_AT_node)
	rowSums(resmod$ML_marginal_prob_each_state_at_branch_bottom_below_node)


	# Change the input tree and geography files
	if (simtr_complete_or_observed == "complete")
		{
		resmod$inputs$trfn = "simtr_complete.newick"
		resmod$inputs$geogfn = "geog_sim_observed.txt"
		} else {
		resmod$inputs$trfn = "simtr_observed.newick"
		resmod$inputs$geogfn = "geog_sim_observed.txt"
		} # END if (simtr_complete_or_observed == "complete")
		
	return(resmod)
	} # END SSEsim_results_processed_into_res



get_dej_params_row_from_SSEsim_results_processed <- function(SSEsim_results_processed)
	{

	d = SSEsim_results_processed$SSEsim_results_raw$SSEsim_inputs$BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"]
	e = SSEsim_results_processed$SSEsim_results_raw$SSEsim_inputs$BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"]
	j = SSEsim_results_processed$SSEsim_results_raw$SSEsim_inputs$BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"]
	brate = SSEsim_results_processed$SSEsim_results_raw$SSEsim_inputs$SSEmodel$brate
	drate = SSEsim_results_processed$SSEsim_results_raw$SSEsim_inputs$SSEmodel$drate
	rangesize_b_exponent = SSEsim_results_processed$SSEsim_results_raw$SSEsim_inputs$SSEmodel$rangesize_b_exponent
	rangesize_d_exponent = SSEsim_results_processed$SSEsim_results_raw$SSEsim_inputs$SSEmodel$rangesize_d_exponent

	dej_params_row = matrix(data=c(d, e, j, brate, drate, rangesize_b_exponent, rangesize_d_exponent), nrow=1)
	dej_params_row = adf2(dej_params_row)
	names(dej_params_row) = c("d", "e", "j", "brate", "drate", "rangesize_b_exponent", "rangesize_d_exponent")
	dej_params_row
	
	return(dej_params_row)
	} # END get_dej_params_row_from_SSEsim_results_processed



SSEsim_results_processed_to_anagenetic_events_table <- function(SSEsim_results_processed, simtr_complete_or_observed="complete", tmptr_table=NULL, convert_clado_events_w1desc=FALSE)
	{
	defaults='
	simtr_complete_or_observed = "complete"
	tmptr_table = NULL
	convert_clado_events_w1desc=FALSE
	'
	# OK, now our events table is set.
	# Translate the anagenetic events table
	# I.e., convert:
	# 
	# 1. SSEsim_results_processed$table_of_range_change_events_translated
	# 
	# ...to...
	# 
	# 2. stochastic_mapping_results$table_w_anagenetic_events

	# Put in the tip states, from:
	if (simtr_complete_or_observed == "complete")
		{
		SSEsim_anagenetic_table = SSEsim_results_processed$table_of_range_change_events_translated
		tmptr = SSEsim_results_processed$simtr
		
		# Make a prt table
		if (is.null(tmptr_table))
			{
			tmptr_table = prt(tmptr, printflag=FALSE)
			} # END if (is.null(tmptr_table))
		
		} else {
		# This will be more complex, since you have to put in missing speciation events
		# and how they changed the range anagenetically
		
		if (convert_clado_events_w1desc == FALSE)
			{
			SSEsim_anagenetic_table = SSEsim_results_processed$table_of_range_change_events_observed
			} else {
			SSEsim_anagenetic_table = SSEsim_results_processed$clado_1desc_to_anagenetic_table
			}
		tmptr = SSEsim_results_processed$simtr_observed

		# Make a prt table
		if (is.null(tmptr_table))
			{
			tmptr_table = prt(tmptr, printflag=FALSE)
			} # END if (is.null(tmptr_table))

		} # if (simtr_complete_or_observed == "complete")
	
	##################################################
	# Fix factor classes - 2018-01-04
	tmp_cls_df = cls.df(SSEsim_anagenetic_table)
	factor_TF = tmp_cls_df$cls_col_list == "factor"
	if (sum(factor_TF) > 0)
		{
		colnums = (1:nrow(tmp_cls_df))[factor_TF]
		for (i in 1:length(colnums))
			{
			SSEsim_anagenetic_table[,colnums[i]] = as.character(SSEsim_anagenetic_table[,colnums[i]])
			} # END for (i in 1:length(colnums))
		} # END if (sum(factor_TF) > 0)
	cls.df(SSEsim_anagenetic_table)
	##################################################
		
	# Setup
	areanames = SSEsim_results_processed$SSEsim_results_raw$areanames
	ranges_list = SSEsim_results_processed$SSEsim_results_raw$ranges_list
	states_list = SSEsim_results_processed$SSEsim_results_raw$states_list
	states_list_0based = states_list
	
	# Check for NA in the anagenetic events table (no anagenetic events)
	if (all(is.na(SSEsim_anagenetic_table)) == TRUE)
		{
		return(NA)
		}
	
	SSEsim_anagenetic_table
	node = SSEsim_anagenetic_table$ape_decnode
	current_rangenum_1based = SSEsim_anagenetic_table$starting_range
	new_rangenum_1based = SSEsim_anagenetic_table$new_range
	current_rangetxt = ranges_list[current_rangenum_1based]
	new_rangetxt = ranges_list[new_rangenum_1based]

	event_type = SSEsim_anagenetic_table$event_type
	event_txt = SSEsim_anagenetic_table$event_txt
	dispersal_to = SSEsim_anagenetic_table$dispersal_to
	extirpation_from = SSEsim_anagenetic_table$extirpation_from
	dispersal_to[is.na(dispersal_to)] = "-"
	extirpation_from[is.na(extirpation_from)] = "-"

	# Label event types as just "d", "e", or "a"
	event_type[event_type == "expansion (d)"] = "d"
	event_type[event_type == "contraction (e)"] = "e"
	event_type[event_type == "range-switching (a)"] = "a"

	# The event times in stochastic simulation are calculated up
	# from the root.  Thus, get root age and subtract
	root_age = get_max_height_tree(tmptr)
	abs_event_time = root_age - (SSEsim_anagenetic_table$t + SSEsim_anagenetic_table$dt)

	# event time is the time above the bottom of the branch
	event_time = tmptr_table$time_bp[node] + tmptr_table$edge.length[node] - abs_event_time
	event_time


	# Get the 1-based numbers for the areas added/lost
	new_area_num_1based = rep(NA, times=length(dispersal_to))
	lost_area_num_1based = rep(NA, times=length(dispersal_to))
	areanums = (1:length(areanames))
	for (i in 1:length(dispersal_to))
		{
		dispersal_area_TF = (areanames == dispersal_to[i])
		if (sum(dispersal_area_TF) == 1)
			{
			new_area_num_1based[i] = areanums[dispersal_area_TF]
			}
		lost_area_TF = (areanames == extirpation_from[i])
		if (sum(lost_area_TF) == 1)
			{
			lost_area_num_1based[i] = areanums[lost_area_TF]
			}
		} # END for (i in 1:length(dispersal_to))


	# Make it like in stochastic mapping
	nodenum_at_top_of_branch = node
	trynum = rep(1, times=length(node))
	brlen = tmptr_table$edge.length[node]

	table_w_anagenetic_events = cbind(nodenum_at_top_of_branch, trynum, brlen, current_rangenum_1based, new_rangenum_1based, current_rangetxt, new_rangetxt, abs_event_time, event_time, event_type, event_txt, new_area_num_1based, lost_area_num_1based, dispersal_to, extirpation_from)
	table_w_anagenetic_events = dfnums_to_numeric(adf2(table_w_anagenetic_events))
	table_w_anagenetic_events
	
	return(table_w_anagenetic_events)
	} # END SSEsim_results_processed_to_anagenetic_events_table





# Convert 
# - SSEsim_results_translated$table_of_cladogenetic_events_translated
# ...to...
# - stochastic_mapping_results$table_cladogenetic_events
SSEsim_results_processed_to_cladogenetic_events_table <- function(SSEsim_results_processed, simtr_complete_or_observed="complete", tmptr_table=NULL)
	{

	# Tree info
	if (simtr_complete_or_observed == "complete")
		{
		tmptr = SSEsim_results_processed$simtr

		# Make a prt table
		if (is.null(tmptr_table))
			{
			tmptr_table = prt(tmptr, printflag=FALSE)
			} # END if (is.null(tmptr_table))

		# Get the table of SSE-simulated cladogenetic events
		SSEsim_clado_table = SSEsim_results_processed$table_of_cladogenetic_events_translated

		} else {
		# This will be more complex, since you have to put in missing speciation events
		# and how they changed the range anagenetically
		tmptr = SSEsim_results_processed$simtr_observed

		# Make a prt table
		if (is.null(tmptr_table))
			{
			tmptr_table = prt(tmptr, printflag=FALSE)
			} # END if (is.null(tmptr_table))

		# Get the table of SSE-simulated cladogenetic events,
		# just those on the observed tree
		SSEsim_clado_table = SSEsim_results_processed$table_of_cladogenetic_events_observed

		} # END if (simtr_complete_or_observed == "complete")

	ntips = length(tmptr$tip.label)
	num_internal_nodes = tmptr$Nnode
	tipnums = 1:ntips
	internal_nodenums = (ntips+1):(ntips+num_internal_nodes)
	all_nodenums = c(tipnums, internal_nodenums)
	numnodes_all = length(all_nodenums)


	table_cladogenetic_events = NULL

	for (nodenum in 1:numnodes_all)
		{
		node = nodenum
	
		node_in_SSE_internal_nodes_TF = SSEsim_clado_table$ape_decnode == nodenum
	
		if (sum(node_in_SSE_internal_nodes_TF) == 1)
			{
			rownum = (1:nrow(SSEsim_clado_table))[node_in_SSE_internal_nodes_TF]
			sampled_states_AT_nodes = SSEsim_clado_table$parent_range[rownum]

		
			# Get nodes at tops of branches descending from node
			decnode_rows_TF = tmptr$edge[,1] == nodenum
			dec_nodenums = tmptr$edge[decnode_rows_TF, 2]
			left_desc_nodes = dec_nodenums[1]	# hope this works
			right_desc_nodes = dec_nodenums[2]	# hope this works
			#print(dec_nodenums)
		
			# Get the parent node of this one
			parent_nodenum = SSEsim_clado_table$ape_ancnode[rownum]
			# This parent had two descendants; so, fill this in later
			sampled_states_AT_brbots = NA
		
			samp_LEFT_dcorner = SSEsim_clado_table$Left_state[rownum]
			samp_RIGHT_dcorner = SSEsim_clado_table$Right_state[rownum]
			
			clado_event_type = SSEsim_clado_table$event_type[rownum]
			clado_event_txt = SSEsim_clado_table$event_txt[rownum]
			clado_dispersal_to = SSEsim_clado_table$dispersal_to[rownum]
			
			} else {
			sampled_states_AT_nodes = NA
			sampled_states_AT_brbots = NA
			left_desc_nodes = NA
			right_desc_nodes = NA
			samp_LEFT_dcorner = NA
			samp_RIGHT_dcorner = NA
			clado_event_type = NA
			clado_event_txt = NA
			clado_dispersal_to = NA
			} # END if (sum(node_in_SSE_internal_nodes_TF) == 1)
		tmprow = c(node, sampled_states_AT_nodes, sampled_states_AT_brbots, left_desc_nodes, right_desc_nodes, samp_LEFT_dcorner, samp_RIGHT_dcorner, clado_event_type, clado_event_txt, clado_dispersal_to)
		table_cladogenetic_events = rbind(table_cladogenetic_events, tmprow)
		} # END for (nodenum in 1:numnodes_all)

	# Now, fill in sampled_states_AT_brbots
	table_cladogenetic_events = adf2(table_cladogenetic_events)
	names(table_cladogenetic_events) = c("node", "sampled_states_AT_nodes", "sampled_states_AT_brbots", "left_desc_nodes", "right_desc_nodes", "samp_LEFT_dcorner", "samp_RIGHT_dcorner", "clado_event_type", "clado_event_txt", "clado_dispersal_to")

	#table_cladogenetic_events[300:320,]

	for (nodenum in 1:numnodes_all)
		{
		# Pass these corner states to the branch_bottoms of the nodes above
	
		# Left
		node_to_change_brbottom_of = as.numeric(table_cladogenetic_events$left_desc_nodes[nodenum])
		if (!is.na(node_to_change_brbottom_of))
			{
			table_cladogenetic_events$sampled_states_AT_brbots[node_to_change_brbottom_of] = table_cladogenetic_events$samp_LEFT_dcorner[nodenum]
			} # END if (!is.na(node_to_change_brbottom_of))
	
		# Right
		node_to_change_brbottom_of = as.numeric(table_cladogenetic_events$right_desc_nodes[nodenum])
		if (!is.na(node_to_change_brbottom_of))
			{
			table_cladogenetic_events$sampled_states_AT_brbots[node_to_change_brbottom_of] = table_cladogenetic_events$samp_RIGHT_dcorner[nodenum]
			} # END if (!is.na(node_to_change_brbottom_of))
		} # for (nodenum in 1:numnodes_all)

	head(table_cladogenetic_events)
	table_cladogenetic_events[300:320,]
	
	
	# Put in the tip states, from:
	if (simtr_complete_or_observed == "complete")
		{
		ancstate_probs = SSEsim_results_processed$ancstate_probs_simfull
		} else {
		ancstate_probs = SSEsim_results_processed$ancstate_probs_simobs
		}

	numstates = ncol(ancstate_probs)
	for (nodenum in 1:ntips)
		{
		statenums = 1:numstates
		tipstate_1based = statenums[ancstate_probs[nodenum,] == 1]
		table_cladogenetic_events$sampled_states_AT_nodes[nodenum] = tipstate_1based
		} # END for (nodenum in 1:ntips)

	table_cladogenetic_events
	return(table_cladogenetic_events)
	} # SSEsim_results_processed_to_cladogenetic_events_table


convert_SSEsim_to_stochastic_mapping_results_format <- function(SSEsim_results_processed, simtr_complete_or_observed="complete", tmptr_table=NULL, include_clado_events_w1desc=TRUE)
	{
	# Tree info
	if (simtr_complete_or_observed == "complete")
		{
		tmptr = SSEsim_results_processed$simtr

		# Make a prt table
		if (is.null(tmptr_table))
			{
			tmptr_table = prt(tmptr, printflag=FALSE)
			} # END if (is.null(tmptr_table))

		} else {
		# This will be more complex, since you have to put in missing speciation events
		# and how they changed the range anagenetically
		tmptr = SSEsim_results_processed$simtr_observed

		# Make a prt table
		if (is.null(tmptr_table))
			{
			tmptr_table = prt(tmptr, printflag=FALSE)
			} # END if (is.null(tmptr_table))

		} # END if (simtr_complete_or_observed == "complete")
	
	
	# Construct the table of cladogenetic events
	table_cladogenetic_events = SSEsim_results_processed_to_cladogenetic_events_table(SSEsim_results_processed=SSEsim_results_processed, simtr_complete_or_observed=simtr_complete_or_observed, tmptr_table=tmptr_table)
tail(table_cladogenetic_events)

	# Save the anagenetic events in a table
	# Get a stochastic-mapping-formatted table of anagenetic events
	table_w_anagenetic_events = SSEsim_results_processed_to_anagenetic_events_table(SSEsim_results_processed, simtr_complete_or_observed=simtr_complete_or_observed, tmptr_table=tmptr_table, convert_clado_events_w1desc=FALSE)
	table_w_anagenetic_events
	
	
	# Get the table of cladogenetic events that left 1 descendant
	if (include_clado_events_w1desc == TRUE)
		{
		table_w_clado_1desc_events = SSEsim_results_processed_to_anagenetic_events_table(SSEsim_results_processed, simtr_complete_or_observed=simtr_complete_or_observed, tmptr_table=tmptr_table, convert_clado_events_w1desc=TRUE)
		
		
		head(table_w_clado_1desc_events)
		head(table_w_anagenetic_events)
		
		# Merge the two anagenetic tables
		if ( (length(table_w_anagenetic_events) == 1) && (is.na(table_w_anagenetic_events) ) == FALSE)
			{
			table_w_anagenetic_events = rbind(table_w_anagenetic_events, table_w_clado_1desc_events)
			} else {
			table_w_anagenetic_events = table_w_clado_1desc_events
			}
		} # END if (include_clado_events_w1desc == TRUE)
	
	# Add the text version of the anagenetic events table to the cladogenesis events table
	# Convert the anagenetic events table into anagenetic_events_txt
	anagenetic_events_txt_below_node = rep(NA, times=nrow(table_cladogenetic_events))
	if ( ((length(table_w_anagenetic_events) == 1) && (is.na(table_w_anagenetic_events))) == FALSE)
		{
		for (i in 1:nrow(table_cladogenetic_events))
			{
			tmp_nodenum = table_cladogenetic_events$node[i]
			anagenetic_events_table_rowsTF = table_w_anagenetic_events$nodenum_at_top_of_branch == tmp_nodenum
			# Correct for NAs, turn them to FALSE
			#anagenetic_events_table_rowsTF[is.na(anagenetic_events_table_rowsTF)] = FALSE
		
			events_table_for_branch = table_w_anagenetic_events[anagenetic_events_table_rowsTF,]
			anagenetic_events_txt_below_node[i] = events_table_into_txt(events_table_for_branch=events_table_for_branch)
			}
		}
	anagenetic_events_txt_below_node[anagenetic_events_txt_below_node == ""] = "none"
	anagenetic_events_txt_below_node[is.na(anagenetic_events_txt_below_node)] = "none"
	anagenetic_events_txt_below_node

	# Store the txt version of anagenetic events
	table_cladogenetic_events = cbind(tmptr_table, table_cladogenetic_events, anagenetic_events_txt_below_node)
	table_cladogenetic_events = adf2(table_cladogenetic_events)
	table_cladogenetic_events = dfnums_to_numeric(table_cladogenetic_events)
	table_cladogenetic_events
	
	# Store and return the output
	SSEsim_converted_to_stochastic_mapping_results_format = NULL
	SSEsim_converted_to_stochastic_mapping_results_format$master_table_cladogenetic_events = table_cladogenetic_events
	SSEsim_converted_to_stochastic_mapping_results_format$table_w_anagenetic_events = table_w_anagenetic_events
	
	return(SSEsim_converted_to_stochastic_mapping_results_format)
	} # END convert_SSEsim_to_stochastic_mapping_results_format





#######################################################
# Convert the node numbers between full and observed trees
# (and edge numbers)
#######################################################

convert_full_simtree_nodenums_to_observed_simtree_nodenums <- function(tmp_events_table, simtr_full_to_observed_node_translation, simtr_complete, simtr_observed, simtr_complete_table=NULL, simtr_observed_table=NULL)
	{
	if (all(is.na(tmp_events_table)) == TRUE)
		{
		# Pass along the NA
		tmp_events_table2 = NA
		return(tmp_events_table2)
		}
	
	# Setup -- you WILL need the node tables
	if (is.null(simtr_complete_table))
		{
		simtr_complete_table = prt(simtr_complete, printflag=FALSE, get_tipnames=TRUE, fossils_older_than=fossils_older_than)
		}
	if (is.null(simtr_observed_table))
		{
		simtr_observed_table = prt(simtr_observed, printflag=FALSE, get_tipnames=TRUE, fossils_older_than=fossils_older_than)
		}

	
	#######################################################
	# Anagenetic events on fulltr nodes that are found in the observed tree
	#######################################################
	full_simtr_ape_ancnode = tmp_events_table$ape_ancnode
	full_simtr_ape_decnode = tmp_events_table$ape_decnode
	full_simtr_sim_edgenum_w_event = tmp_events_table$sim_edgenum_w_event
	
	rownums_in_conversion_table = match(x=full_simtr_ape_decnode, table=simtr_full_to_observed_node_translation$full)
	rownums_in_conversion_table
	
	# Nodenums in the observed tree
	ape_decnode = simtr_full_to_observed_node_translation$obs[rownums_in_conversion_table]
	
	# Get the ancestor nodes in the observed tree
	edgerow_nums = match(x=ape_decnode, table=simtr_observed$edge[,2])
	
	# NA is presumably the root node, which has no ancestors
	edgerow_nums
	
	# Ancestor nodenums in observed tree
	ape_ancnode = simtr_observed$edge[edgerow_nums,1]
	
	
	# Convert the edgenums
	sim_edgenum_w_event = edgerow_nums
	
	# Save the results
	tmp_events_table2 = tmp_events_table
	tmp_events_table2$ape_ancnode = ape_ancnode
	tmp_events_table2$ape_decnode = ape_decnode
	tmp_events_table2$sim_edgenum_w_event = sim_edgenum_w_event
	
	# Add the full_simtr columns
	tmp_events_table2 = cbind(full_simtr_ape_ancnode, full_simtr_ape_decnode, full_simtr_sim_edgenum_w_event, tmp_events_table2)


	#######################################################
	# Anagenetic events on fulltr nodes that are just clado_1desc in the observed tree
	# 2014-06-03_NJM fix
	#######################################################
	for (i in 1:nrow(tmp_events_table2))
		{
		# Skip nodes where the decnode in the observed tree has already been found
		if (!is.na(tmp_events_table2$ape_decnode[i]))
			{
			next()
			}
		
		# Get the nodenum in the full sim tree
		full_simtr_ape_decnode = tmp_events_table2$full_simtr_ape_ancnode[i]

		# Find the node above that is in the observed tree
		# 1. Get the descending tips in the full tree
		tipnames_in_full_tree_str = simtr_complete_table$tipnames[full_simtr_ape_decnode]
		tipnames_in_full_tree = strsplit(tipnames_in_full_tree_str, split=",")[[1]]
		
		# 2. Reduce this list to the list of what's left in the observed tree
		tips_in_observed_tr_TF = tipnames_in_full_tree %in% simtr_observed$tip.label
		tipnames_in_full_tree_reduced_to_obs = tipnames_in_full_tree[tips_in_observed_tr_TF]
		
		# 3. Match this list to some node in the observed tree table
		tipnames_in_full_tree_reduced_to_obs_string = paste(sort(tipnames_in_full_tree_reduced_to_obs), collapse=",", sep="")
		tipnames_in_full_tree_reduced_to_obs_string
		
		
		# Match this string to a string in the observed tree
		matching_rownum = match(x=tipnames_in_full_tree_reduced_to_obs_string, table=simtr_observed_table$tipnames)

		if (is.na(matching_rownum))
			{
			errortxt = paste("ERROR in convert_full_simtree_nodenums_to_observed_simtree_nodenums(): You shouldn't get NA from searching node ", full_simtr_ape_decnode, "'s full tree tipnames against the observed tree.\n\n", sep="")
			stop(errortxt)
			}

		if (length(matching_rownum) != 1)
			{
			errortxt = paste("ERROR in convert_full_simtree_nodenums_to_observed_simtree_nodenums(): You should only get 1 hit from searching node ", full_simtr_ape_decnode, "'s full tree tipnames against the observed tree.\n\n", sep="")
			stop(errortxt)
			}
		
		observed_tree_row = simtr_observed_table[matching_rownum,]
		ape_decnode = observed_tree_row$node
		ape_ancnode = observed_tree_row$ancestor
		sim_edgenum_w_event = observed_tree_row$parent_br
		
		tmp_events_table2$ape_decnode[i] = ape_decnode
		tmp_events_table2$ape_ancnode[i] = ape_ancnode
		tmp_events_table2$sim_edgenum_w_event[i] = sim_edgenum_w_event
		
		} # for (i in 1:nrow(tmp_events_table2))
	# Now the anagenetic below nodes that are only clado_1desc should have
	# appropriate observed tree nodenums
	
	return(tmp_events_table2)
	} # END convert_full_simtree_nodenums_to_observed_simtree_nodenums





#######################################################
# Convert the table of cladogenetic events with only 
# 1 descendant into an anagenetic events table
#######################################################
convert_clado_1desc_to_anagenetic <- function(clado_events_w1desc_table, simtr_full_to_observed_node_translation, simtr_complete, simtr_observed, simtr_complete_table=NULL, simtr_observed_table=NULL, fossils_older_than=0.001)
	{
	
	if (is.null(simtr_complete_table))
		{
		simtr_complete_table = prt(simtr_complete, printflag=FALSE, get_tipnames=TRUE, fossils_older_than=fossils_older_than)
		}
	if (is.null(simtr_observed_table))
		{
		simtr_observed_table = prt(simtr_observed, printflag=FALSE, get_tipnames=TRUE, fossils_older_than=fossils_older_than)
		}
	
	
	
	
	# Get the list of tipnames for each node in the observed tree
	tipnames_list_for_each_node_in_observed_tree_strings = simtr_observed_table$tipnames
	tipnames_list_for_each_node_in_observed_tree_list = NULL
	numspecies_for_each_node_in_observed_tree = NULL
	for (j in 1:length(tipnames_list_for_each_node_in_observed_tree_strings))
		{
		tmp_tipnames = strsplit(tipnames_list_for_each_node_in_observed_tree_strings, split=",")[[1]]
		tipnames_list_for_each_node_in_observed_tree_list[[j]] = tmp_tipnames
		
		numspecies_for_each_node_in_observed_tree = c(numspecies_for_each_node_in_observed_tree, length(tmp_tipnames))
		}
	
	numrows = nrow(clado_events_w1desc_table)
	if (numrows == 0)
		{
		clado_events_w1desc_table2 = NA
		return(clado_events_w1desc_table2)
		} else {
		clado_events_w1desc_table2 = NULL
		} # END if (numrows == 0)
	
	for (i in 1:numrows)
		{
		full_simtr_ape_ancnode = clado_events_w1desc_table$ape_ancnode[i]
		full_simtr_ape_decnode = clado_events_w1desc_table$ape_decnode[i]
		complete_simtr_row = simtr_complete_table[full_simtr_ape_decnode, ]
		
		# Find the node above that is in the observed tree
		# 1. Get the descending tips in the full tree
		tipnames_in_full_tree_str = simtr_complete_table$tipnames[full_simtr_ape_decnode]
		tipnames_in_full_tree = strsplit(tipnames_in_full_tree_str, split=",")[[1]]
		
		# 2. Reduce this list to the list of what's left in the observed tree
		tips_in_observed_tr_TF = tipnames_in_full_tree %in% simtr_observed$tip.label
		tipnames_in_full_tree_reduced_to_obs = tipnames_in_full_tree[tips_in_observed_tr_TF]
		
		# 3. Match this list to some node in the observed tree table
		tipnames_in_full_tree_reduced_to_obs_string = paste(sort(tipnames_in_full_tree_reduced_to_obs), collapse=",", sep="")
		tipnames_in_full_tree_reduced_to_obs_string
		
		
		# Match this string to a string in the observed tree
		matching_rownum = match(x=tipnames_in_full_tree_reduced_to_obs_string, table=simtr_observed_table$tipnames)
		
		if (is.na(matching_rownum))
			{
			errortxt = paste("ERROR in convert_clado_1desc_to_anagenetic(): You shouldn't get NA from searching node ", full_simtr_ape_decnode, "'s full tree tipnames against the observed tree.\n\n", sep="")
			stop(errortxt)
			}

		if (length(matching_rownum) != 1)
			{
			errortxt = paste("ERROR in convert_clado_1desc_to_anagenetic(): You should only get 1 hit from searching node ", full_simtr_ape_decnode, "'s full tree tipnames against the observed tree.\n\n", sep="")
			stop(errortxt)
			}
		
		# Otherwise you have a match in the observed tree
		observed_tree_row = simtr_observed_table[matching_rownum,]
		ape_decnode = observed_tree_row$node
		ape_ancnode = observed_tree_row$ancestor
		sim_edgenum_w_event = observed_tree_row$parent_br
		
		if (is.na(complete_simtr_row$edge.length) == TRUE)
			{
			complete_simtr_branch_bottom_abs_time = complete_simtr_row$time_bp + 0
			} else {
			complete_simtr_branch_bottom_abs_time = complete_simtr_row$time_bp + complete_simtr_row$edge.length
			}
			
		
		
		
		# Put the age of the simtr_full nodes as the event times
		# The times_bp remain the same
		clado_event_w1desc_abs_time_bp = complete_simtr_row$time_bp
		time_above_observed_root = get_max_height_tree(simtr_observed) - clado_event_w1desc_abs_time_bp
		
		# 't' is actually THE ABSOLUTE HEIGHT ABOVE THE ROOT
		t = time_above_observed_root
		
		# How far above the bottom of this observed branch, was
		# the cladogenesis event with 1 descendant?
		# dt is really "how long until the next event", which we don't know
		# for these cladogenesis events
		#dt = obs_simtr_t - clado_event_w1desc_abs_time_bp
		dt = 0
		
		# Which of the daughter nodenums in the complete simtree
		# survived to be sampled
		complete_simtr_daughter_nodenums = complete_simtr_row$daughter_nds[[1]]
		
		#############################################################
		# Check the daughters of this node in the complete tree to
		# see which survived in the observed tree
		#############################################################
		# Check the left node
		# Find the node above that is in the observed tree
		# 1. Get the descending tips in the full tree
		tipnames_in_full_tree_str = simtr_complete_table$tipnames[complete_simtr_daughter_nodenums[1]]
		tipnames_in_full_tree = strsplit(tipnames_in_full_tree_str, split=",")[[1]]
		
		# 2. Reduce this list to the list of what's left in the observed tree
		tips_in_observed_tr_TF = tipnames_in_full_tree %in% simtr_observed$tip.label
		tipnames_in_full_tree_reduced_to_obs = tipnames_in_full_tree[tips_in_observed_tr_TF]
		
		# 3. Match this list to some node in the observed tree table
		tipnames_in_full_tree_reduced_to_obs_string = paste(sort(tipnames_in_full_tree_reduced_to_obs), collapse=",", sep="")
		tipnames_in_full_tree_reduced_to_obs_string
		
		matching_rownum1 = match(x=tipnames_in_full_tree_reduced_to_obs_string, table=simtr_observed_table$tipnames)
		#############################################################


		#############################################################
		# Check the right node
		# Find the node above that is in the observed tree
		# 1. Get the descending tips in the full tree
		tipnames_in_full_tree_str = simtr_complete_table$tipnames[complete_simtr_daughter_nodenums[2]]
		tipnames_in_full_tree = strsplit(tipnames_in_full_tree_str, split=",")[[1]]
		
		# 2. Reduce this list to the list of what's left in the observed tree
		tips_in_observed_tr_TF = tipnames_in_full_tree %in% simtr_observed$tip.label
		tipnames_in_full_tree_reduced_to_obs = tipnames_in_full_tree[tips_in_observed_tr_TF]
		
		# 3. Match this list to some node in the observed tree table
		tipnames_in_full_tree_reduced_to_obs_string = paste(sort(tipnames_in_full_tree_reduced_to_obs), collapse=",", sep="")
		tipnames_in_full_tree_reduced_to_obs_string
		
		matching_rownum2 = match(x=tipnames_in_full_tree_reduced_to_obs_string, table=simtr_observed_table$tipnames)
		#############################################################
		
		
		daughter_survived_nodenums = c(matching_rownum1, matching_rownum2)
		daughter_survived_TF = !is.na(daughter_survived_nodenums)
		
		if (sum(daughter_survived_TF) != 1)
			{
			errortxt = paste("\n\nERROR in convert_clado_1desc_to_anagenetic: complete simtree node #", full_simtr_ape_decnode, " should have 1 surviving daughter\nin the observed tree, if the input was really a 'clado_events_w1desc_table'.\n\nHowever, the number of suriviving daughters at this node was ", sum(daughter_survived_TF), ".\n\n", sep="")
			stop(errortxt)
			}
		
		descendant_state_after_cladogenesis_in_complete_simtr = c(clado_events_w1desc_table$Left_state[i], clado_events_w1desc_table$Right_state[i])
		descendant_state_after_cladogenesis_observed = descendant_state_after_cladogenesis_in_complete_simtr[daughter_survived_TF]
		
		
		starting_range = clado_events_w1desc_table$parent_range[i]
		new_range = descendant_state_after_cladogenesis_observed
		
		event_type = clado_events_w1desc_table$event_type[i]
		event_txt = clado_events_w1desc_table$event_txt[i]
		dispersal_to = clado_events_w1desc_table$dispersal_to[i]
		extirpation_from = NA


		tmprow_clado_event_converted_to_anagenetic = c(full_simtr_ape_ancnode, full_simtr_ape_decnode, ape_ancnode, ape_decnode, sim_edgenum_w_event, t, dt, starting_range, new_range, event_type, event_txt, dispersal_to, extirpation_from)		
		clado_events_w1desc_table2 = rbind(clado_events_w1desc_table2, tmprow_clado_event_converted_to_anagenetic)
		
		clado_events_w1desc_table2
		
		} # END for (i in 1:numrows)
	
	
	clado_events_w1desc_table2 = adf2(clado_events_w1desc_table2)
	names(clado_events_w1desc_table2) = c("full_simtr_ape_ancnode", "full_simtr_ape_decnode ", "ape_ancnode", "ape_decnode", "sim_edgenum_w_event", "t", "dt", "starting_range", "new_range", "event_type", "event_txt", "dispersal_to", "extirpation_from")
	clado_events_w1desc_table2 = dfnums_to_numeric(clado_events_w1desc_table2)
	
	return(clado_events_w1desc_table2)
	} # END convert_clado_1desc_to_anagenetic



#######################################################
# Check/fix d events
#######################################################
# Note: If, during SSE simulation, the first event on a
# a branch which is later merged in the observed tree
# is a "d" or "e", the SSE simulation code seems
# to save its time correctly (as '$t', absolute time above
# root), however, it seems to place its branch position
# as the ancestral branch below.  
#
# This may be happening during translation from the 
# full simtr to the reduced, observed tree when the tree
# events tables are converted to stochastic mapping format
# for painting with paint_stochastic maps.
#
# The resulting plots have an obvious error where d events 
# plot above the branch tops of the branches ancestral
# to where they should plot.
# 
# Until I figure out the original labeling issue, it is
# simplest to just use this function to re-assign the 
# ape_decnode of these nodes.  This involves deciding
# which of the two branches they belong on, which is 
# helped by the fact that the states should
# match, and a d event cannot be followed by a y
# event.
#######################################################

correct_de_events_on_observed_tree <- function(table_of_range_change_events_observed, clado_events_w1desc_table2, simtr_observed_table, printflag=FALSE)
	{

	if (printflag)
		{
		cat("\n\nCorrecting d/e events placed 'below' ancestral node rather than descendant node:\n")
		} # END if (printflag)

	# Cut out d events with t (age above observed tree root) less than 0
	table_of_range_change_events_observed = table_of_range_change_events_observed[table_of_range_change_events_observed$t > 0,]

	# Check for d events with t above the time_bp of the node they are on
	anagenetic_event_obs_tree_nodenums = table_of_range_change_events_observed$ape_decnode

	d_event_too_high_TF = simtr_observed_table$node_ht[anagenetic_event_obs_tree_nodenums] <= table_of_range_change_events_observed$t
	cbind(simtr_observed_table$node_ht[anagenetic_event_obs_tree_nodenums], table_of_range_change_events_observed$t, d_event_too_high_TF)

	# For these, the ape_decnode is wrong and is referring to the ancestor
	# Get the correct branch decnode and use that
	# (This traces back to some difficult issue with 
	#  labeling in the original simulation)
	ds_to_move = table_of_range_change_events_observed[d_event_too_high_TF,]
	
	
	if (nrow(ds_to_move) > 0)
		{
		# Hopefully clado_events_w1desc_table2 will never be NA when nrow(ds_to_move) > 0 
		for (m in 1:nrow(ds_to_move))
			{
			d_to_move = ds_to_move[m,]
			dec_nodenum_to_change = d_to_move$ape_decnode
		
			# Get the daughter nodenums
			daughter_dec_nodenums = simtr_observed_table$daughter_nds[dec_nodenum_to_change][[1]]
			
			
			# cladogenetic events on left daughter:
			if (!is.na(clado_events_w1desc_table2))
				{
				TF = clado_events_w1desc_table2$ape_decnode == daughter_dec_nodenums[1]
				tmprownums = (1:nrow(clado_events_w1desc_table2))[TF]
				tmp_events_table = clado_events_w1desc_table2[tmprownums,]
		
				if (nrow(tmp_events_table) != 0)
					{
					# Cladogenetic events on the daughter branch, sorted
					tmp_events_table = tmp_events_table[order(tmp_events_table$t), ]
		
					# Is the first event after the problematic event?
					age_left_1st_event = tmp_events_table$t[1]
					left_TF1 = age_left_1st_event > d_to_move$t
		
					# Is the first event a possible cladogenetic range reduction?
					left_TF2 = tmp_events_table$event_type[1] != "sympatry (y)"
		
					# Does the new range of the ancestral event match the 
					# current range of the desc. event?
					left_TF3 = d_to_move$new_range == tmp_events_table$starting_range[1]
		
					left_TF_count = (left_TF1 + left_TF2 + left_TF3)
				
					if (printflag)
						{
						cat("\nLeft:\n")
						cat(c(left_TF1, left_TF2, left_TF3))
						cat("\n")
						} # END if (printflag)

					} else {
					left_TF_count = 0
					age_left_1st_event = 0.2
					}
		
	
	
		
				# cladogenetic events on right daughter:
				TF = clado_events_w1desc_table2$ape_decnode == daughter_dec_nodenums[2]
				tmprownums = (1:nrow(clado_events_w1desc_table2))[TF]
				tmp_events_table = clado_events_w1desc_table2[tmprownums,]

				if (nrow(tmp_events_table) != 0)
					{
					# Cladogenetic events on the daughter branch, sorted
					tmp_events_table = tmp_events_table[order(tmp_events_table$t), ]
		
					# Is the first event after the problematic event?
					age_right_1st_event = tmp_events_table$t[1]
					right_TF1 = age_right_1st_event > d_to_move$t
		
					# Is the first event a possible cladogenetic range reduction?
					right_TF2 = tmp_events_table$event_type[1] != "sympatry (y)"
		
					# Does the new range of the ancestral event match the 
					# current range of the desc. event?
					right_TF3 = d_to_move$new_range == tmp_events_table$starting_range[1]

					right_TF_count = (right_TF1 + right_TF2 + right_TF3)

					if (printflag)
						{
						cat("\nRight:\n")
						cat(c(right_TF1, right_TF2, right_TF3))
						cat("\n")
						} # END if (printflag)
		
					} else {
					right_TF_count = 0
					age_left_1st_event = 0.1
					} # END if (nrow(tmp_events_table) != 0)
				} else {
				left_TF_count = 0
				age_left_1st_event = 0.2
				right_TF_count = 0
				age_left_1st_event = 0.1
				} # END if (!is.na(clado_events_w1desc_table2))
		
		
			# Pick one of the daughters
			if (left_TF_count > right_TF_count)
				{
				if (printflag)
					{
					cat("Correcting a d event label, moving to LEFT desc. branch.\n\n")
					} # END printflag
				table_of_range_change_events_observed$ape_decnode[d_event_too_high_TF][m] = daughter_dec_nodenums[1]
				
				# Get the ape_ancnode and edge number
				table_of_range_change_events_observed$ape_ancnode[d_event_too_high_TF][m] = simtr_observed_table$ancestor[daughter_dec_nodenums[1]]
				table_of_range_change_events_observed$sim_edgenum_w_event[d_event_too_high_TF][m] = simtr_observed_table$parent_br[daughter_dec_nodenums[1]]
				} # END if (left_TF_count > right_TF_count)

			if (left_TF_count < right_TF_count)
				{
				if (printflag)
					{
					cat("Correcting a d event label, moving to RIGHT desc. branch.\n\n")
					} # END printflag
				table_of_range_change_events_observed$ape_decnode[d_event_too_high_TF][m] = daughter_dec_nodenums[2]

				# Get the ape_ancnode and edge number
				table_of_range_change_events_observed$ape_ancnode[d_event_too_high_TF][m] = simtr_observed_table$ancestor[daughter_dec_nodenums[2]]
				table_of_range_change_events_observed$sim_edgenum_w_event[d_event_too_high_TF][m] = simtr_observed_table$parent_br[daughter_dec_nodenums[2]]
				} # END if (left_TF_count < right_TF_count)

		
			# Warning -- algorithm can't decide for sure
			if (left_TF_count == right_TF_count)
				{
				if (age_left_1st_event >= age_right_1st_event)
					{
					# Pick the older one
					if (printflag)
						{
						cat("Correcting a d event label, tie score, picking left branch event as this is later.\n\n")
						} # END if (printflag)
					table_of_range_change_events_observed$ape_decnode[d_event_too_high_TF][m] = daughter_dec_nodenums[1]

					# Get the ape_ancnode and edge number
					table_of_range_change_events_observed$ape_ancnode[d_event_too_high_TF][m] = simtr_observed_table$ancestor[daughter_dec_nodenums[1]]
					table_of_range_change_events_observed$sim_edgenum_w_event[d_event_too_high_TF][m] = simtr_observed_table$parent_br[daughter_dec_nodenums[1]]

					} else {
					if (printflag)
						{
						cat("Correcting a d event label, tie score, picking right branch event as this is later.\n\n")
						} # END if (printflag)
					table_of_range_change_events_observed$ape_decnode[d_event_too_high_TF][m] = daughter_dec_nodenums[2]

					# Get the ape_ancnode and edge number
					table_of_range_change_events_observed$ape_ancnode[d_event_too_high_TF][m] = simtr_observed_table$ancestor[daughter_dec_nodenums[2]]
					table_of_range_change_events_observed$sim_edgenum_w_event[d_event_too_high_TF][m] = simtr_observed_table$parent_br[daughter_dec_nodenums[2]]
					} # END if (age_left_1st_event >= age_right_1st_event)
				} # END if (left_TF_count == right_TF_count)

			# Warning
# 			if ((left_TF == FALSE) && (right_TF == FALSE))
# 				{
# 				if (age_left_1st_event >= age_right_1st_event)
# 					{
# 					# Pick the older one
# 					cat("Correcting a d event label, picking2 left branch.\n\n")
# 					table_of_range_change_events_observed$ape_decnode[d_event_too_high_TF][m] = daughter_dec_nodenums[1]
# 					} else {
# 					cat("Correcting a d event label, picking2 right branch.\n\n")
# 					table_of_range_change_events_observed$ape_decnode[d_event_too_high_TF][m] = daughter_dec_nodenums[2]
# 					}
# 				} # END if ((left_TF == FALSE) && (right_TF == FALSE))
			} # END for (m in 1:nrow(ds_to_move))
		} # END if (nrow(ds_to_move) > 0)

	# Double-check
	anagenetic_event_obs_tree_nodenums = table_of_range_change_events_observed$ape_decnode
	d_event_too_high_TF = simtr_observed_table$node_ht[anagenetic_event_obs_tree_nodenums] <= table_of_range_change_events_observed$t
	cbind(simtr_observed_table$node_ht[anagenetic_event_obs_tree_nodenums], table_of_range_change_events_observed$t, d_event_too_high_TF)
		
	return(table_of_range_change_events_observed)
	} # END correct_de_events_on_observed_tree




#######################################################
# Fix node numbering for observed matrices in the observed tree
#######################################################

convert_simtr_observed_events_to_prelim_stochastic_mapping_format <- function(SSEsim_results_processed, fossils_older_than = 0.001, simtr_complete_table=NULL, simtr_observed_table=NULL, printflag=FALSE)
	{
	# Fix node numbering for observed matrices in the observed tree

	defaults='
	# SSEsim_results_processed = orig_SSEsim_results_processed
	
	fossils_older_than = 0.001
	' # END defauls
	
	if (printflag)
		{
		cat("\n\nConverting SSEsim events into a format readable into convert_SSEsim_to_stochastic_mapping_results_format()...\n\n")
		cat("Processing (inputting simtr_complete_table/simtr_observed_table will speed this up)...\n\n")
		}
	
	# Preliminary
	SSEsim_results_processed
	names(SSEsim_results_processed)

	# Observed cladogenetic events
	dim(SSEsim_results_processed$table_of_cladogenetic_events_observed)

	# Observed anagenetic events
	dim(SSEsim_results_processed$table_of_range_change_events_observed)

	# Cladogenetic events converted to anagenetic in the observed tree
	#dim(SSEsim_results_processed$clado_events_converted_to_anagenetic)

	# Translation between node numbers
	dim(SSEsim_results_processed$simtr_full_to_observed_node_translation)

	# All nodes not sampled in the observed tree
	dim(SSEsim_results_processed$table_of_cladogenetic_events_extinct)

	# Node node sampled in the observed tree, but which had a descendant
	dim(SSEsim_results_processed$table_of_cladogenetic_events_w_1_descendant)


	# Get the tree tables
	#fossils_older_than = 0.001
	fossils_older_than = fossils_older_than
	simtr_observed = SSEsim_results_processed$simtr_observed
	simtr_complete = SSEsim_results_processed$simtr
	
	if (is.null(simtr_observed_table))
		{
		simtr_observed_table = prt(simtr_observed, printflag=FALSE, get_tipnames=TRUE, fossils_older_than=fossils_older_than)
		}
	if (is.null(simtr_complete_table))
		{
		simtr_complete_table = prt(simtr_complete, printflag=FALSE, get_tipnames=TRUE, fossils_older_than=fossils_older_than)
		}
	
	# Convert node/edge numbers of observed cladogenetic events
	tmp_events_table = SSEsim_results_processed$table_of_cladogenetic_events_observed
	
	simtr_full_to_observed_node_translation = SSEsim_results_processed$simtr_full_to_observed_node_translation
	
	# Error fix, if making the simtr_full_to_observed_node_translation
	# was skipped for some reason
	if (is.null(SSEsim_results_processed$simtr_full_to_observed_node_translation))
		{
		simtr_full_to_observed_node_translation = get_simtr_full_to_observed_node_translation(simtr_complete, simtr_observed, simtr_observed_table)
		SSEsim_results_processed$simtr_full_to_observed_node_translation = simtr_full_to_observed_node_translation
		} # if (is.null(SSEsim_results_processed$simtr_full_to_observed_node_translation))
		
	table_of_cladogenetic_events_observed = convert_full_simtree_nodenums_to_observed_simtree_nodenums(tmp_events_table, simtr_full_to_observed_node_translation, simtr_complete=simtr_complete, simtr_observed=simtr_observed, simtr_complete_table=simtr_complete_table, simtr_observed_table=simtr_observed_table)

	
	# Update SSEsim results:
	#SSEsim_results_processed$table_of_cladogenetic_events_observed_OLD = SSEsim_results_processed$table_of_cladogenetic_events_observed
	SSEsim_results_processed$table_of_cladogenetic_events_observed = table_of_cladogenetic_events_observed


	#######################################################
	# Convert node/edge numbers of observed clado_1desc events
	#######################################################

	# Convert the cladogenetic events with 1 descendant into anagenetic events
	table_of_cladogenetic_events_w_1_descendant = SSEsim_results_processed$table_of_cladogenetic_events_w_1_descendant
	clado_events_w1desc_table = table_of_cladogenetic_events_w_1_descendant
	simtr_full_to_observed_node_translation = SSEsim_results_processed$simtr_full_to_observed_node_translation
	
	
	if (length(clado_events_w1desc_table) == 0)
		{
		clado_events_w1desc_table2 = NA
		} else {
		if ( ( (length(clado_events_w1desc_table) == 1) && (is.na(clado_events_w1desc_table)) ) == FALSE)
			{
			clado_events_w1desc_table2 = convert_clado_1desc_to_anagenetic(clado_events_w1desc_table, simtr_full_to_observed_node_translation, simtr_complete, simtr_observed, simtr_complete_table, simtr_observed_table)

			# Cut the cladogenetic-1-descendant events below the observed tree root node
			# (which will be a lot of nodes!)
			clado_events_w1desc_table2 = clado_events_w1desc_table2[!is.na(clado_events_w1desc_table2$sim_edgenum_w_event),]
			head(clado_events_w1desc_table2)

			# Remove commas from cladogenetic events_txt
			clado_events_w1desc_table2$event_txt = gsub(pattern=",", replacement="_", x=clado_events_w1desc_table2$event_txt)

			clado_events_w1desc_table2$dispersal_to[clado_events_w1desc_table2$dispersal_to == ""] = "-"
			clado_events_w1desc_table2$dispersal_to[is.na(clado_events_w1desc_table2$dispersal_to)] = "-"
			clado_events_w1desc_table2$extirpation_from[clado_events_w1desc_table2$extirpation_from == ""] = "-"
			clado_events_w1desc_table2$extirpation_from[is.na(clado_events_w1desc_table2$extirpation_from)] = "-"
			} else {
			clado_events_w1desc_table2 = NA
			} # END if (!is.na(clado_events_w1desc_table))
		} # END if (length(clado_events_w1desc_table) == 0)
		
	# Save in SSEsim results
	SSEsim_results_processed$clado_1desc_to_anagenetic_table = clado_events_w1desc_table2



	#######################################################
	# Convert node/edge numbers of observed anagenetic events
	#######################################################
	tmp_events_table = SSEsim_results_processed$table_of_range_change_events_observed
	
	if (all(is.na(tmp_events_table)) == FALSE)
		{
		simtr_full_to_observed_node_translation = SSEsim_results_processed$simtr_full_to_observed_node_translation
		table_of_range_change_events_observed = convert_full_simtree_nodenums_to_observed_simtree_nodenums(tmp_events_table, simtr_full_to_observed_node_translation, simtr_complete=simtr_complete, simtr_observed=simtr_observed, simtr_complete_table=simtr_complete_table, simtr_observed_table=simtr_observed_table)

		#######################################################
		# Check/fix d events
		#######################################################
		table_of_range_change_events_observed2 = correct_de_events_on_observed_tree(table_of_range_change_events_observed=table_of_range_change_events_observed, clado_events_w1desc_table2=clado_events_w1desc_table2, simtr_observed_table=simtr_observed_table, printflag=printflag)
		head(table_of_range_change_events_observed2)
		
		} else {
		table_of_range_change_events_observed = NA
		table_of_range_change_events_observed2 = NA
		}

	SSEsim_results_processed$table_of_range_change_events_observed = table_of_range_change_events_observed2
	
	return(SSEsim_results_processed)
	} # END convert_simtr_observed_events_to_prelim_stochastic_mapping_format
	



