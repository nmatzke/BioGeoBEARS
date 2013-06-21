test='
require("ape")
require("rexpokit")
require("cladoRcpp")
require(BioGeoBEARS)
source("/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_classes_v1.R", chdir = TRUE)
source("/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_univ_model_v1.R", chdir = TRUE)
extdata_dir = system.file("extdata/", package="BioGeoBEARS")
setwd(extdata_dir)

# LAGRANGE run
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object
BioGeoBEARS_model_object=BioGeoBEARS_run_object$BioGeoBEARS_model_object

run1 = bears_optim_run(BioGeoBEARS_run_object=BioGeoBEARS_run_object)

# DECj run
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.1
run2 = bears_optim_run(BioGeoBEARS_run_object=BioGeoBEARS_run_object)

run1$optim_result
run2$optim_result

run1$output$BioGeoBEARS_model_object@params_table
run2$output$BioGeoBEARS_model_object@params_table


# DECjv run
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.1
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.1
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv-v"

run3 = bears_optim_run(BioGeoBEARS_run_object=BioGeoBEARS_run_object)

run1$optim_result
run2$optim_result
run3$optim_result

run1$output$BioGeoBEARS_model_object@params_table
run2$output$BioGeoBEARS_model_object@params_table
run3$output$BioGeoBEARS_model_object@params_table





#######################################################
# Do some simulations
#######################################################
# Simulate under run1 model
BioGeoBEARS_model_object = run1$output$BioGeoBEARS_model_object


trfn="Psychotria_5.2.newick"
phy = read.tree(trfn)

# Get geographic ranges at tips
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=BioGeoBEARS_run_object$geogfn)


# Get the list of geographic areas
areas = getareas_from_tipranges_object(tipranges)
areas_list = seq(0, length(areas)-1, 1)		# 0-base indexes

# Change the names to tipranges@df:
names(tipranges@df) = areas_list

maxareas = length(areas_list)
states_list = rcpp_areas_list_to_states_list(areas=areas_list, maxareas=maxareas, include_null_range=TRUE)


# Set the dispersal and extinction rate
d = BioGeoBEARS_model_object@params_table["d","est"]
e = BioGeoBEARS_model_object@params_table["e","est"]

# Equal dispersal in all directions (unconstrained)
# Equal extinction probability for all areas
distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

dmat = matrix(d, nrow=length(areas), ncol=length(areas))
elist = rep(e, length(areas))

# Set up the instantaneous rate matrix (Q matrix)
force_sparse = FALSE
Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)

# Cladogenic model
j = BioGeoBEARS_model_object@params_table["j","est"]
ysv = BioGeoBEARS_model_object@params_table["ysv","est"]
v = BioGeoBEARS_model_object@params_table["v","est"]
ys = BioGeoBEARS_model_object@params_table["ys","est"]
sum_SPweights = ys + j + v

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
states_indices = states_list
states_indices[1] = NULL	# shorten the states_indices by 1 (cutting the null range state from the speciation matrix)
spPmat_inputs$l = states_indices
spPmat_inputs$s = ys
spPmat_inputs$v = v
spPmat_inputs$j = j
spPmat_inputs$y = ys
spPmat_inputs$dmat = distances_mat
spPmat_inputs$maxent01s_param = maxent01s_param
spPmat_inputs$maxent01v_param = maxent01v_param
spPmat_inputs$maxent01j_param = maxent01j_param
spPmat_inputs$maxent01y_param = maxent01y_param

			l = spPmat_inputs$l		# states_indices
			s = spPmat_inputs$s
			v = spPmat_inputs$v
			j = spPmat_inputs$j
			y = spPmat_inputs$y
			
			dmat = spPmat_inputs$dmat
			
			# Take the max of the indices of the possible areas, and add 1
			numareas = max(unlist(spPmat_inputs$l), na.rm=TRUE) + 1
			
			maxent01s_param = spPmat_inputs$maxent01s_param
			maxent01v_param = spPmat_inputs$maxent01v_param
			maxent01j_param = spPmat_inputs$maxent01j_param
			maxent01y_param = spPmat_inputs$maxent01y_param
			
			maxent01s = relative_probabilities_of_subsets(max_numareas=numareas, maxent_constraint_01=maxent01s_param, NA_val=0)
			maxent01v = relative_probabilities_of_vicariants(max_numareas=numareas, maxent_constraint_01v=maxent01v_param, NA_val=0)
			maxent01j = relative_probabilities_of_subsets(max_numareas=numareas, maxent_constraint_01=maxent01j_param, NA_val=0)
			maxent01y = relative_probabilities_of_subsets(max_numareas=numareas, maxent_constraint_01=maxent01y_param, NA_val=0)

			# You really need a list of sizes here:
			
			# Matrix of probs for each ancsize
			maxprob_as_function_of_ancsize_and_decsize = mapply(FUN=max, maxent01s, maxent01v, maxent01j, maxent01y, MoreArgs=list(na.rm=TRUE))
			maxprob_as_function_of_ancsize_and_decsize = matrix(data=maxprob_as_function_of_ancsize_and_decsize, nrow=nrow(maxent01s), ncol=ncol(maxent01s))
			maxprob_as_function_of_ancsize_and_decsize[maxprob_as_function_of_ancsize_and_decsize > 0] = 1
			maxprob_as_function_of_ancsize_and_decsize[maxprob_as_function_of_ancsize_and_decsize <= 0] = 0
			
			# Now, go through, and make a list of the max minsize for each decsize
			max_minsize_as_function_of_ancsize = apply(X=maxprob_as_function_of_ancsize_and_decsize, MARGIN=1, FUN=maxsize)


tmpca_1 = rep(1, length(states_list)-1)
tmpcb_1 = rep(1, length(states_list)-1)
COO_weights_columnar = rcpp_calc_anclikes_sp_COOweights_faster(Rcpp_leftprobs=tmpca_1, Rcpp_rightprobs=tmpcb_1, l=spPmat_inputs$l, s=s, v=v, j=j, y=y, dmat=dmat, maxent01s=maxent01s, maxent01v=maxent01v, maxent01j=maxent01j, maxent01y=maxent01y, max_minsize_as_function_of_ancsize=max_minsize_as_function_of_ancsize, printmat=FALSE)

index_Qmat_0based_of_starting_state = sample(1:(length(states_list)-1), size=1)
simulate_biogeog_history(phy, Qmat, COO_probs_columnar=COO_weights_columnar, index_Qmat_0based_of_starting_state)


LAGRANGE_inf = NULL
DECj_inf = NULL


MLstate_accuracy_2param = NULL
MLprobs_accuracy_2param = NULL
MLstate_accuracy_3param = NULL
MLprobs_accuracy_3param = NULL

for (i in 1:100)
	{
	LAGRANGE_sim = simulate_biogeog_history(phy, Qmat, COO_probs_columnar=COO_weights_columnar, index_Qmat_0based_of_starting_state)
	
	simgeog_fn = simulated_indexes_to_tipranges_file(simulated_states_by_node=LAGRANGE_sim, areas_list, states_list, trfn, out_geogfn="lagrange_area_data_file.data")
	
	# LAGRANGE run
	BioGeoBEARS_run_object = define_BioGeoBEARS_run()
	BioGeoBEARS_run_object$geogfn = simgeog_fn
	tmpinf2 = bears_optim_run(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
	siminf2 = tmpinf2$optim_result
	LAGRANGE_inf = rbind(LAGRANGE_inf, siminf2)

	# DECj run
	BioGeoBEARS_run_object = define_BioGeoBEARS_run()
	BioGeoBEARS_run_object$geogfn = simgeog_fn
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.1

	tmpinf3 = bears_optim_run(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
	siminf3 = tmpinf3$optim_result
	DECj_inf = rbind(DECj_inf, siminf3)
	

	
	
	# 3-parameter inference
	simstates = LAGRANGE_sim[(1+length(phy$tip.label)):(length(phy$tip.label)+phy$Nnode)]
	infstates = unlist(get_ML_states(relprobs_matrix=tmpinf3$ML_marginal_prob_each_state_at_branch_top_AT_node)[(1+length(phy$tip.label)):(length(phy$tip.label)+phy$Nnode)])
	simstates
	infstates
	
	matchTF = mapply(FUN="==", simstates, infstates)
	accuracy = sum(matchTF) / length(matchTF)
	MLstate_accuracy_3param = c(MLstate_accuracy_3param, accuracy)
	
	# Accuracy of ancprobs at nodes
	infprobs = unlist(get_ML_probs(relprobs_matrix=tmpinf3$ML_marginal_prob_each_state_at_branch_top_AT_node)[(1+length(phy$tip.label)):(length(phy$tip.label)+phy$Nnode)])
	infprobs
	
	prob_diffs = rep(0, phy$Nnode)
	# add the prob_diffs when same
	prob_diffs[matchTF == TRUE] = 1-infprobs[matchTF == TRUE]
	# add the prob_diffs when different
	prob_diffs[matchTF == FALSE] = infprobs[matchTF == FALSE]
	
	MLprobs_accuracy_3param = rbind(MLprobs_accuracy_3param, prob_diffs)



	# Accuracy of ancstates at nodes

	# 2-parameter inference
	simstates = LAGRANGE_sim[(1+length(phy$tip.label)):(length(phy$tip.label)+phy$Nnode)]
	infstates = unlist(get_ML_states(relprobs_matrix=tmpinf2$ML_marginal_prob_each_state_at_branch_top_AT_node)[(1+length(phy$tip.label)):(length(phy$tip.label)+phy$Nnode)])
	simstates
	infstates
	
	matchTF = mapply(FUN="==", simstates, infstates)
	accuracy = sum(matchTF) / length(matchTF)
	MLstate_accuracy_2param = c(MLstate_accuracy_2param, accuracy)
	
	# Accuracy of ancprobs at nodes
	infprobs = unlist(get_ML_probs(relprobs_matrix=tmpinf2$ML_marginal_prob_each_state_at_branch_top_AT_node)[(1+length(phy$tip.label)):(length(phy$tip.label)+phy$Nnode)])
	infprobs
	
	prob_diffs = rep(0, phy$Nnode)
	# add the prob_diffs when same
	prob_diffs[matchTF == TRUE] = 1-infprobs[matchTF == TRUE]
	# add the prob_diffs when different
	prob_diffs[matchTF == FALSE] = infprobs[matchTF == FALSE]
	
	MLprobs_accuracy_2param = rbind(MLprobs_accuracy_2param, prob_diffs)

	}


simp_MLstate_accuracy_2param = MLstate_accuracy_2param
simp_MLprobs_accuracy_2param = MLprobs_accuracy_2param
simp_MLstate_accuracy_3param = MLstate_accuracy_3param
simp_MLprobs_accuracy_3param = MLprobs_accuracy_3param



LAGRANGE_inf
DECj_inf










#######################################################
# Do some simulations
#######################################################
# Simulate under run2 model
BioGeoBEARS_model_object = run2$output$BioGeoBEARS_model_object


trfn="Psychotria_5.2.newick"
phy = read.tree(trfn)

# Get geographic ranges at tips
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=BioGeoBEARS_run_object$geogfn)


# Get the list of geographic areas
areas = getareas_from_tipranges_object(tipranges)
areas_list = seq(0, length(areas)-1, 1)		# 0-base indexes

# Change the names to tipranges@df:
names(tipranges@df) = areas_list

maxareas = length(areas_list)
states_list = rcpp_areas_list_to_states_list(areas=areas_list, maxareas=maxareas, include_null_range=TRUE)


# Set the dispersal and extinction rate
d = BioGeoBEARS_model_object@params_table["d","est"]
e = BioGeoBEARS_model_object@params_table["e","est"]

# Equal dispersal in all directions (unconstrained)
# Equal extinction probability for all areas
distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

dmat = matrix(d, nrow=length(areas), ncol=length(areas))
elist = rep(e, length(areas))

# Set up the instantaneous rate matrix (Q matrix)
force_sparse = FALSE
Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)

# Cladogenic model
j = BioGeoBEARS_model_object@params_table["j","est"]
ysv = BioGeoBEARS_model_object@params_table["ysv","est"]
v = BioGeoBEARS_model_object@params_table["v","est"]
ys = BioGeoBEARS_model_object@params_table["ys","est"]
sum_SPweights = ys + j + v

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
states_indices = states_list
states_indices[1] = NULL	# shorten the states_indices by 1 (cutting the null range state from the speciation matrix)
spPmat_inputs$l = states_indices
spPmat_inputs$s = ys
spPmat_inputs$v = v
spPmat_inputs$j = j
spPmat_inputs$y = ys
spPmat_inputs$dmat = distances_mat
spPmat_inputs$maxent01s_param = maxent01s_param
spPmat_inputs$maxent01v_param = maxent01v_param
spPmat_inputs$maxent01j_param = maxent01j_param
spPmat_inputs$maxent01y_param = maxent01y_param

			l = spPmat_inputs$l		# states_indices
			s = spPmat_inputs$s
			v = spPmat_inputs$v
			j = spPmat_inputs$j
			y = spPmat_inputs$y
			
			dmat = spPmat_inputs$dmat
			
			# Take the max of the indices of the possible areas, and add 1
			numareas = max(unlist(spPmat_inputs$l), na.rm=TRUE) + 1
			
			maxent01s_param = spPmat_inputs$maxent01s_param
			maxent01v_param = spPmat_inputs$maxent01v_param
			maxent01j_param = spPmat_inputs$maxent01j_param
			maxent01y_param = spPmat_inputs$maxent01y_param
			
			maxent01s = relative_probabilities_of_subsets(max_numareas=numareas, maxent_constraint_01=maxent01s_param, NA_val=0)
			maxent01v = relative_probabilities_of_vicariants(max_numareas=numareas, maxent_constraint_01v=maxent01v_param, NA_val=0)
			maxent01j = relative_probabilities_of_subsets(max_numareas=numareas, maxent_constraint_01=maxent01j_param, NA_val=0)
			maxent01y = relative_probabilities_of_subsets(max_numareas=numareas, maxent_constraint_01=maxent01y_param, NA_val=0)

			# You really need a list of sizes here:
			
			# Matrix of probs for each ancsize
			maxprob_as_function_of_ancsize_and_decsize = mapply(FUN=max, maxent01s, maxent01v, maxent01j, maxent01y, MoreArgs=list(na.rm=TRUE))
			maxprob_as_function_of_ancsize_and_decsize = matrix(data=maxprob_as_function_of_ancsize_and_decsize, nrow=nrow(maxent01s), ncol=ncol(maxent01s))
			maxprob_as_function_of_ancsize_and_decsize[maxprob_as_function_of_ancsize_and_decsize > 0] = 1
			maxprob_as_function_of_ancsize_and_decsize[maxprob_as_function_of_ancsize_and_decsize <= 0] = 0
			
			# Now, go through, and make a list of the max minsize for each decsize
			max_minsize_as_function_of_ancsize = apply(X=maxprob_as_function_of_ancsize_and_decsize, MARGIN=1, FUN=maxsize)


tmpca_1 = rep(1, length(states_list)-1)
tmpcb_1 = rep(1, length(states_list)-1)
COO_weights_columnar = rcpp_calc_anclikes_sp_COOweights_faster(Rcpp_leftprobs=tmpca_1, Rcpp_rightprobs=tmpcb_1, l=spPmat_inputs$l, s=s, v=v, j=j, y=y, dmat=dmat, maxent01s=maxent01s, maxent01v=maxent01v, maxent01j=maxent01j, maxent01y=maxent01y, max_minsize_as_function_of_ancsize=max_minsize_as_function_of_ancsize, printmat=FALSE)

index_Qmat_0based_of_starting_state = sample(1:(length(states_list)-1), size=1)
simulate_biogeog_history(phy, Qmat, COO_probs_columnar=COO_weights_columnar, index_Qmat_0based_of_starting_state)


set.seed(12345)

LAGRANGE_inf_3param = NULL
DECj_inf_3param = NULL

MLstate_accuracy_2param = NULL
MLprobs_accuracy_2param = NULL
MLstate_accuracy_3param = NULL
MLprobs_accuracy_3param = NULL

for (i in 1:100)
	{
	LAGRANGE_sim = simulate_biogeog_history(phy, Qmat, COO_probs_columnar=COO_weights_columnar, index_Qmat_0based_of_starting_state)
	
	simgeog_fn = simulated_indexes_to_tipranges_file(simulated_states_by_node=LAGRANGE_sim, areas_list, states_list, trfn, out_geogfn="lagrange_area_data_file.data")
	
	# LAGRANGE run
	BioGeoBEARS_run_object = define_BioGeoBEARS_run()
	BioGeoBEARS_run_object$geogfn = simgeog_fn
	tmpinf2 = bears_optim_run(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
	siminf2 = tmpinf2$optim_result
	LAGRANGE_inf_3param = rbind(LAGRANGE_inf_3param, siminf2)

	# DECj run
	BioGeoBEARS_run_object = define_BioGeoBEARS_run()
	BioGeoBEARS_run_object$geogfn = simgeog_fn
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.1

	tmpinf3 = bears_optim_run(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
	siminf3 = tmpinf3$optim_result
	DECj_inf_3param = rbind(DECj_inf_3param, siminf3)
	

	
	
	# 3-parameter inference
	simstates = LAGRANGE_sim[(1+length(phy$tip.label)):(length(phy$tip.label)+phy$Nnode)]
	infstates = unlist(get_ML_states(relprobs_matrix=tmpinf3$ML_marginal_prob_each_state_at_branch_top_AT_node)[(1+length(phy$tip.label)):(length(phy$tip.label)+phy$Nnode)])
	simstates
	infstates
	
	matchTF = mapply(FUN="==", simstates, infstates)
	accuracy = sum(matchTF) / length(matchTF)
	MLstate_accuracy_3param = c(MLstate_accuracy_3param, accuracy)
	
	# Accuracy of ancprobs at nodes
	infprobs = unlist(get_ML_probs(relprobs_matrix=tmpinf3$ML_marginal_prob_each_state_at_branch_top_AT_node)[(1+length(phy$tip.label)):(length(phy$tip.label)+phy$Nnode)])
	infprobs
	
	prob_diffs = rep(0, phy$Nnode)
	# add the prob_diffs when same
	prob_diffs[matchTF == TRUE] = 1-infprobs[matchTF == TRUE]
	# add the prob_diffs when different
	prob_diffs[matchTF == FALSE] = infprobs[matchTF == FALSE]
	
	MLprobs_accuracy_3param = rbind(MLprobs_accuracy_3param, prob_diffs)



	# Accuracy of ancstates at nodes

	# 2-parameter inference
	simstates = LAGRANGE_sim[(1+length(phy$tip.label)):(length(phy$tip.label)+phy$Nnode)]
	infstates = unlist(get_ML_states(relprobs_matrix=tmpinf2$ML_marginal_prob_each_state_at_branch_top_AT_node)[(1+length(phy$tip.label)):(length(phy$tip.label)+phy$Nnode)])
	simstates
	infstates
	
	matchTF = mapply(FUN="==", simstates, infstates)
	accuracy = sum(matchTF) / length(matchTF)
	MLstate_accuracy_2param = c(MLstate_accuracy_2param, accuracy)
	
	# Accuracy of ancprobs at nodes
	infprobs = unlist(get_ML_probs(relprobs_matrix=tmpinf2$ML_marginal_prob_each_state_at_branch_top_AT_node)[(1+length(phy$tip.label)):(length(phy$tip.label)+phy$Nnode)])
	infprobs
	
	prob_diffs = rep(0, phy$Nnode)
	# add the prob_diffs when same
	prob_diffs[matchTF == TRUE] = 1-infprobs[matchTF == TRUE]
	# add the prob_diffs when different
	prob_diffs[matchTF == FALSE] = infprobs[matchTF == FALSE]
	
	MLprobs_accuracy_2param = rbind(MLprobs_accuracy_2param, prob_diffs)
	
	
	}


MLstate_accuracy_2param
MLstate_accuracy_3param

simp_MLstate_accuracy_2param
simp_MLstate_accuracy_3param


mean(MLstate_accuracy_2param)
mean(MLstate_accuracy_3param)

mean(simp_MLstate_accuracy_2param)
mean(simp_MLstate_accuracy_3param)



MLprobs_accuracy_2param
MLprobs_accuracy_3param


mean(MLprobs_accuracy_2param)
mean(MLprobs_accuracy_3param)

mean(simp_MLprobs_accuracy_2param)
mean(simp_MLprobs_accuracy_3param)


LAGRANGE_inf
DECj_inf
LAGRANGE_inf_3param
DECj_inf_3param

dim(LAGRANGE_inf)
dim(DECj_inf)
dim(LAGRANGE_inf_3param)
dim(DECj_inf_3param)



mean(as.numeric(LAGRANGE_inf$fvalues))
mean(as.numeric(DECj_inf$fvalues))
mean(as.numeric(LAGRANGE_inf_3param$fvalues))
mean(as.numeric(DECj_inf_3param$fvalues))


# Mean difference
LnL_diff_sim2param = as.numeric(DECj_inf$fvalues) - as.numeric(LAGRANGE_inf$fvalues)
LnL_diff_sim3param = as.numeric(DECj_inf_3param$fvalues) - as.numeric(LAGRANGE_inf_3param$fvalues)

mean(LnL_diff_sim2param)
mean(LnL_diff_sim3param)

par(mfrow=c(2,1))
hist(LnL_diff_sim2param, breaks=50, xlim=c(0,20), col="red", main="simulation under 2-parameter model", xlab="log-likelihood advantage of 3-parameter inference")
hist(LnL_diff_sim3param, breaks=50, xlim=c(0,20), col="blue", main="simulation under 3-parameter model", xlab="log-likelihood advantage of 3-parameter inference")


# Likelihood ratio test
par(mfrow=c(2,1))
pval_cutoff = 0.05

LRT_results_2param = mapply(FUN=lrttest, LnL_1=as.numeric(DECj_inf$fvalues), LnL_2=as.numeric(LAGRANGE_inf$fvalues), numparams1=3, numparams2=2, returnwhat="pval")
LRT_results_2param
hist(LRT_results_2param, breaks=50, xlim=c(0,1), col="red", main="Likelihood ratio test, simulation under 2-parameter model", xlab="p-value")

# false positives (rejecting null hypothesis that the two models fit the same)
false_positive_rate = sum(LRT_results_2param <= pval_cutoff)
false_positive_rate
false_positive_rate = false_positive_rate / length(LRT_results_2param)
false_positive_rate


LRT_results_3param = mapply(FUN=lrttest, LnL_1=as.numeric(DECj_inf_3param$fvalues), LnL_2=as.numeric(LAGRANGE_inf_3param$fvalues), numparams1=3, numparams2=2, returnwhat="pval")
LRT_results_3param
hist(LRT_results_3param, breaks=50, xlim=c(0,1), col="blue", main="Likelihood ratio test, simulation under 3-parameter model", xlab="p-value")

# false_negative_rate (failing to reject the null hypothesis, when the two models DONT fit the same)
false_negative_rate = sum(LRT_results_3param > pval_cutoff)
false_negative_rate
false_negative_rate = false_negative_rate / length(LRT_results_3param)
false_negative_rate





' # End junk

