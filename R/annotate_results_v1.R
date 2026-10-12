
#######################################################
# Annotate BioGeoBEARS results
#######################################################
cft_ancstates_res <- function(res, sepval="", numbers_below_this_get_scientific=0.0001, numdigits_for_superlow_scientific=1, numbers_above_this_get_scientific=10000000, numdigits_for_superhigh_scientific=2, numdigits_inbetween_have_fixed_digits=4, digits=NULL)
	{
	setup='
	BioGeoBEARS_run_object = define_BioGeoBEARS_run()
	res = bears_optim_run(BioGeoBEARS_run_object=BioGeoBEARS_run_object, skip_optim=TRUE, skip_optim_option="return_all")
	sepval=""	
	numbers_below_this_get_scientific=0.0001
	numdigits_for_superlow_scientific=1
	numbers_above_this_get_scientific=10000000
	numdigits_for_superhigh_scientific=2
	numdigits_inbetween_have_fixed_digits=4
	digits=NULL

	ancstates_df = cft_ancstates_res(res, sepval=sepval, numbers_below_this_get_scientific=numbers_below_this_get_scientific, numdigits_for_superlow_scientific=numdigits_for_superlow_scientific, numbers_above_this_get_scientific=numbers_above_this_get_scientific, numdigits_for_superhigh_scientific=numdigits_for_superhigh_scientific, numdigits_inbetween_have_fixed_digits=numdigits_inbetween_have_fixed_digits, digits=digits)
	'
	
	res2 = annotate_BGB_results(res, sepval=sepval)
	ancstates = res2$ML_marginal_prob_each_state_at_branch_top_AT_node
	col_headers = colnames(ancstates)
	ancstates_df = cft(input_table=ancstates, numbers_below_this_get_scientific=numbers_below_this_get_scientific, numdigits_for_superlow_scientific=numdigits_for_superlow_scientific, numbers_above_this_get_scientific=numbers_above_this_get_scientific, numdigits_for_superhigh_scientific=numdigits_for_superhigh_scientific, numdigits_inbetween_have_fixed_digits=numdigits_inbetween_have_fixed_digits, digits=digits, list_of_colnames=col_headers)
	return(ancstates_df)
	} # END cft_ancstates_res <- function(res, sepval="")


annotate_BGB_results <- function(res, sepval="")
	{
	setup='
	BioGeoBEARS_run_object = define_BioGeoBEARS_run()
	res = bears_optim_run(BioGeoBEARS_run_object=BioGeoBEARS_run_object, skip_optim=TRUE, skip_optim_option="return_all")
	sepval=""
	ranges_list = get_ranges_list_from_BGB_results_object(res, sepval=sepval)
	res2 = annotate_BGB_results(res, sepval=sepval)
	'
	
	
	nodenums = 1:length(res$computed_likelihoods_at_each_node)
	ranges_list = get_ranges_list_from_BGB_results_object(res, sepval=sepval)
	res2 = res
	
	names(res2$computed_likelihoods_at_each_node) = nodenums
	names(res2$relative_probs_of_each_state_at_bottom_of_root_branch) = ranges_list

	colnames(res2$condlikes_of_each_state) = ranges_list
	row.names(res2$condlikes_of_each_state) = nodenums

	colnames(res2$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS) = ranges_list
	row.names(res2$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS) = nodenums

	colnames(res2$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS) = ranges_list
	row.names(res2$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS) = nodenums

	colnames(res2$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS) = ranges_list
	row.names(res2$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS) = nodenums

	colnames(res2$relative_probs_of_each_state_at_branch_top_AT_node_UPPASS) = ranges_list
	row.names(res2$relative_probs_of_each_state_at_branch_top_AT_node_UPPASS) = nodenums

	colnames(res2$ML_marginal_prob_each_state_at_branch_bottom_below_node) = ranges_list
	row.names(res2$ML_marginal_prob_each_state_at_branch_bottom_below_node) = nodenums

	colnames(res2$ML_marginal_prob_each_state_at_branch_top_AT_node) = ranges_list
	row.names(res2$ML_marginal_prob_each_state_at_branch_top_AT_node) = nodenums

	return(res2)	
	} # END  <- function(res, sepval="")

# Make a text-based list of ranges, with separator sep=""
get_ranges_list_from_BGB_run_object <- function(BioGeoBEARS_run_object=define_BioGeoBEARS_run(), sepval="")
	{
	setup='
	BioGeoBEARS_run_object=define_BioGeoBEARS_run()
	sepval=""
	ranges_list = get_ranges_list_from_BGB_run_object(BioGeoBEARS_run_object=BioGeoBEARS_run_object, sepval=sepval)
	'
	
	txt = paste0("get_ranges_list_from_BGB_run_object(), in order to get a res$inputs object with\nres$inputs$all_geog_states_list_usually_inferred_from_areas_maxareas, is running:\nres = bears_optim_run(BioGeoBEARS_run_object=BioGeoBEARS_run_object, skip_optim=TRUE, skip_optim_option='return_all')")
	cat("\n")
	cat(txt)
	cat("\n")
	
	res = bears_optim_run(BioGeoBEARS_run_object=BioGeoBEARS_run_object, skip_optim=TRUE, skip_optim_option="return_all")

	ranges_list = get_ranges_list_from_BGB_results_object <- function(res, sepval=sepval)
	
	return(ranges_list)
	} # END get_ranges_list_from_BGB_run_object <- function(BioGeoBEARS_run_object=define_BioGeoBEARS_run(), sepval="")


get_ranges_list_from_BGB_results_object <- function(res, sepval="")
	{
	setup='
	BioGeoBEARS_run_object = define_BioGeoBEARS_run()
	res = bears_optim_run(BioGeoBEARS_run_object=BioGeoBEARS_run_object, skip_optim=TRUE, skip_optim_option="return_all")
	ranges_list = get_ranges_list_from_BGB_results_object(res, sepval="")
	
	'
	
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=res$inputs$geogfn)
	tipranges
	
	area_names = names(tipranges@df)
	area_names
	
	all_states_areanums01 = res$inputs$all_geog_states_list_usually_inferred_from_areas_maxareas
	all_states_areanums01
	
	# Make the list of ranges
	ranges_list = NULL
	for (i in 1:length(all_states_areanums01))
			{    
			if ( (length(all_states_areanums01[[i]]) == 1) && (is.na(all_states_areanums01[[i]])) )
					{
					tmprange = "_"
					} else {
					tmprange = paste0(area_names[all_states_areanums01[[i]]+1], collapse=sepval)
					}
			ranges_list = c(ranges_list, tmprange)
			}
	
	# Look at the ranges list
	ranges_list
	
	return(ranges_list)	
	} # END get_ranges_list_from_BGB_results_object


