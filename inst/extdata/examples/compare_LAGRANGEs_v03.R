
#######################################################
# Source code
#######################################################
library(LaplacesDemon)
library(BioGeoBEARS)
source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_basics_v1.R', chdir = TRUE)
source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_generics_v1.R', chdir = TRUE)
source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_classes_v1.R', chdir = TRUE)
source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_univ_model_v1.R', chdir = TRUE)
source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_readwrite_v1.R', chdir = TRUE)
source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_models_v1.R', chdir = TRUE)
source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/calc_loglike_sp_v01.R', chdir = TRUE)
source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_LaplacesDemon_v1.R', chdir = TRUE)
#source('/Dropbox/_njm/_biogeog_sim_utils_v1.R', chdir = TRUE)


#######################################################
# File locations
#######################################################
# In package
extdata_dir = system.file("extdata", package="BioGeoBEARS")
# development
extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"


#######################################################
# Input files (tree and tipdata)
#######################################################
trfn = paste(addslash(extdata_dir), "Psychotria_5.2.newick", sep="")
tr = read.tree(trfn)
geogfn = paste(addslash(extdata_dir), "Psychotria_geog.data", sep="")
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges = order_tipranges_by_tree_tips(tipranges, tr)
tiprange_names = tipranges_to_area_strings(tipranges=tipranges)
tiprange_names

#######################################################
# Reading Python LAGRANGE results
#######################################################
wd = slashslash(paste(extdata_dir, "/examples/Psychotria_M0/LAGRANGE_python/", sep=""))
setwd(wd)
outfn = "Psychotria_5.2_demo.results.txt"
results_dir = wd
new_splits_fn = TRUE
new_states_fn = FALSE
filecount=0

summstats = parse_lagrange_python_output(outfn, results_dir, new_splits_fn = TRUE, new_states_fn = FALSE, filecount=0)
summstats


#######################################################
# OK now read the Python LAGRANGE version
#######################################################
splits_fn = "Psychotria_5.2_demo.results_splits00001.txt"
moref(splits_fn)

splits_LGpy = LGpy_splits_fn_to_table(splits_fn)
MLsplits_LGpy = LGpy_MLsplit_per_node(splits_LGpy)
#MLsplits_LGpy = map_LGpy_MLsplits_to_tree(MLsplits_LGpy, tr, tiprange_names)
MLsplits_LGpy = map_LG_MLsplits_to_tree(MLsplits_LGpy, tr, tiprange_names, type="python")


#######################################################
# OK now read the C++ LAGRANGE version
#######################################################
# wd = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/Psychotria_M0/LAGRANGE_python/"
# setwd(wd)
outfn = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/Psychotria_M0/LAGRANGE_C++/Psychotria_M0_lgcpp_out.txt"
results_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/Psychotria_M0/LAGRANGE_C++/"
new_splits_fn = TRUE
new_states_fn = TRUE
filecount=0

summstats = parse_lagrange_output(outfn, results_dir, new_splits_fn, new_states_fn, filecount)
summstats

splits_fn = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/Psychotria_M0/LAGRANGE_C++/Psychotria_M0_lgcpp_out_splits00001.txt"

splits_LGcpp = LGcpp_splits_fn_to_table(splits_fn)
MLsplits_LGcpp = LGpy_MLsplit_per_node(splits_LGcpp)
MLsplits_LGcpp = map_LG_MLsplits_to_tree(MLsplits_LGcpp, tr, tiprange_names, removechar="_", type="C++")
MLsplits_LGcpp








#######################################################
# Put the splits on the corners
#######################################################
sourcedir = '/Dropbox/_njm/'
source3 = '_R_tree_functions_v1.R'
source(paste(sourcedir, source3, sep=""))

# Get corner coordinates
corners_list = corner_coords(tr)
leftcorns = corners_list$leftcorns
rightcorns = corners_list$rightcorns

# Plot splits on corners
# Ensure correct order
MLsplits_LGcpp = order_LGnodes(MLsplits_LGcpp, removechar="_", type="C++")
MLsplits_LGcpp

plot(tr, label.offset=0.15)
cornerlabels(text=MLsplits_LGcpp$leftBB, coords=leftcorns, bg="green3")
cornerlabels(text=MLsplits_LGcpp$rightBB, coords=rightcorns, bg="green3")
tiplabels(tiprange_names)
axisPhylo()
mtext(text="million years ago", side=1, line=2)









#######################################################
# Get the states from LAGRANGE C++
#######################################################
source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_readwrite_v1.R', chdir = TRUE)

# wd = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/Psychotria_M0/LAGRANGE_python/"
# setwd(wd)
outfn = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/Psychotria_M0/LAGRANGE_C++/Psychotria_M0_lgcpp_out.txt"
results_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/Psychotria_M0/LAGRANGE_C++/"
new_splits_fn = TRUE
new_states_fn = TRUE
filecount=0

summstats = parse_lagrange_output(outfn, results_dir, new_splits_fn, new_states_fn, filecount)
summstats

states_fn = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/Psychotria_M0/LAGRANGE_C++/Psychotria_M0_lgcpp_out_states00001.txt"

states_LGcpp = LGcpp_states_fn_to_table(states_fn)
MLstates_LGcpp = LGcpp_MLstate_per_node(states_LGcpp)
order_LGnodes(MLstates_LGcpp, removechar=NULL, type="C++", type2="states")
MLstates_LGcpp = map_LG_MLstates_to_tree(MLstates_LGcpp, tr, tiprange_names, removechar="_", type="C++")
MLstates_LGcpp





#######################################################
# Compare to states from BioGeoBEARS 2-parameter analysis
#######################################################
bears_output = bears_2param_standard_fast(trfn=trfn, geogfn=geogfn)
bears_output


# Extract key parameters and LnL values
fn_name = "bears_2param_standard_fast"
tmp_params = get_infparams_optimx_nosim(bears_output, fn_name)
tmp_LnL = get_inf_LgL_etc_optimx(bears_output)
results_row = cbind(fn_name, t(tmp_params), tmp_LnL)
results_row



#######################################################
# Make plots
#######################################################

# Get ranges from geogfn
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

# Get the default areas
areas = getareas_from_tipranges_object(tipranges)
areas

# Name the areas (Kauai, Oahu, Maui-Nui, Hawaii)
areas = c("K","O","M","H")

# Make states and state indices
states = areas_list_to_states_list_old(areas=areas)
states
states2 = rcpp_areas_list_to_states_list(areas=areas)
states2






#######################################################
#######################################################
# Plot the various ancestral state estimates
#######################################################
#######################################################
piecols = c("black", 
rgb(blue=100, green=0, red=0, maxColorValue=100), 		# blue
rgb(blue=100, green=100, red=0, maxColorValue=100), 	# cyan
rgb(blue=0, green=100, red=100, maxColorValue=100), 	# yellow
rgb(blue=0, green=0, red=100, maxColorValue=100), 		# red

rgb(blue=100, green=50, red=0, maxColorValue=100), 
rgb(blue=75, green=50, red=50, maxColorValue=100), 
rgb(blue=75, green=0, red=50, maxColorValue=100), 

rgb(blue=0, green=100, red=50, maxColorValue=100), 
rgb(blue=0, green=50, red=50, maxColorValue=100), 

rgb(blue=0, green=50, red=100, maxColorValue=100), 


rgb(blue=50, green=25, red=25, maxColorValue=100), 
rgb(blue=50, green=50, red=25, maxColorValue=100), 
rgb(blue=25, green=25, red=50, maxColorValue=100), 

rgb(blue=25, green=50, red=50, maxColorValue=100), 

rgb(blue=100, green=100, red=100, maxColorValue=100))


doPDFs = FALSE
pdffn = "compare_anc_range_reconstructions_v01.pdf"
if (doPDFs) 
	{
	pdf(pdffn, height=11, width=8)
	}


par(mfrow=c(2,1))




###########################################################
# 2-parameter model (as in LAGRANGE)
# (Ree & Smith 2008)
###########################################################
results_object = bears_output
relprobs_matrix = results_object$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS
MLstates_0based_indexes = get_ML_states(relprobs_matrix)

MLstates = list()
for (i in 1:length(MLstates_0based_indexes))
	{
	MLstates[[i]] = states[[ MLstates_0based_indexes[[i]]+1 ]]
	}

MLstates

MLstates_txt = unlist(lapply(X=MLstates, FUN=paste, sep="", collapse=""))
MLstates_txt

ntips = length(tr$tip.label)
MLstates_txt_tips = MLstates_txt[1:ntips]
MLstates_txt_nodes = MLstates_txt[(1+ntips):(length(MLstates_txt))]
MLstates_txt_nodes


tr2 = reorder(tr, "pruningwise")
plot(tr2, label.offset=0.15)
nodelabels(MLstates_txt_nodes)
tiplabels(MLstates_txt_tips)
title("Best states based on downpass cond. likes. at branch tops with\nLAGRANGE 2-parameter model (not real inference)")

# Pie charts
MLstates_probs = get_ML_probs(relprobs_matrix=results_object$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS)
MLstates_probs = unlist(MLstates_probs)
MLstates_probs_other = 1-MLstates_probs
probs = cbind(MLstates_probs, MLstates_probs_other)

plot(tr2, label.offset=0.15)
#tiplabels(pie=relprobs_matrix[1:ntips,], piecol=piecols)
#nodelabels(pie=relprobs_matrix[(1+ntips):(length(MLstates_txt)),], piecol=piecols)
tiplabels(pie=probs[1:ntips,])
nodelabels(pie=probs[(1+ntips):(length(MLstates_txt)),])





###########################################################
# 2-parameter model (as in LAGRANGE)
# (Ree & Smith 2008)
###########################################################
results_object = bears_output
relprobs_matrix = results_object$ML_marginal_prob_each_state_at_branch_top_AT_node
MLstates_0based_indexes = get_ML_states(relprobs_matrix)

MLstates = list()
for (i in 1:length(MLstates_0based_indexes))
	{
	MLstates[[i]] = states[[ MLstates_0based_indexes[[i]]+1 ]]
	}

MLstates

MLstates_txt = unlist(lapply(X=MLstates, FUN=paste, sep="", collapse=""))
MLstates_txt

ntips = length(tr$tip.label)
MLstates_txt_tips = MLstates_txt[1:ntips]
MLstates_txt_nodes = MLstates_txt[(1+ntips):(length(MLstates_txt))]
MLstates_txt_nodes


tr2 = reorder(tr, "pruningwise")
plot(tr2, label.offset=0.15)
nodelabels(MLstates_txt_nodes)
tiplabels(MLstates_txt_tips)
title("ML marginal ancestral ranges at branch tops\nwith LAGRANGE 2-parameter model")

# Pie charts
MLstates_probs = get_ML_probs(relprobs_matrix=results_object$ML_marginal_prob_each_state_at_branch_top_AT_node)
MLstates_probs = unlist(MLstates_probs)
MLstates_probs_other = 1-MLstates_probs
probs = cbind(MLstates_probs, MLstates_probs_other)

plot(tr2, label.offset=0.15)
#tiplabels(pie=relprobs_matrix[1:ntips,], piecol=piecols)
#nodelabels(pie=relprobs_matrix[(1+ntips):(length(MLstates_txt)),], piecol=piecols)
tiplabels(pie=probs[1:ntips,])
nodelabels(pie=probs[(1+ntips):(length(MLstates_txt)),])


#######################################################
# Plot LAGRANGE C++ again
#######################################################
map_LG_MLstates_to_tree(MLstates_LGcpp, tr, tiprange_names, removechar="_", type="C++")






# Compare probs:
MLstates_LGcpp_justMLprob = MLstates_LGcpp$relprob[order(MLstates_LGcpp$Rnodes)]
nums = length(tr$tip.label)+1 : tr$Nnode
MLstates_bears_2param = probs[nums,1]

#plot(MLstates_LGcpp_justMLprob, MLstates_bears, xlim=c(0,1), ylim=c(0,1))










#######################################################
# Let's try joint splits/states
#######################################################
# Get ranges from geogfn
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

# Get the default areas
areas = getareas_from_tipranges_object(tipranges)
areas

# Name the areas (Kauai, Oahu, Maui-Nui, Hawaii)
areas = c("K","O","M","H")

# Make states and state indices
states = areas_list_to_states_list_old(areas=areas)
states
states2 = rcpp_areas_list_to_states_list(areas=areas)
states2



# Go down tree and get the downpass marginal probs above
# each node; then multiply to get the downpass marginal
# splits at each node
phy = tr

phy2 <- reorder(phy, "pruningwise")
numtips = length(phy2$tip.label)
tipnums <- 1:numtips
num_internal_nodes = phy2$Nnode
i = 1
edges_to_visit = seq(from=1, by=2, length.out=num_internal_nodes)
edges_to_visit

# Numstates during cladogenesis does NOT include the null range
numstates = numstates_from_numareas(numareas=length(areas), include_null_range=FALSE)
joint_uppass_likelihoods_of_splits_just_above_nodes = matrix(NA, nrow=num_internal_nodes, ncol=numstates*numstates)
joint_downpass_likelihoods_of_splits_just_above_nodes = matrix(NA, nrow=num_internal_nodes, ncol=numstates*numstates)
joint_marginal_splitprobs_at_internal_nodes = matrix(NA, nrow=num_internal_nodes, ncol=numstates*numstates)

# We can also do the ancstates, DEPENDENT on the splits
# Numstates during anagenesis DOES include the null range
numstates = numstates_from_numareas(numareas=length(areas), include_null_range=TRUE)
joint_uppass_likelihoods_of_states_dependent_on_splits_above = matrix(NA, nrow=num_internal_nodes, ncol=numstates)
joint_downpass_likelihoods_of_states_dependent_on_splits_above = matrix(NA, nrow=num_internal_nodes, ncol=numstates)
joint_marginal_states_dependent_on_splits_above = matrix(NA, nrow=num_internal_nodes, ncol=numstates)

dim(joint_splitprobs_at_internal_nodes)


#######################################################
#######################################################
# THIS IS A DOWNPASS FROM THE TIPS TO THE ROOT
#######################################################
#######################################################
for (i in edges_to_visit)
	{
	# First edge visited is i
	#print(i)
	
	# Its sister is j 
	j <- i + 1
	#print(j)

	# Get the node numbers at the tips of these two edges		
	left_desc_nodenum <- phy2$edge[i, 2]
	right_desc_nodenum <- phy2$edge[j, 2]
	
	# And for the ancestor edge (i or j shouldn't matter, should produce the same result!!!)
	anc <- phy2$edge[i, 1]
	
	# Here, we are just doing internal nodes, so
	anc = anc - numtips

	# Get the cladogenesis COO matrix
	# We can also do the ancstates, DEPENDENT on the splits
	# Using COOmat
	spPmat_inputs = bears_output$spPmat_inputs
	l = spPmat_inputs$l		# states_indices
	s = spPmat_inputs$s
	v = spPmat_inputs$v
	j = spPmat_inputs$j
	y = spPmat_inputs$y
	dmat = spPmat_inputs$dmat
	maxent01s = spPmat_inputs$maxent01s_param
	maxent01v = spPmat_inputs$maxent01v_param
	maxent01j = spPmat_inputs$maxent01j_param
	maxent01y = spPmat_inputs$maxent01y_param

	# Same either way for rowsums (just placeholders), but we need real probs 
	# for rcpp_calc_splitlikes_using_COOweights_columnar
	numstates = numstates_from_numareas(numareas=length(areas), include_null_range=FALSE)
	leftprobs = rep(1, numstates)
	rightprobs = rep(1, numstates)
	
	COO_weights_columnar = NULL
	
	COO_weights_columnar = rcpp_calc_anclikes_sp_COOweights_faster(Rcpp_leftprobs=leftprobs,
		Rcpp_rightprobs=rightprobs, l=l, s=s, v=v, j=j, y=y,
		dmat=dmat, maxent01s=maxent01s, maxent01v=maxent01v,
		maxent01j=maxent01j, maxent01y=maxent01y,
		max_minsize_as_function_of_ancsize = NULL,
		printmat = TRUE)
	COO_weights_columnar
	
	COO_weights_columnar_rowsums = rcpp_calc_rowsums_for_COOweights_columnar(COO_weights_columnar,
    numstates)
    
    
    



	
	# Get the relative probabilities (scaled conditional downpass likelihoods)
	# of each state on the downpass
	
	# add [-1] to not use null range in cladogenesis model
	ancnodes_left_daughter_downpass_relprob = bears_output$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[left_desc_nodenum,][-1]
	ancnodes_right_daughter_downpass_relprob = bears_output$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[right_desc_nodenum,][-1]

	# Do a fast cross-product with cladoRcpp:::rcpp_mult2probvect(a=ca, b=cb)
	require(cladoRcpp)
	tmp = rcpp_mult2probvect(a=ancnodes_left_daughter_downpass_relprob, b=ancnodes_right_daughter_downpass_relprob)
	# This does sum to 1
	
	joint_downpass_likelihoods_of_splits_just_above_nodes[anc,] = tmp

	# Same either way for rowsums (just placeholders), but we need real probs 
	# for rcpp_calc_splitlikes_using_COOweights_columnar
	leftprobs = ancnodes_left_daughter_downpass_relprob
	rightprobs = ancnodes_right_daughter_downpass_relprob

    stateprobs = rcpp_calc_splitlikes_using_COOweights_columnar(Rcpp_leftprobs=leftprobs,
    Rcpp_rightprobs=rightprobs, COO_weights_columnar, Rsp_rowsums=COO_weights_columnar_rowsums)
    stateprobs[is.nan(stateprobs)] = 0.0
    
    # Normalize & store
    stateprobs = stateprobs / sum(stateprobs, na.rm=TRUE)
    joint_downpass_likelihoods_of_states_dependent_on_splits_above[anc, ] = c(0,stateprobs)


	
	
	# Also gather the uppass likelihoods the same way
	ancnodes_left_daughter_uppass_relprob = bears_output$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS[left_desc_nodenum,][-1]
	ancnodes_right_daughter_uppass_relprob = bears_output$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS[right_desc_nodenum,][-1]

	# Do a fast cross-product with cladoRcpp:::rcpp_mult2probvect(a=ca, b=cb)
	require(cladoRcpp)
	tmp = rcpp_mult2probvect(a=ancnodes_left_daughter_uppass_relprob, b=ancnodes_right_daughter_uppass_relprob)
	# This does sum to 1
	
	joint_uppass_likelihoods_of_splits_just_above_nodes[anc,] = tmp
	
	# Same either way for rowsums (just placeholders), but we need real probs 
	# for rcpp_calc_splitlikes_using_COOweights_columnar
	leftprobs = ancnodes_left_daughter_uppass_relprob
	rightprobs = ancnodes_right_daughter_uppass_relprob

    stateprobs = rcpp_calc_splitlikes_using_COOweights_columnar(Rcpp_leftprobs=leftprobs,
    Rcpp_rightprobs=rightprobs, COO_weights_columnar, Rsp_rowsums=COO_weights_columnar_rowsums)
    stateprobs[is.nan(stateprobs)] = 0.0
    
    # Normalize & store
    stateprobs = stateprobs / sum(stateprobs, na.rm=TRUE)
    joint_uppass_likelihoods_of_states_dependent_on_splits_above[anc, ] = c(0,stateprobs)


	
	
	# And, the marginal is just the uppass * the downpass
	tmp = joint_downpass_likelihoods_of_splits_just_above_nodes[anc,] * joint_uppass_likelihoods_of_splits_just_above_nodes[anc,]
	tmp = tmp / sum(tmp, na.rm=TRUE)
	joint_marginal_splitprobs_at_internal_nodes[anc,] = tmp
	
	# For states also
	tmp = joint_downpass_likelihoods_of_states_dependent_on_splits_above[anc,] * joint_uppass_likelihoods_of_states_dependent_on_splits_above[anc,]
	tmp = tmp / sum(tmp, na.rm=TRUE)
	joint_marginal_states_dependent_on_splits_above[anc,] = tmp
    
    
	}

# Add the tip likelihoods, from e.g. the standard downpass likelihoods
joint_downpass_likelihoods_of_states_dependent_on_splits_above = rbind(bears_output$condlikes_of_each_state, joint_downpass_likelihoods_of_states_dependent_on_splits_above)

joint_uppass_likelihoods_of_states_dependent_on_splits_above = rbind(bears_output$condlikes_of_each_state, joint_uppass_likelihoods_of_states_dependent_on_splits_above)

joint_marginal_states_dependent_on_splits_above = rbind(bears_output$ML_marginal_prob_each_state_at_branch_top_AT_node, joint_marginal_states_dependent_on_splits_above)



###########################################################
# Plot splitprob-based states (marginal)
###########################################################
#results_object = bears_output
relprobs_matrix = joint_marginal_states_dependent_on_splits_above
MLstates_0based_indexes = get_ML_states(relprobs_matrix)

MLstates = list()
for (i in 1:length(MLstates_0based_indexes))
	{
	MLstates[[i]] = states[[ MLstates_0based_indexes[[i]]+1 ]]
	}

MLstates

MLstates_txt = unlist(lapply(X=MLstates, FUN=paste, sep="", collapse=""))
MLstates_txt

ntips = length(tr$tip.label)
MLstates_txt_tips = MLstates_txt[1:ntips]
MLstates_txt_nodes = MLstates_txt[(1+ntips):(length(MLstates_txt))]
MLstates_txt_nodes

tr2 = reorder(tr, "pruningwise")
plot(tr2, label.offset=0.15)
nodelabels(MLstates_txt_nodes)
tiplabels(MLstates_txt_tips)
title("Best states based on marginal joint splitprobs above node\nwith LAGRANGE 2-parameter model")

# Pie charts
MLstates_probs = get_ML_probs(relprobs_matrix=relprobs_matrix)
MLstates_probs = unlist(MLstates_probs)
MLstates_probs_other = 1-MLstates_probs
probs = cbind(MLstates_probs, MLstates_probs_other)

plot(tr2, label.offset=0.15)
#tiplabels(pie=relprobs_matrix[1:ntips,], piecol=piecols)
#nodelabels(pie=relprobs_matrix[(1+ntips):(length(MLstates_txt)),], piecol=piecols)
tiplabels(pie=probs[1:ntips,])
nodelabels(pie=probs[(1+ntips):(length(MLstates_txt)),])




# Compare probs (C++ Lagrange states, and BEARS 2-param marginal states)
plot(MLstates_LGcpp_justMLprob, MLstates_bears_2param, xlim=c(0,1), ylim=c(0,1))



# Compare probs:
#MLstates_LGcpp_justMLprob = MLstates_LGcpp$relprob[order(MLstates_LGcpp$Rnodes)]
nums = length(tr$tip.label)+1 : tr$Nnode
MLstates_bears_2param_states_dependent_on_joint_splits = probs[nums,1]
plot(MLstates_LGcpp_justMLprob, MLstates_bears_2param_states_dependent_on_joint_splits, xlim=c(0,1), ylim=c(0,1))















###########################################################
# Plot splitprob-based states (downpass likelihoods)
###########################################################
#results_object = bears_output
relprobs_matrix = joint_downpass_likelihoods_of_states_dependent_on_splits_above
MLstates_0based_indexes = get_ML_states(relprobs_matrix)

MLstates = list()
for (i in 1:length(MLstates_0based_indexes))
	{
	MLstates[[i]] = states[[ MLstates_0based_indexes[[i]]+1 ]]
	}

MLstates

MLstates_txt = unlist(lapply(X=MLstates, FUN=paste, sep="", collapse=""))
MLstates_txt

ntips = length(tr$tip.label)
MLstates_txt_tips = MLstates_txt[1:ntips]
MLstates_txt_nodes = MLstates_txt[(1+ntips):(length(MLstates_txt))]
MLstates_txt_nodes

tr2 = reorder(tr, "pruningwise")
plot(tr2, label.offset=0.15)
nodelabels(MLstates_txt_nodes)
tiplabels(MLstates_txt_tips)
title("Best states based on marginal joint splitprobs above node\nwith LAGRANGE 2-parameter model")

# Pie charts
MLstates_probs = get_ML_probs(relprobs_matrix=relprobs_matrix)
MLstates_probs = unlist(MLstates_probs)
MLstates_probs_other = 1-MLstates_probs
probs = cbind(MLstates_probs, MLstates_probs_other)

plot(tr2, label.offset=0.15)
#tiplabels(pie=relprobs_matrix[1:ntips,], piecol=piecols)
#nodelabels(pie=relprobs_matrix[(1+ntips):(length(MLstates_txt)),], piecol=piecols)
tiplabels(pie=probs[1:ntips,])
nodelabels(pie=probs[(1+ntips):(length(MLstates_txt)),])




# Compare probs (C++ Lagrange states, and BEARS 2-param marginal states)
plot(MLstates_LGcpp_justMLprob, MLstates_bears_2param, xlim=c(0,1), ylim=c(0,1))



# Compare probs:
#MLstates_LGcpp_justMLprob = MLstates_LGcpp$relprob[order(MLstates_LGcpp$Rnodes)]
nums = length(tr$tip.label)+1 : tr$Nnode
MLstates_bears_2param_states_dependent_on_joint_splits = probs[nums,1]
plot(MLstates_LGcpp_justMLprob, MLstates_bears_2param_states_dependent_on_joint_splits, xlim=c(0,1), ylim=c(0,1))



# None of these really line up.  Let's try optimization, conditional on a fixed state at a node.

# Fix the ancestral node (node 20) to state 5 (K_O), to see if we can get what
# LAGRANGE C++ gives us for this state:
# 18	K_O	0.516739	35.2022
# 18	K	0.278394	35.8207
# 18	K_O_M	0.0730309	37.1588
source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_models_v1.R', chdir = TRUE)
source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/calc_loglike_sp_v01.R', chdir = TRUE)

# Works!
fixnode = 20
fixlikes = rep(0, numstates)
fixlikes[5] = 1
optim_result = bears_2param_standard_fast_fixnode(trfn=trfn, geogfn=geogfn, fixnode=fixnode, fixlikes=fixlikes)
optim_result

bears_2param_standard_fast_result = bears_2param_standard_fast(trfn=trfn, geogfn=geogfn)
bears_2param_standard_fast_result = bears_2param_standard_fast_result$optim_result

# Fix all states at the nodes
# assumes numstates is WITHOUT null range

#######################################################
# This is quite slow, requires numnodes * numstates ML runs
#######################################################
runthis = FALSE

# Marginal probs of all states is this many ML runs
# 18 internal nodes
# 15 states
nn = 18
ns = 15
nn * ns

nums = length(tr$tip.label)+1 : tr$Nnode
numstates = numstates_from_numareas(numareas=length(areas), include_null_range=FALSE)
if (runthis)
	{
	optim_results = NULL
	for (nodenum in nums)
		{
		for (i in 1:numstates)
		#for (i in 1:2)
			{
			stateindex_0based_with_NULL_range = i
			statenum = stateindex_0based_with_NULL_range
	
			tmpstr = paste("\n\nStarting ML with node #", nodenum, " fixed to state ", statenum, "...\n\n", sep="")
			cat(tmpstr)
	
			
			fixnode = nodenum
			fixlikes = rep(0, numstates)
			fixlikes[i] = 1
			optim_result = bears_2param_standard_fast_fixnode(trfn=trfn, geogfn=geogfn, fixnode=fixnode, fixlikes=fixlikes)
			optim_result2 = c(fixnode, statenum, unlist(optim_result))
			
			optim_results = rbind(optim_results, optim_result2)
			}
		}
	
	optim_results = adf2(optim_results)
	names(optim_results) = c("fixnode", "statenum", "d", "e", "LnL", "method", "fns", "grs", "conv", "KKT1", "KKT2", "xtimes")
	optim_results = dfnums_to_numeric(optim_results)
	optim_results
	
	save(optim_results, file="/Dropbox/_njm/__packages/BioGeoBEARS_setup/_examples/Psychotria_2param_MLstates_local.Rdata")
	# end runthis
	} else {
	# Load optim_results
	load(file="/Dropbox/_njm/__packages/BioGeoBEARS_setup/_examples/Psychotria_2param_MLstates_local.Rdata")
	}


global_ML_LnL = bears_2param_standard_fast_result$fvalues[[1]]
relprobs_under_local_ML = exp(optim_results$LnL) / exp(global_ML_LnL)
relprob = relprobs_under_local_ML

colnums = (1:ncol(optim_results))
LnLcol = colnums[names(optim_results) == "LnL"]
optim_results2 = cbind(optim_results[, 1:LnLcol], relprob, optim_results[, (LnLcol+1):ncol(optim_results)])

# Lagrange node 18
subtable = optim_results2[optim_results2$fixnode==20, ]
subtable[rev(order(subtable$relprob)), ]

# Lagrange node 1
subtable = optim_results2[optim_results2$fixnode==27, ]
subtable[rev(order(subtable$relprob)), ]

# Lagrange node 11
subtable = optim_results2[optim_results2$fixnode==32, ]
subtable[rev(order(subtable$relprob)), ]

get_lagrange_nodenums(tr)

# Compare to C++ LAGRANGE:
# 1	M	0.98823	34.5538
# 2	O	0.982444	34.5597
# 3	O_M	0.980038	34.5621
# 4	M	0.999563	34.5424
# 5	O_M	0.955127	34.5879
# 6	O_M	0.549586	35.1406
# 6	O	0.416161	35.4187
# 7	K	0.990151	34.5519
# 8	K_O	0.480417	35.2751
# 8	K_O_M	0.425713	35.396
# 8	K_M	0.0699342	37.2022
# 9	M_H	0.538402	35.1611
# 9	K_M_H	0.0811531	37.0534
# 10	O_M_H	0.174998	36.2849
# 10	O_M	0.142277	36.4919
# 10	O_H	0.137926	36.523
# 10	O	0.095075	36.8951
# 10	K	0.079093	37.0791
# 10	K_O	0.0662898	37.2557
# 10	M	0.0545041	37.4514
# 10	H	0.052816	37.4829
# 10	K_M	0.0413575	37.7275
# 10	K_H	0.0409942	37.7363
# 10	K_O_M	0.0321924	37.978
# 10	K_O_H	0.0316803	37.994
# 11	K	0.221002	36.0516
# 11	K_O	0.189177	36.207
# 11	K_M	0.112554	36.7263
# 11	K_O_M	0.111695	36.7339
# 11	K_O_M_H	0.110042	36.7489
# 11	K_H	0.109137	36.7571
# 11	K_O_H	0.107115	36.7758
# 11	K_M_H	0.033549	37.9367
# 12	K	0.634435	34.997
# 12	K_O	0.104353	36.8019
# 13	K	0.671399	34.9404
# 13	K_O	0.134266	36.5499
# 14	K	0.992332	34.5497
# 15	K	0.800718	34.7642
# 16	K	0.968912	34.5736
# 17	K_O	0.856548	34.6968
# 18	K_O	0.516739	35.2022
# 18	K	0.278394	35.8207
# 18	K_O_M	0.0730309	37.1588

global_LnL = -34.542
LnL_1 = -35.2022
LnL_2 = -35.8207
LnL_3 = -37.1588

exp(LnL_1) / exp(global_LnL)
exp(LnL_2) / exp(global_LnL)
exp(LnL_3) / exp(global_LnL)


#######################################################
# Basically, the above works!  Sometimes things close in
# probability have slightly differing orders, but this is
# flat likelihood surface...
#######################################################








#######################################################
# Do Local ML state estimates using 3-param model (Sloooooow!)
#######################################################
source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_models_v1.R', chdir = TRUE)
source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/calc_loglike_sp_v01.R', chdir = TRUE)

# Works!
fixnode = 20
fixlikes = rep(0, numstates)
fixlikes[5] = 1
optim_result = bears_3param_standard_fast_fixnode(trfn=trfn, geogfn=geogfn, fixnode=fixnode, fixlikes=fixlikes)
optim_result


bears_3param_standard_fast_result = bears_3param_standard_fast(trfn=trfn, geogfn=geogfn)
bears_3param_standard_fast_result = bears_3param_standard_fast_result$optim_result

#######################################################
# This is quite slow, requires numnodes * numstates ML runs
#######################################################
runthis = TRUE

# Marginal probs of all states is this many ML runs
# 18 internal nodes
# 15 states
nn = 18
ns = 15
nn * ns

nums = length(tr$tip.label)+1 : tr$Nnode
numstates = numstates_from_numareas(numareas=length(areas), include_null_range=FALSE)
if (runthis)
	{
	optim_results = NULL
	for (nodenum in nums)
		{
		for (i in 1:numstates)
		#for (i in 1:2)
			{
			stateindex_0based_with_NULL_range = i
			statenum = stateindex_0based_with_NULL_range
	
			tmpstr = paste("\n\nStarting ML with node #", nodenum, " fixed to state ", statenum, "...\n\n", sep="")
			cat(tmpstr)
	
			
			fixnode = nodenum
			fixlikes = rep(0, numstates)
			fixlikes[i] = 1
			optim_result = bears_3param_standard_fast_fixnode(trfn=trfn, geogfn=geogfn, fixnode=fixnode, fixlikes=fixlikes)
			optim_result2 = c(fixnode, statenum, unlist(optim_result))
			
			optim_results = rbind(optim_results, optim_result2)
			}
		}
	
	optim_results = adf2(optim_results)
	names(optim_results) = c("fixnode", "statenum", "d", "e", "j", "LnL", "method", "fns", "grs", "conv", "KKT1", "KKT2", "xtimes")
	optim_results = dfnums_to_numeric(optim_results)
	optim_results
	
	save(optim_results, file="/Dropbox/_njm/__packages/BioGeoBEARS_setup/_examples/Psychotria_3param_MLstates_local.Rdata")
	# end runthis
	} else {
	# Load optim_results
	load(file="/Dropbox/_njm/__packages/BioGeoBEARS_setup/_examples/Psychotria_3param_MLstates_local.Rdata")
	}




# Get the statenames
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
areas = getareas_from_tipranges_object(tipranges)
states = areas_list_to_states_list_old(areas=areas)
statenames = sapply(FUN=paste, states, sep="", collapse="")
statenames = statenames[-1]
statenames

global_LnL = bears_3param_standard_fast_result$fvalues[[1]]
MLlocal_ancstates_to_relprob_states <- function(tr, optim_results, global_LnL, addnull=TRUE)
	{
	#get_APE_nodenums(tr)
	cat("\nGlobal LnL is: ", global_LnL, "\n", sep="")
	
	relprobs = exp(optim_results$LnL) / exp(global_LnL)
	
	local_ML_relprobs_matrix = matrix(data=relprobs, nrow=tr$Nnode, ncol=length(statenames), byrow=TRUE)
	local_ML_relprobs_matrix
	
	# rowSums sometimes slightly above 1, normalize
	local_ML_relprobs_matrix = local_ML_relprobs_matrix / rowSums(local_ML_relprobs_matrix)
	rowSums(local_ML_relprobs_matrix)
	
	if (addnull == TRUE)
		{
		zeros_col = matrix(data=0, ncol=1, nrow=tr$Nnode)
		local_ML_relprobs_matrix = cbind(zeros_col, local_ML_relprobs_matrix)
		}
	
	return(local_ML_relprobs_matrix)
	}


local_ML_relprobs_matrix = MLlocal_ancstates_to_relprob_states(tr, optim_results, global_LnL, addnull=TRUE)
tip_condlikes_of_data_on_each_state = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges, phy=tr)
local_ML_relprobs_matrix = rbind(tip_condlikes_of_data_on_each_state, local_ML_relprobs_matrix)

relprobs_matrix = local_ML_relprobs_matrix
MLstates_0based_indexes = get_ML_states(relprobs_matrix)



MLstates = list()
for (i in 1:length(MLstates_0based_indexes))
	{
	MLstates[[i]] = states[[ MLstates_0based_indexes[[i]]+1 ]]
	}

MLstates

MLstates_txt = unlist(lapply(X=MLstates, FUN=paste, sep="", collapse=""))
MLstates_txt

ntips = length(tr$tip.label)
MLstates_txt_tips = MLstates_txt[1:ntips]
MLstates_txt_nodes = MLstates_txt[(1+ntips):(length(MLstates_txt))]
MLstates_txt_nodes

tr2 = reorder(tr, "pruningwise")
plot(tr2, label.offset=0.15)
nodelabels(MLstates_txt_nodes)
tiplabels(MLstates_txt_tips)
title("Best states based on marginal joint stateprobs at node\nwith BEARS 3-parameter model")

# Pie charts
MLstates_probs = get_ML_probs(relprobs_matrix=relprobs_matrix)
MLstates_probs = unlist(MLstates_probs)
MLstates_probs_other = 1-MLstates_probs
probs = cbind(MLstates_probs, MLstates_probs_other)

plot(tr2, label.offset=0.15)
#tiplabels(pie=relprobs_matrix[1:ntips,], piecol=piecols)
#nodelabels(pie=relprobs_matrix[(1+ntips):(length(MLstates_txt)),], piecol=piecols)
tiplabels(pie=probs[1:ntips,])
nodelabels(pie=probs[(1+ntips):(length(MLstates_txt)),])




tr2 = reorder(tr, "pruningwise")
plot(tr2, label.offset=0.15)
nodelabels(MLstates_txt_nodes)
tiplabels(MLstates_txt_tips)
title("Best states based on marginal joint stateprobs at node\nwith BEARS 3-parameter model")

# Pie charts
MLstates_probs = get_ML_probs(relprobs_matrix=relprobs_matrix)
MLstates_probs = unlist(MLstates_probs)
MLstates_probs_other = 1-MLstates_probs
probs = cbind(MLstates_probs, MLstates_probs_other)

plot(tr2, label.offset=0.15)
tiplabels(pie=relprobs_matrix[1:ntips,], piecol=piecols)
nodelabels(pie=relprobs_matrix[(1+ntips):(length(MLstates_txt)),], piecol=piecols)















#######################################################
# Correct LAGRANGE node ordering!
#######################################################
downpass_node_matrix = get_lagrange_nodenums(tr)
downpass_node_matrix = downpass_node_matrix[order(downpass_node_matrix[,1]), ]

plot(tr)
nodelabels(node=20:37, downpass_node_matrix[,1])
tiplabels(1:19)

# THIS WORKS
plot(tr)
nodelabels(node=20:37, downpass_node_matrix[,2])
tiplabels(1:19)





plot(tr)
nodelabels(text=MLsplits_LGcpp$LGnodes)
tiplabels()
title("LAGRANGE (C++) numbers")



plot(tr)
nodelabels(node=20:37, 38-postorder_table$nodenums)
tiplabels(1:19)


plot(tr)
nodelabels(node=20:37, 19-postorder_table$internal_nodenums)
tiplabels(1:19)

plot(tr)
nodelabels(node=20:37, 19-postorder_table$postorder)
tiplabels(1:19)

plot(tr)
nodelabels(node=20:37, 38-postorder_table$DIVA_postorder)
tiplabels(1:19)

tr4 = as(tr, "phylo4")
tr4 = reorder(tr4, "preorder")
plot(as(tr4, "phylo"))
attr(as(tr4, "phylo"),"order")
nodelabels()
tiplabels(1:19)

tr4 = as(tr, "phylo4")
tr4 = reorder(tr4, "postorder")
plot(as(tr4, "phylo"))
attr(as(tr4, "phylo"),"order")
nodelabels()
tiplabels(1:19)

plot(tr)
nodelabels(node=20:37, postorder_table$nodenums)
tiplabels(1:19)
attr(tr,"order")



sourcedir = '/Dropbox/_njm/'
source3 = '_R_tree_functions_v1.R'
source(paste(sourcedir, source3, sep=""))

tr = read.tree(trfn)
attr(tr,"order")

node_numbers = get_lagrange_nodenums(tr)
node_numbers

tr = reorder(tr, "pruningwise")

plot(tr)
nodelabels(text=node_numbers[,1], node=20:37)
tiplabels(1:19)

plot(tr)
nodelabels(text=node_numbers[,2], node=20:37)
tiplabels(1:19)

plot(tr)
nodelabels(text=node_numbers[order(node_numbers[,1]),1], node=20:37)
tiplabels(1:19)

plot(tr)
nodelabels(text=node_numbers[order(node_numbers[,2]),1], node=20:37)
tiplabels(1:19)








sourcedir = '/Dropbox/_njm/'
source3 = '_R_tree_functions_v1.R'
source(paste(sourcedir, source3, sep=""))

tr = read.tree(trfn)
attr(tr,"order")

res = get_lagrange_nodenums2(tr)
node_numbers = res[[1]]
node_numbers

tr2 = res[[2]]
#tr = reorder(tr, "pruningwise")

plot(tr2)
nodelabels(text=node_numbers[,1], node=20:37)
tiplabels(1:19)

plot(tr2)
nodelabels(text=node_numbers[,2], node=20:37)
tiplabels(1:19)

plot(tr2)
nodelabels(text=node_numbers[order(node_numbers[,1]),1], node=20:37)
tiplabels(1:19)

plot(tr2)
nodelabels(text=node_numbers[order(node_numbers[,2]),1], node=20:37)
tiplabels(1:19)




















#tmporder = order(postorder_table$internal_nodenums)
#tmporder = order(postorder_table$nodenums)
tmporder = postorder_table$internal_nodenums
tmporder = order(postorder_table$postorder)
tmporder = rev(match(x=MLsplits_LGcpp$nodenum_LGcpp, postorder_table$internal_nodenums))
postorder_table$postorder[order(postorder_table$internal_nodenums)]
#tmporder = order(MLsplits_LGcpp$nodenum_ORD1)
#tmporder = order(rev(MLsplits_LGcpp$nodenum_ORD1))
tmporder = order(rev(postorder_table$internal_nodenums))
tmporder = order(postorder_table$nodenums)

par(mfrow=c(2,1))
plot(tr)
nodelabels(node=20:37, MLsplits_LGcpp$splits[tmporder])
tiplabels(tiprange_names)
title("LAGRANGE (C++) ML splits")

plot(tr)
pievals = as.matrix(MLsplits_LGcpp[,c("relprob","relprob2")])
nodelabels(node=20:37, pie=pievals[tmporder,], piecol=c("blue", "white"))
tiplabels(tiprange_names)
title("LAGRANGE (C++) split probs")




