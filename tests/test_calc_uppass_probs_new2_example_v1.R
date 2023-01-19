

#######################################################
# EXAMPLE TO CHECK calc_uppass_probs_new2
#######################################################
calc_uppass_probs_new2_example <- function()
	{
# dput(resDEC)

# dput(resDEC)

resDEC = structure(list(computed_likelihoods_at_each_node = c(1, 1, 1, 
1, 0.36502788215878, 0.227783575654047, 0.136167332530107), relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS = structure(c(0, 
0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 6.42693382782259e-14, 7.46190547654875e-14, 
3.13642570191625e-13, 1, 0, 1, 1, 0.757828232249181, 0.344749222130095, 
3.13642570191625e-13, 0, 0, 0, 0, 0.242171767750754, 0.65525077786983, 
0.999999999999373), dim = c(7L, 4L)), condlikes_of_each_state = structure(c(0, 
0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2.3460100439447e-14, 1.69969951064079e-14, 
4.27078721508805e-14, 1, 0, 1, 1, 0.276628434658051, 0.0785282105207443, 
4.27078721508805e-14, 0, 0, 0, 0, 0.0883994475007057, 0.149255365133286, 
0.136167332530022), dim = c(7L, 4L)), relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS = structure(c(0, 
0, 0, 0, NA, 0, 0, 5.22737626427221e-14, 0.999999999998895, 2.16444459101218e-13, 
5.04411011774348e-13, NA, 0.0576312779227834, 0.0806194459381856, 
0.999999999998895, 5.22737626427221e-14, 0.999999999997567, 0.999999999995991, 
NA, 0.342775485350427, 0.0806194459381856, 1.05227375864476e-12, 
1.05227375864476e-12, 2.21644445110031e-12, 3.5044109997644e-12, 
NA, 0.599593236726789, 0.838761108123629), dim = c(7L, 4L)), 
    relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS = structure(c(0, 
    0, 0, 0, NA, 0, 0, 0.477617887308721, 0.0266167777965505, 
    0.308323346796742, 0.292600516988572, NA, 0.125000000000644, 
    0.046615618818038, 0.261191056346013, 0.946766444406928, 
    0.649573122214116, 0.616448394398207, NA, 0.749999999999466, 
    0.90676876236405, 0.261191056345266, 0.0266167777965214, 
    0.0421035309891421, 0.090951088613221, NA, 0.12499999999989, 
    0.0466156188179121), dim = c(7L, 4L)), relative_probs_of_each_state_at_branch_top_AT_node_UPPASS = structure(c(7.02704868858637e-13, 
    9.25815984473407e-13, 1.73460095194895e-12, 2.35254814927231e-12, 
    0.25, 8.32240548158005e-13, 9.06794445714384e-13, 0.431710549603533, 
    0.0240584452060588, 0.251901392427722, 0.216078384422621, 
    0.25, 0.112985338561828, 0.0421350517952136, 0.236086079443338, 
    0.855765818076806, 0.530703807120855, 0.455232186573681, 
    0.25, 0.67791203136619, 0.819612604898638, 0.332203370953129, 
    0.120175736717135, 0.217394800451423, 0.328689429003699, 
    0.25, 0.209102630071982, 0.138252343306148), dim = c(7L, 
    4L)), ML_marginal_prob_each_state_at_branch_bottom_below_node = structure(c(0, 
    0, 0, 0, NA, 0, 0, 9.55885872371469e-14, 0.999999999997088, 
    1.02736516865682e-13, 2.3942137600093e-13, NA, 0.0212357703981079, 
    0.0324086154040224, 0.999999999998852, 1.85939277741373e-12, 
    0.999999999999754, 0.999999999999244, NA, 0.757828224601766, 
    0.630413600097194, 1.05227375864171e-12, 1.05227375864171e-12, 
    1.43663791560109e-13, 5.17042461743951e-13, NA, 0.220936005000126, 
    0.337177784498784), dim = c(7L, 4L)), ML_marginal_prob_each_state_at_branch_top_AT_node = structure(c(0, 
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 6.42693382782259e-14, 2.27415872607374e-14, 
    9.55885855108166e-14, 1, 0, 1, 1, 0.757828232249181, 0.630413602214152, 
    1.85939274383416e-12, 0, 0, 0, 0, 0.242171767750754, 0.369586397785825, 
    0.999999999998045), dim = c(7L, 4L)), relative_probs_of_each_state_at_bottom_of_root_branch = c(0, 
    6.42693382782259e-14, 0.757828232249181, 0.242171767750754
    ), total_loglikelihood = -4.48101163256014, inputs = list(
        geogfn = "geog.data", trfn = "tree.newick", abbr = "default", 
        description = "defaults", BioGeoBEARS_model_object = new("BioGeoBEARS_model", 
            params_table = structure(list(type = c("free", "free", 
            "fixed", "fixed", "fixed", "fixed", "fixed", "fixed", 
            "fixed", "3-j", "ysv*2/3", "ysv*1/3", "ysv*1/3", 
            "ysv*1/3", "fixed", "mx01", "mx01", "mx01", "mx01", 
            "fixed", "fixed", "fixed", "fixed"), init = c(0.101055674351941, 
            1e-12, 0, 1, 0, 0, 1, 0, 0, 2.99999, 1.99999, 1, 
            1, 1, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 0.5, 0.1, 
            1, 0), min = c(1e-12, 1e-12, 1e-12, 1e-12, -2.5, 
            -10, -10, -10, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 
            1e-05, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 
            0.005, 0.005, 0.005), max = c(4.999999999999, 4.999999999999, 
            4.999999999999, 0.999999999999, 2.5, 10, 10, 10, 
            2.99999, 3, 2, 1, 1, 1, 0.9999, 0.9999, 0.9999, 0.9999, 
            0.9999, 0.9999, 0.995, 0.995, 0.995), est = c(0.101055674351941, 
            1e-12, 0, 1, 0, 0, 1, 0, 0, 3, 2, 1, 1, 1, 1e-04, 
            1e-04, 1e-04, 1e-04, 1e-04, 0.5, 0.1, 1, 0), note = c("works", 
            "works", "works", "non-stratified only", "works", 
            "works", "works", "works", "works", "works", "works", 
            "works", "works", "works", "works", "works", "works", 
            "works", "works", "no", "yes", "yes", "yes"), desc = c("anagenesis: rate of 'dispersal' (range expansion)", 
            "anagenesis: rate of 'extinction' (range contraction)", 
            "anagenesis: rate of range-switching (i.e. for a standard char.)", 
            "anagenesis: exponent on branch lengths", "exponent on distance (modifies d, j, a)", 
            "exponent on environmental distance (modifies d, j, a)", 
            "exponent on manual dispersal multipliers (modifies d, j, a)", 
            "anagenesis: exponent on extinction risk with area (modifies e)", 
            "cladogenesis: relative per-event weight of jump dispersal", 
            "cladogenesis: y+s+v", "cladogenesis: y+s", "cladogenesis: relative per-event weight of sympatry (range-copying)", 
            "cladogenesis: relative per-event weight of subset speciation", 
            "cladogenesis: relative per-event weight of vicariant speciation", 
            "cladogenesis: controls range size of smaller daughter", 
            "cladogenesis: controls range size of smaller daughter", 
            "cladogenesis: controls range size of smaller daughter", 
            "cladogenesis: controls range size of smaller daughter", 
            "cladogenesis: controls range size of smaller daughter", 
            "root: controls range size probabilities of root", 
            "mean frequency of truly sampling OTU of interest", 
            "detection probability per true sample of OTU of interest", 
            "false detection of OTU probability per true taphonomic control sample"
            )), row.names = c("d", "e", "a", "b", "x", "n", "w", 
            "u", "j", "ysv", "ys", "y", "s", "v", "mx01", "mx01j", 
            "mx01y", "mx01s", "mx01v", "mx01r", "mf", "dp", "fdp"
            ), class = "data.frame")), timesfn = NA, distsfn = NA, 
        dispersal_multipliers_fn = NA, area_of_areas_fn = NA, 
        areas_allowed_fn = NA, areas_adjacency_fn = NA, detects_fn = NA, 
        controls_fn = NA, max_range_size = 2, force_sparse = FALSE, 
        use_detection_model = FALSE, print_optim = TRUE, printlevel = 0, 
        on_NaN_error = -1e+50, wd = "/GitHub/PhyBEARS.jl/test/apes_SSE", 
        num_cores_to_use = 1, cluster_already_open = FALSE, use_optimx = TRUE, 
        rescale_params = FALSE, return_condlikes_table = TRUE, 
        calc_TTL_loglike_from_condlikes_table = TRUE, calc_ancprobs = TRUE, 
        speedup = TRUE, include_null_range = TRUE, useAmbiguities = FALSE, 
        min_branchlength = 1e-06, allow_null_tips = FALSE, all_geog_states_list_usually_inferred_from_areas_maxareas = list(
            NA, 0L, 1L, 0:1)), outputs = new("BioGeoBEARS_model", 
        params_table = structure(list(type = c("free", "free", 
        "fixed", "fixed", "fixed", "fixed", "fixed", "fixed", 
        "fixed", "3-j", "ysv*2/3", "ysv*1/3", "ysv*1/3", "ysv*1/3", 
        "fixed", "mx01", "mx01", "mx01", "mx01", "fixed", "fixed", 
        "fixed", "fixed"), init = c(0.101055674351941, 1e-12, 
        0, 1, 0, 0, 1, 0, 0, 2.99999, 1.99999, 1, 1, 1, 1e-04, 
        1e-04, 1e-04, 1e-04, 1e-04, 0.5, 0.1, 1, 0), min = c(1e-12, 
        1e-12, 1e-12, 1e-12, -2.5, -10, -10, -10, 1e-05, 1e-05, 
        1e-05, 1e-05, 1e-05, 1e-05, 1e-04, 1e-04, 1e-04, 1e-04, 
        1e-04, 1e-04, 0.005, 0.005, 0.005), max = c(4.999999999999, 
        4.999999999999, 4.999999999999, 0.999999999999, 2.5, 
        10, 10, 10, 2.99999, 3, 2, 1, 1, 1, 0.9999, 0.9999, 0.9999, 
        0.9999, 0.9999, 0.9999, 0.995, 0.995, 0.995), est = c(0.101055674351941, 
        1e-12, 0, 1, 0, 0, 1, 0, 0, 3, 2, 1, 1, 1, 1e-04, 1e-04, 
        1e-04, 1e-04, 1e-04, 0.5, 0.1, 1, 0), note = c("works", 
        "works", "works", "non-stratified only", "works", "works", 
        "works", "works", "works", "works", "works", "works", 
        "works", "works", "works", "works", "works", "works", 
        "works", "no", "yes", "yes", "yes"), desc = c("anagenesis: rate of 'dispersal' (range expansion)", 
        "anagenesis: rate of 'extinction' (range contraction)", 
        "anagenesis: rate of range-switching (i.e. for a standard char.)", 
        "anagenesis: exponent on branch lengths", "exponent on distance (modifies d, j, a)", 
        "exponent on environmental distance (modifies d, j, a)", 
        "exponent on manual dispersal multipliers (modifies d, j, a)", 
        "anagenesis: exponent on extinction risk with area (modifies e)", 
        "cladogenesis: relative per-event weight of jump dispersal", 
        "cladogenesis: y+s+v", "cladogenesis: y+s", "cladogenesis: relative per-event weight of sympatry (range-copying)", 
        "cladogenesis: relative per-event weight of subset speciation", 
        "cladogenesis: relative per-event weight of vicariant speciation", 
        "cladogenesis: controls range size of smaller daughter", 
        "cladogenesis: controls range size of smaller daughter", 
        "cladogenesis: controls range size of smaller daughter", 
        "cladogenesis: controls range size of smaller daughter", 
        "cladogenesis: controls range size of smaller daughter", 
        "root: controls range size probabilities of root", "mean frequency of truly sampling OTU of interest", 
        "detection probability per true sample of OTU of interest", 
        "false detection of OTU probability per true taphonomic control sample"
        )), row.names = c("d", "e", "a", "b", "x", "n", "w", 
        "u", "j", "ysv", "ys", "y", "s", "v", "mx01", "mx01j", 
        "mx01y", "mx01s", "mx01v", "mx01r", "mf", "dp", "fdp"
        ), class = "data.frame")), optim_result = structure(list(
        p1 = 0.101055674351941, p2 = 1e-12, value = -4.48101163256014, 
        fevals = 39, gevals = NA_real_, niter = NA_real_, convcode = 0, 
        kkt1 = FALSE, kkt2 = TRUE, xtime = 0.354), row.names = "bobyqa", class = c("optimx", 
    "data.frame"), details = structure(list("bobyqa", c(0.0313402484276759, 
    0), structure(c(-65.7496196774338, 0, 0, 0), dim = c(2L, 
    2L)), c(-65.7496196774338, 0), "none"), dim = c(1L, 5L), dimnames = list(
        "bobyqa", c("method", "ngatend", "nhatend", "hev", "message"
        ))), maximize = TRUE, npar = 2L, follow.on = FALSE)), class = "calc_loglike_sp_results")

trfn = "/GitHub/PhyBEARS.jl/test/apes_SSE/tree.newick"
moref(trfn)

tree_string = "(((chimp:1,human:1):1,gorilla:2):1,orang:3);"
tr = read.tree(file="", text=tree_string)
trtable = prt(tr, printflag=FALSE)

numstates = 4
tmpres = get_Qmat_COOmat_from_res(resDEC, numstates=numstates, include_null_range=TRUE, max_range_size=NULL, timeperiod_i=1)

probs_ancstate = rep(0.25, numstates)
COO_weights_columnar = tmpres$COO_weights_columnar
include_null_range = TRUE
left_branch_downpass_likes = rep(1, numstates)
right_branch_downpass_likes = rep(1, numstates)
Rsp_rowsums = tmpres$Rsp_rowsums


clado_table = cbind(tmpres$COO_weights_columnar[[1]]+1+include_null_range, tmpres$COO_weights_columnar[[2]]+1+include_null_range, tmpres$COO_weights_columnar[[3]]+1+include_null_range, tmpres$COO_weights_columnar[[4]])

probs = rep(0, nrow(clado_table))
for (i in 1:length(probs))
	probs[i] = clado_table[i,4] / Rsp_rowsums[clado_table[i,1]-include_null_range]
end
clado_table = cbind(clado_table, probs)
clado_table_df = adf2(clado_table)
names(clado_table_df) = c("i", "j", "k", "wt", "prob")
clado_table_df

# equal ancstate probs
calc_uppass_probs_new2(probs_ancstate, COO_weights_columnar, numstates, include_null_range = include_null_range, 
    left_branch_downpass_likes = NULL, right_branch_downpass_likes = NULL, 
    Rsp_rowsums = NULL)

# $condprob_each_split_scenario_df2
#            [,1]
# [1,] 0.33333333
# [2,] 0.33333333
# [3,] 0.05555556
# [4,] 0.05555556
# [5,] 0.05555556
# [6,] 0.05555556
# [7,] 0.05555556
# [8,] 0.05555556
# 
# $relprobs_just_after_speciation_UPPASS_Left
# [1] 0.0000000 0.4444444 0.4444444 0.1111111
# 
# $relprobs_just_after_speciation_UPPASS_Right
# [1] 0.0000000 0.4444444 0.4444444 0.1111111



# Ancestral node states

# Ancestral root node number: 5
# Left branch node number: 6  (internal)
# Right branch node number: 4 (a tip, orangutan)

numstates = ncol(resDEC$ML_marginal_prob_each_state_at_branch_top_AT_node)
rootnode = 5
#probs_ancstate = resDEC$ML_marginal_prob_each_state_at_branch_top_AT_node[5,]
probs_ancstate = rep(1/numstates, numstates)
calc_uppass_probs_new2(probs_ancstate, COO_weights_columnar, numstates, include_null_range = include_null_range, 
    left_branch_downpass_likes = NULL, right_branch_downpass_likes = NULL, 
    Rsp_rowsums = NULL)

# $condprob_each_split_scenario_df2
#              [,1]
# [1,] 6.426934e-14
# [2,] 7.578282e-01
# [3,] 4.036196e-02
# [4,] 4.036196e-02
# [5,] 4.036196e-02
# [6,] 4.036196e-02
# [7,] 4.036196e-02
# [8,] 4.036196e-02
# 
# $relprobs_just_after_speciation_UPPASS_Left
# [1] 0.00000000 0.08072392 0.83855215 0.08072392
# 
# $relprobs_just_after_speciation_UPPASS_Right
# [1] 0.00000000 0.08072392 0.83855215 0.08072392

# Calculate uppass probabilities for RIGHT branch, the corner below node 4 (sister branch is node 6)
left_branch_downpass_likes = resDEC$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[6,]
right_branch_downpass_likes = rep(1.0, numstates)
tmpres = calc_uppass_probs_new2(probs_ancstate, COO_weights_columnar, numstates, include_null_range = include_null_range, 
    left_branch_downpass_likes=left_branch_downpass_likes, right_branch_downpass_likes=right_branch_downpass_likes, 
    Rsp_rowsums = NULL)

tmpres$relprobs_just_after_speciation_UPPASS_Right

calc = tmpres$relprobs_just_after_speciation_UPPASS_Right * resDEC$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[4,]
calc / sum(calc)

resDEC$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS[4,]
resDEC$ML_marginal_prob_each_state_at_branch_bottom_below_node[4,]

# > tmpres$relprobs_just_after_speciation_UPPASS_Right
# [1] 0.00000000 0.29260052 0.61644839 0.09095109
# > 
# > calc = tmpres$relprobs_just_after_speciation_UPPASS_Right * resDEC$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[4,]
# > calc / sum(calc)
# [1] 0.000000e+00 2.394214e-13 1.000000e+00 5.170425e-13
# > 
# > resDEC$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS[4,]
# [1] 0.00000000 0.29260052 0.61644839 0.09095109
# > resDEC$ML_marginal_prob_each_state_at_branch_bottom_below_node[4,]
# [1] 0.000000e+00 2.394214e-13 1.000000e+00 5.170425e-13



# Calculate uppass probabilities for LEFT branch, the corner below node 6 (sister branch is node 4)
left_branch_downpass_likes = rep(1.0, numstates)
right_branch_downpass_likes = resDEC$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[4,]
tmpres = calc_uppass_probs_new2(probs_ancstate, COO_weights_columnar, numstates, include_null_range = include_null_range, 
    left_branch_downpass_likes=left_branch_downpass_likes, right_branch_downpass_likes=right_branch_downpass_likes, 
    Rsp_rowsums = NULL)

tmpres$relprobs_just_after_speciation_UPPASS_Left

calc = tmpres$relprobs_just_after_speciation_UPPASS_Left * resDEC$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[6,]
calc / sum(calc)

resDEC$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS[6,]
resDEC$ML_marginal_prob_each_state_at_branch_bottom_below_node[6,]


# > tmpres$relprobs_just_after_speciation_UPPASS_Left
# [1] 0.000 0.125 0.750 0.125
# > 
# > calc = tmpres$relprobs_just_after_speciation_UPPASS_Left * resDEC$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[6,]
# > calc / sum(calc)
# [1] 0.00000000 0.02123577 0.75782822 0.22093601
# > 
# > resDEC$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS[6,]
# [1] 0.000 0.125 0.750 0.125
# > resDEC$ML_marginal_prob_each_state_at_branch_bottom_below_node[6,]
# [1] 0.00000000 0.02123577 0.75782822 0.22093601







# Check uppass

mats = get_Qmat_COOmat_from_res(resDEC, numstates=numstates, include_null_range=TRUE, max_range_size=NULL, timeperiod_i=1)
mats = tmpres$Qmat

u0 = c(0.0, 0.125, 0.75, 0.125)
u0 %*% expm()
	} # END calc_uppass_probs_new2_example <- function()
	
	
	