
#######################################################
# Source code
#######################################################
library(LaplacesDemon)
library(BioGeoBEARS)
source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_generics_v1.R', chdir = TRUE)
source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_classes_v1.R', chdir = TRUE)
source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_univ_model_v1.R', chdir = TRUE)
source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_readwrite_v1.R', chdir = TRUE)
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

wd = slashslash(paste(extdata_dir, "/examples/Psychotria_M0/LAGRANGE_python/", sep=""))
setwd(wd)


#######################################################
# Setup for reading Python/C++ LAGRANGE results
#######################################################
outfn = "Psychotria_5.2_demo.results.txt"
results_dir = wd
new_splits_fn = TRUE
new_states_fn = FALSE
filecount=0

summstats = parse_lagrange_python_output(outfn, results_dir, new_splits_fn = TRUE, new_states_fn = FALSE, filecount=0)
summstats




#######################################################
# Input files (tree and tipdata)
#######################################################
trfn = "../Psychotria_5.2.newick"

# Plot
geogfn = paste(addslash(extdata_dir), "Psychotria_geog.data", sep="")
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tiparea_names = c("K", "O", "M", "H")
names(tipranges@df) = tiparea_names
tipranges

tiprangesTF = tipranges@df == 1



# ORDER THE FUCKING TIPNAMES
tiprange_names = apply(X=tiprangesTF, 1, getname, tiparea_names)
matchnums = get_indices_where_list1_occurs_in_list2(list1=tr$tip.label, list2=names(tiprange_names))
tiprange_names = tiprange_names[matchnums]
names(tiprange_names) ==  tr$tip.label

# 
# 
# 
# sourcedir = '/Dropbox/_njm/'
# source3 = '_R_tree_functions_v1.R'
# source(paste(sourcedir, source3, sep=""))
# 





#######################################################
# OK now read the Python LAGRANGE version
#######################################################

splits_fn = "Psychotria_5.2_demo.results_splits00001.txt"
moref(splits_fn)

splits = LGpy_splits_fn_to_table(splits_fn)
splits
class(splits)


MLsplits_LGpy = LGpy_MLsplit_per_node(splits)
MLsplits_LGpy
tmpnames = names(MLsplits_LGpy)

# Order in LAGRANGE order
tr = read.tree(trfn)
downpass_node_matrix = get_lagrange_nodenums(tr, option=2)
downpass_node_matrix[,2] = 1:18
#downpass_node_matrix = downpass_node_matrix[order(downpass_node_matrix[,1]), ]

MLsplits_LGpy = MLsplits_LGpy[order(MLsplits_LGpy$nodenum_LGpy),]
MLsplits_LGpy

MLsplits_LGpy = cbind(MLsplits_LGpy, downpass_node_matrix)
names(MLsplits_LGpy) = c(tmpnames, "Rnodes", "LGnodes")
MLsplits_LGpy

MLsplits_LGpy = MLsplits_LGpy[order(MLsplits_LGpy$Rnodes), ]
MLsplits_LGpy
#cls.df(MLsplits_LGpy)


par(mfrow=c(2,1))
plot(tr)
nodelabels(text=MLsplits_LGpy$splits, node=20:37)
tiplabels(tiprange_names)
title("LAGRANGE (python) ML splits")

plot(tr)
pievals = as.matrix(MLsplits_LGpy[,c("relprob","relprob2")])
nodelabels(pie=pievals, piecol=c("blue", "white"))
tiplabels(tiprange_names)
title("LAGRANGE (python) split probs")



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
tmpnames = names(MLsplits_LGcpp)

# Order in LAGRANGE order
tr = read.tree(trfn)
downpass_node_matrix = get_lagrange_nodenums(tr, option=1)

MLsplits_LGcpp = MLsplits_LGcpp[order(MLsplits_LGcpp$nodenum_LGcpp),]
MLsplits_LGcpp

MLsplits_LGcpp = cbind(MLsplits_LGcpp, downpass_node_matrix)
MLsplits_LGcpp
names(MLsplits_LGcpp) = c(tmpnames, "Rnodes", "LGnodes")
MLsplits_LGcpp = MLsplits_LGcpp[order(MLsplits_LGcpp$Rnodes), ]
MLsplits_LGcpp



plot(tr)
nodelabels(text=MLsplits_LGcpp$splits)
tiplabels(tiprange_names)
title("LAGRANGE (C++) ML splits")

plot(tr)
pievals = as.matrix(MLsplits_LGcpp[,c("relprob","relprob2")])
nodelabels(pie=pievals, piecol=c("blue", "white"))
tiplabels(tiprange_names)
title("LAGRANGE (C++) split probs")






plot(tr)
nodelabels(text=MLsplits_LGcpp$LGnodes)
tiplabels()
title("LAGRANGE (C++) numbers")





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





#######################################################
# 
#######################################################
#' Title
#' 
#' Description
#' 
#' Details
#' 
#' 
#' 
#' Equation
#' 
#' The formula for the number of geographic states, based on the number of areas (\emph{N}),
#' is the sum of \emph{N} choose \emph{k}, from \emph{k}=1 to \emph{m}
#' (maximum range size) \deqn{s = \sum_{k=1}^{m}{N\choose k}}{s = sum(k=1...m)(N choose k)}
#' 
#' This equation assumes that the null range (a species lives in no areas, i.e. is extinct)
#' is not allowed. In the LAGRANGE program of \cite{ReeSmith2008}), the null range is included
#' in the transition matrix, and thus this is one more state.  This situation is represented in 
#' \code{numstates_from_numareas} by setting \code{include_null_range=TRUE}.
#'
#'
#' @param 
#' @param 
#' @param 
#' @return \code{} 
#' @export
#' @seealso \code{\link[stats]{convolve}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#'	 @cite FosterIdiots
#' @examples
#' 
LGcpp_splits_fn_to_table <- function(splits_fn)
	{
	splits = read.table(splits_fn)
	names(splits) = c("nodenum_LGpy", "splits", "relprob", "LnL")
	
	#splits = splits[, c("nodenum_LGpy", "splits", "LnL", "relprob")]
	splits$LnL = -1 * splits$LnL
	
	# Split the splits into left- and right- branch bottoms
	
	leftright = t(sapply(X=splits$splits, FUN=strsplit2, split="\\|"))
	row.names(leftright)=NULL
	leftright = as.data.frame(leftright)
	
	# "BB" means "branch bottom"
	names(leftright) = c("leftBB", "rightBB")
	leftright
	
	
	# Re-do nodenums, from bottom up
	nodenum_LGpy = splits$nodenum_LGpy
	uniq_nodenum_LGpy = rev(sort(unique(nodenum_LGpy)))
	nodenum_new = 1:length(uniq_nodenum_LGpy)
	
	nodenum_ORD1 = rep(NA, length(nodenum_LGpy))
	
	for (u in 1:length(uniq_nodenum_LGpy))
		{
		# Find matches to Python nodenums
		matchTF = nodenum_LGpy == uniq_nodenum_LGpy[u]
		
		# Produce new nodenums
		nodenum_ORD1[matchTF] = nodenum_new[u]
		}
	
	
	splits = cbind(nodenum_ORD1, splits[,c("nodenum_LGpy", "splits")], leftright, splits[,c("LnL", "relprob")])
	
	return(splits)
	}
