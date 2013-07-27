
require("ape")
require("rexpokit")
require("cladoRcpp")



#######################################################
# Reading/writing DIVA and LAGRANGE files
#######################################################


#######################################################
# get_pruningwise_nodenums
#######################################################
#' Get internal node numbers in pruningwise order
#' 
#' There are many ways of numbering nodes in a tree.  This returns a matrix
#' containing (column 1) R's native internal numbering scheme, and (column 2)
#' the node numbers in a pruningwise downpass.  Note that this is different
#' from \code{LAGRANGE}'s downpass ordering (see \code{\link{get_lagrange_nodenums}}).
#' 
#' @param tr A \code{\link[ape]{phylo}} tree object
#' @return \code{node_numbers_matrix} A matrix of node numbers
#' @export
#' @seealso \code{\link{get_lagrange_nodenums}}, \code{\link{prt}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
#' tmpdir = paste(extdata_dir, 
#' "/examples/Psychotria_M0/LGcpp/Psychotria_5.2.newick", sep="")
#' trfn = np(slashslash(tmpdir))
#' tr = read.tree(trfn)
#' node_numbers_matrix = get_pruningwise_nodenums(tr)
#' node_numbers_matrix
#' 
get_pruningwise_nodenums <- function(tr)
	{
	# Assumes tree is cladewise
	attr(tr,"order")
	
	# The edge matrix has ancestor (col 1) and descendant (col 2) nodes
	
	tr = reorder(tr, "pruningwise")
	attr(tr,"order")
	tr$edge
	
	#tr$edge[,2] = rev(tr$edge[,2])
	
	
	# Perhaps, the pairs of nodes are numbered in order
	i = 1
	edges_to_visit = seq(from=1, by=2, length.out=tr$Nnode)

	
	node_numbers_matrix = matrix(NA, nrow=tr$Nnode, ncol=2)
	node_rownum = 0
	
	for (i in edges_to_visit)
		{
		# Its sister is j 
		j <- i + 1
		
		node_rownum = node_rownum+1
		# Their LCA is the left node in the edge matrix
		node_numbers_matrix[node_rownum, 1] = tr$edge[node_rownum,1]
		
		# Another node number is the node_rownum
		node_numbers_matrix[node_rownum, 2] = node_rownum
		}
	
	return(node_numbers_matrix)
	}




# Get internal node numbers
#######################################################
# get_APE_nodenums
#######################################################
#' Get R internal node numbers
#' 
#' Utility function
#' 
#' @param tr A \code{\link[ape]{phylo}} tree object
#' @return \code{nodenums} A list of node numbers
#' @export
#' @seealso \code{\link{get_lagrange_nodenums}}, \code{\link{prt}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
#' 
get_APE_nodenums <- function(tr)
	{
	nodenums = length(tr$tip.label)+1 : tr$Nnode
	return(nodenums)
	}






#######################################################
# get_lagrange_nodenums
#######################################################
#' Get internal node numbers in LAGRANGE's downpass order
#' 
#' There are many ways of numbering nodes in a tree.  This returns a matrix
#' containing (column 1) R's native internal numbering scheme, and (column 2)
#' the node numbers in the downpass numbering used by C++ \code{LAGRANGE}, in particular in 
#' their .bgkey output file.  Note that this is different
#' from \code{\link[ape]{ape}}'s \code{pruningwise} downpass ordering (see \code{\link{get_pruningwise_nodenums}}).
#' 
#' The python version of LAGRANGE labels internal nodes differently (sigh), but they are
#' in the same order at least, so can just be renumbered from 1 to \code{tr$Nnode} to get them
#' to match the C++ \code{LAGRANGE} node numbering.
#' 
#' DIVA has yet a different node numbering scheme; see \code{\link{postorder_nodes_phylo4_return_table}}
#' 
#' @param tr A \code{\link[ape]{phylo}} tree object
#' @return \code{downpass_node_matrix} A matrix of node numbers
#' @export
#' @seealso \code{\link{get_pruningwise_nodenums}}, \code{\link{prt}}, \code{\link{postorder_nodes_phylo4_return_table}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
#' tmppath = paste(extdata_dir, 
#' "/examples/Psychotria_M0/LGcpp/Psychotria_5.2.newick", sep="")
#' trfn = np(slashslash(tmppath))
#' tr = read.tree(trfn)
#' downpass_node_matrix = get_lagrange_nodenums(tr)
#' downpass_node_matrix
#' 
#' 
#' downpass_node_matrix = get_lagrange_nodenums(tr)
#' downpass_node_matrix = downpass_node_matrix[order(downpass_node_matrix[,2]), ]
#' plot(tr)
#' nodelabels(node=20:37, downpass_node_matrix[,1])
#' tiplabels(1:19)
#' 
#' plot(tr)
#' nodelabels(node=20:37, downpass_node_matrix[,2])
#' tiplabels(1:19)
#' 
#' downpass_node_matrix = get_lagrange_nodenums(tr)
#' downpass_node_matrix = downpass_node_matrix[order(downpass_node_matrix[,1]), ]
#' plot(tr)
#' nodelabels(node=20:37, downpass_node_matrix[,1])
#' tiplabels(1:19)
#' 
#' # THIS WORKS
#' plot(tr)
#' nodelabels(node=20:37, downpass_node_matrix[,2])
#' tiplabels(1:19)
#' 
get_lagrange_nodenums <- function(tr)
	{
	defaults='
	downpass_node_matrix = get_lagrange_nodenums2(tr) 
	downpass_node_matrix
	plot(tr)
	tiplabels()
	nodelabels(text=downpass_node_matrix[,2], node=downpass_node_matrix[,1])
	'
	# Assumes tree is cladewise
	attr(tr,"order")
	
	# The edge matrix has ancestor (col 1) and descendant (col 2) nodes
	
	tr = reorder(tr, "pruningwise")
	attr(tr,"order")
	tr$edge
	
	#tr$edge[,2] = rev(tr$edge[,2])
	
	
	# Perhaps, the pairs of nodes are numbered in order
	i = 1
	edges_to_visit = seq(from=1, by=2, length.out=tr$Nnode)

	
	downpass_node_matrix = matrix(0, nrow=tr$Nnode, ncol=2)

	# Its sister is j 
	j <- i + 1

	anc = tr$edge[i,1]
	
	downpass_node_matrix = add_to_downpass_labels(tr, downpass_node_matrix, currnode=anc) 

	
	return(downpass_node_matrix)
	}


#######################################################
# LGpy_splits_fn_to_table
#######################################################
#' Get the ML splits per node, from Python LAGRANGE output
#' 
#' Python LAGRANGE outputs a list of splits and split probabilities for each node. This function
#' converts them to a table.
#' 
#' LAGRANGE outputs just the splits making up the top 95% of the probability, or 15 states, 
#' whichever comes first.
#' 
#' See \code{\link{LGpy_MLsplit_per_node}} for choosing the single ML split at each node, and 
#' see \code{\link{get_lagrange_nodenums}} for connecting these node numbers to APE node numbers.
#' 
#' @param splits_fn The filename of a Python LAGRANGE output file.
#' @return \code{splits} A data.frame containing the node numbers, splits, and split probabilities.
#' @export
#' @seealso \code{\link{get_lagrange_nodenums}}, \code{\link{LGpy_MLsplit_per_node}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @examples
#' test=1
#' 
LGpy_splits_fn_to_table <- function(splits_fn)
	{
	splits = read.table(splits_fn)
	names(splits) = c("nodenum_LGpy", "splits", "LnL", "relprob")
	
	# Split the splits into left- and right- branch bottoms
	
	leftright = t(sapply(X=splits$splits, FUN=strsplit2, split="\\|"))
	row.names(leftright)=NULL
	#leftright = as.data.frame(leftright)
	leftright = as.data.frame(leftright[,c(2,1)])

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

#######################################################
# LGcpp_splits_fn_to_table2
#######################################################
#' Get the ML splits per node, from Python LAGRANGE output
#' 
#' Python LAGRANGE outputs a list of splits and split probabilities for each node. This function
#' converts them to a table.
#' 
#' LAGRANGE outputs just the splits making up the top 95% of the probability, or 15 states, 
#' whichever comes first.
#' 
#' See \code{\link{LGpy_MLsplit_per_node}} for choosing the single ML split at each node, and 
#' see \code{\link{get_lagrange_nodenums}} for connecting these node numbers to APE node numbers.
#' 
#' @param splits_fn The filename of a Python LAGRANGE output file.
#' @return \code{splits} A data.frame containing the node numbers, splits, and split probabilities.
#' @export
#' @seealso \code{\link{get_lagrange_nodenums}}, \code{\link{LGpy_MLsplit_per_node}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @examples
#' test=1
#' 
LGcpp_splits_fn_to_table2 <- function(splits_fn)
	{
	splits = read.table(splits_fn)
	names(splits) = c("nodenum_LGcpp", "splits", "LnL", "relprob")
	
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



# Get the ML split
#######################################################
# LGpy_MLsplit_per_node
#######################################################
#' Get the ML splits per node, from a splits table
#' 
#' Given a table of splits probabilities from either \code{\link{LGpy_splits_fn_to_table}} or 
#' \code{\link{LGcpp_splits_fn_to_table}}, get the ML state for each node.
#' 
#' See \code{\link{get_lagrange_nodenums}} for connecting these node numbers to APE node numbers.
#' 
#' @param splits A data.frame containing the node numbers, splits, and split probabilities.
#' @return \code{MLsplits} A data.frame containing the node numbers, ML splits, and split probabilities.
#' @export
#' @seealso \code{\link{get_lagrange_nodenums}}, \code{\link{LGpy_splits_fn_to_table}}, \code{\link{LGcpp_splits_fn_to_table}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @examples
#' test=1
#' 
LGpy_MLsplit_per_node <- function(splits)
	{
	unique_nodenums = unique(splits$nodenum_ORD1)
	
	# Make table
	MLsplits = matrix(data=NA, nrow=length(unique_nodenums), ncol=ncol(splits))
	
	for (i in 1:length(unique_nodenums))
		{
		#print(i)
		# Subset the table
		tmptable = splits[splits$nodenum_ORD1==unique_nodenums[i],]
		
		# Find the (first) maxrow
		max_LnL = max(tmptable$LnL)
		nums = 1:nrow(tmptable)
		maxTF = tmptable$LnL == max_LnL
		num = nums[maxTF][1]
		tmprow = tmptable[num, ]
		
		MLsplits[i, ] = unlist(tmprow)
		}
	
	MLsplits = as.data.frame(MLsplits)
	names(MLsplits) = names(splits)
	MLsplits = unlist_df4(MLsplits)
	MLsplits = dfnums_to_numeric(MLsplits)
	
	# Get the probability of ALL of the the non-ML states
	relprob2 = 1-MLsplits$relprob
	MLsplits = cbind(MLsplits, relprob2)
	
	# Order according to ORD1 (hopefully R order)
	#MLsplits = MLsplits[order(MLsplits$nodenum_ORD1), ]

	
	return(MLsplits)
	}








# Get the ML state
#######################################################
# LGcpp_MLstate_per_node
#######################################################
#' Get the ML states per node, from a states table
#' 
#' Given a table of states probabilities from either \code{\link{LGcpp_states_fn_to_table}} or 
#' \code{\link{LGcpp_states_fn_to_table}}, get the ML state for each node.
#' 
#' See \code{\link{get_lagrange_nodenums}} for connecting these node numbers to APE node numbers.
#' 
#' @param states A data.frame containing the node numbers, states, and state probabilities.
#' @return \code{MLstates} A data.frame containing the node numbers, ML states, and state probabilities.
#' @export
#' @seealso \code{\link{get_lagrange_nodenums}}, \code{\link{LGcpp_states_fn_to_table}}, \code{\link{LGcpp_states_fn_to_table}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @examples
#' test=1
#' 
LGcpp_MLstate_per_node <- function(states)
	{
	unique_nodenums = unique(states$nodenum_ORD1)
	
	# Make table
	MLstates = matrix(data=NA, nrow=length(unique_nodenums), ncol=ncol(states))
	
	for (i in 1:length(unique_nodenums))
		{
		#print(i)
		# Subset the table
		tmptable = states[states$nodenum_ORD1==unique_nodenums[i],]
		
		# Find the (first) maxrow
		max_LnL = max(tmptable$LnL)
		nums = 1:nrow(tmptable)
		maxTF = tmptable$LnL == max_LnL
		num = nums[maxTF][1]
		tmprow = tmptable[num, ]
		
		MLstates[i, ] = unlist(tmprow)
		}
	
	MLstates = as.data.frame(MLstates)
	names(MLstates) = names(states)
	MLstates = unlist_df4(MLstates)
	MLstates = dfnums_to_numeric(MLstates)
	
	# Get the probability of ALL of the the non-ML states
	relprob2 = 1-MLstates$relprob
	MLstates = cbind(MLstates, relprob2)
	
	# Order according to ORD1 (hopefully R order)
	#MLstates = MLstates[order(MLstates$nodenum_ORD1), ]

	
	return(MLstates)
	}






#######################################################
# get_sister_node
#######################################################
#' Get the node sister to two nodes
#' 
#' Input two sister nodes, returns their "aunt".  Assumes a binary tree.
#' 
#' @param tr A \code{\link[ape]{phylo}} tree object.
#' @param nodepair A vector (length 2) with the node numbers of two nodes/tips.
#' @return \code{moms_sister} The aunt node.
#' @export
#' @seealso \code{\link{add_to_downpass_labels}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
#' 
get_sister_node <- function(tr, nodepair)
	{
	# Sort the pair
	nodepair = sort(nodepair)
	
	# get their moms
	anc1TF = tr$edge[,2] == nodepair[1]
	anc1 = tr$edge[anc1TF,1]
	anc2TF = tr$edge[,2] == nodepair[2]
	anc2 = tr$edge[anc2TF,1]
	
	# Check that they have the same mom
	if (anc1 != anc2)
		{
		stop("ERROR! Your nodes in nodepair are not sisters.")
		} else {
		anc = anc1
		}
	
	# Now, get mom's mom, i.e. the grandmom
	grandancTF = tr$edge[,2] == anc
	grandanc = tr$edge[grandancTF,1]
	
	# Get the moms_sister
	moms_and_sister = tr$edge[grandanc,2]
	moms_sister_TF = moms_and_sister != anc
	
	moms_sister = moms_and_sister[moms_sister_TF]
	return(moms_sister)
	}





#######################################################
# get_all_daughter_tips_of_a_node
#######################################################
#' Get all the daughter tips of a node
#' 
#' Like it says. Utility function.
#' 
#' @param nodenum The node to find
#' @param t A \code{\link[ape]{phylo}} tree object.
#' @return \code{temp_tips} The list of daughter tipnodes
#' @export
#' @seealso \code{\link{add_to_downpass_labels}}, \code{\link[ape]{extract.clade}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
#' 
get_all_daughter_tips_of_a_node <- function(nodenum, t)
	{
	subtree = extract.clade(t, nodenum)
	temp_tips = subtree$tip.label
	return(temp_tips)
	}





#######################################################
# add_to_downpass_labels
#######################################################
#' Iterate up and down a tree in C++ LAGRANGE downpass order
#' 
#' This is the utility function for \code{\link{get_lagrange_nodenums}}, which traces a 
#' tree down and up in C++ LAGRANGE's downpass order.
#'
#' This returns a matrix
#' containing (column 1) R's native internal numbering scheme, and (column 2)
#' the node numbers in a LAGRANGE downpass.  Note that this is different
#' from \code{LAGRANGE}'s downpass ordering (see \code{\link{get_lagrange_nodenums}}).
#' 
#' @param tr A \code{\link[ape]{phylo}} tree object.
#' @param downpass_node_matrix A matrix (\code{tr$Nnode} rows, 2 columns). Column 1 has R's native internal numbering scheme, 
#' and column 2 has the node numbers in a LAGRANGE downpass.
#' @param currnode The current node being viewed
#' @return \code{downpass_node_matrix} A matrix containing node numbers.
#' @export
#' @seealso \code{\link{get_lagrange_nodenums}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
#' 
add_to_downpass_labels <- function(tr, downpass_node_matrix, currnode) 
	{
	defaults='
	downpass_node_matrix = add_to_downpass_labels(tr, downpass_node_matrix, currnode=anc) 
	downpass_node_matrix
	'

	# Is the current node a tip? -- shouldn't happen!
	if (currnode <= length(tr$tip.label))
		{
		stop("ERROR: add_to_downpass_labels() should never have currnode be a tip!")
		}
	
	# Is the current node beneath tips?
	# Or, beneath nodes that have been named?
	daughtersTF = tr$edge[,1] == currnode
	daughters = tr$edge[daughtersTF,2]
	
	node_below_tips_TF = daughters <= length(tr$tip.label)
	node_below_names_TF = daughters %in% downpass_node_matrix[,1]
	sumTFs = (node_below_tips_TF + node_below_names_TF) > 0
	
	#print(downpass_node_matrix)
	
	# If the node is below tips/named things, then add
	if (all(sumTFs))
		{
		new_downpass_nodenum = max(downpass_node_matrix[,2], na.rm=TRUE) + 1
		unfilled_TF = downpass_node_matrix[,1] == 0
		
		# If you're full up, just return
		if (sum(unfilled_TF) == 0)
			{
			return(downpass_node_matrix)
			}
		

		unfilled_nums = (1:nrow(downpass_node_matrix))[unfilled_TF]
		
		# Fill first empty row
		downpass_node_matrix[unfilled_nums[1], ] = c(currnode, new_downpass_nodenum)
		
		# Then move DOWN
		
		# UNLESS you're at the root
		currnode_in_nodes_TF = currnode %in% tr$edge[,2]
		if (currnode_in_nodes_TF == FALSE)
			{
			# You're at the root
			return(downpass_node_matrix)
			}
		
		ancTF = tr$edge[,2] == currnode
		newanc = tr$edge[ancTF,1][1]
		downpass_node_matrix = add_to_downpass_labels(tr, downpass_node_matrix, currnode=newanc) 
		
		} else {
		# If not, then iterate UP to unfilled sisters
		
		for (i in 1:length(daughters[sumTFs==FALSE]))
			{
			daughter = daughters[sumTFs==FALSE][i]
			downpass_node_matrix = add_to_downpass_labels(tr, downpass_node_matrix, currnode=daughter) 
			}
		}
	return(downpass_node_matrix)
	}






#######################################################
# get_MLsplitprobs_from_results
#######################################################
#' Extract the ML probs for the base of each branch above a split
#' 
#' This function takes a BioGeoBEARS results_object from a ML search,
#' extracts the downpass and uppass likelihoods of the data for each
#' possible state at the base of each left and right branch, and 
#' produces the ML ancestral split estimates for the bottom of each branch.
#' 
#' @param results_object The results from a BioGeoBEARS ML search.
#' @return results_object with results_object$ML_marginal_prob_each_split_at_branch_bottom_BELOW_node added
#' @export
#' @seealso \code{\link{get_lagrange_nodenums}}, \code{\link{LGpy_splits_fn_to_table}}, \code{\link{LGcpp_splits_fn_to_table}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @examples
#' test=1
#' 
get_MLsplitprobs_from_results <- function(results_object)
	{
	#######################################################
	# Get the marginal probs of the splits (global ML probs, not local)
	# (These are marginal, rather than joint probs; but not local optima)
	#######################################################
	
	ML_marginal_prob_each_split_at_branch_bottom_BELOW_node = results_object$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS * results_object$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS

	results_object$ML_marginal_prob_each_split_at_branch_bottom_BELOW_node = ML_marginal_prob_each_split_at_branch_bottom_BELOW_node / rowSums(ML_marginal_prob_each_split_at_branch_bottom_BELOW_node)
	
	results_object$ML_marginal_prob_each_split_at_branch_bottom_BELOW_node
	rowSums(results_object$ML_marginal_prob_each_split_at_branch_bottom_BELOW_node)
		
	return(results_object)
	}




#######################################################
# get_leftright_nodes_matrix_from_results
#######################################################
#' Make a table of the Right and Left nodes descending from each node
#' 
#' This table shows the Right, then Left, descendant nodenums for each node. This
#' gets used later to plot splits at corners.
#' 
#' @param tr An ape phylo object
#' @param results_object The results from a BioGeoBEARS ML search.
#' @param nodes A list of internal node numbers for tree \code{tr}.
#' @return leftright_nodes_matrix A table with the Right, the Left, nodes
#' @export
#' @seealso \code{\link{get_lagrange_nodenums}}, \code{\link{LGpy_splits_fn_to_table}}, \code{\link{LGcpp_splits_fn_to_table}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @examples
#' test=1
#' 
get_leftright_nodes_matrix_from_results <- function(tr, results_object, nodes)
	{
	#######################################################
	# Get the marginal probs of the splits (global ML probs, not local)
	# (These are marginal, rather than joint probs; but not local optima)
	#######################################################
		
	# Make a splits table
	tr_table = prt(tr, printflag=FALSE)
	daughter_nds = tr_table$daughter_nds
	daughter_nds = daughter_nds[nodes]
	daughter_nds
	leftright_nodes_matrix = matrix(data=unlist(daughter_nds), ncol=2, byrow=TRUE)
	leftright_nodes_matrix = as.data.frame(leftright_nodes_matrix)
	names(leftright_nodes_matrix) = c("right", "left")
			
	return(leftright_nodes_matrix)
	}




#######################################################
# LGcpp_splits_fn_to_table
#######################################################
#' Get the ML splits per node, from C++ LAGRANGE output
#' 
#' C++ \code{LAGRANGE} outputs a list of splits and split probabilities for each node. This function
#' converts them to a table.
#' 
#' \code{LAGRANGE} outputs just the splits making up the top 95% of the probability, or 15 states, 
#' whichever comes first.
#' 
#' See \code{\link{LGpy_MLsplit_per_node}} for choosing the single ML split at each node, and 
#' see \code{\link{get_lagrange_nodenums}} for connecting these node numbers to APE node numbers.
#' 
#' @param splits_fn The filename of a C++ \code{LAGRANGE} output file.
#' @return \code{splits} A data.frame containing the node numbers, splits, and split probabilities.
#' @export
#' @seealso \code{\link{get_lagrange_nodenums}}, \code{\link{LGpy_MLsplit_per_node}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @examples
#' test=1
#' 
#' # splits_fn = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/
#' # examples/Psychotria_M0/LAGRANGE_C++/Psychotria_M0_lgcpp_out_splits00001.txt"
#' # LGcpp_splits_fn_to_table(splits_fn)
LGcpp_splits_fn_to_table <- function(splits_fn)
	{
	splits = read.table(splits_fn)
	names(splits) = c("nodenum_LGcpp", "splits", "relprob", "LnL")
	splits = splits[ , c("nodenum_LGcpp", "splits", "LnL", "relprob")]
	splits$LnL = -1 * splits$LnL
	
	# Split the splits into left- and right- branch bottoms
	leftright = t(sapply(X=splits$splits, FUN=strsplit2, split="\\|"))
	row.names(leftright)=NULL
	leftright = as.data.frame(leftright[,c(2,1)])
	
	# "BB" means "branch bottom"
	names(leftright) = c("leftBB", "rightBB")
	leftright
	
	
	# Re-do nodenums, from bottom up
	nodenum_LGcpp = splits$nodenum_LGcpp
	uniq_nodenum_LGcpp = rev(sort(unique(nodenum_LGcpp)))
	nodenum_new = 1:length(uniq_nodenum_LGcpp)
	
	nodenum_ORD1 = rep(NA, length(nodenum_LGcpp))
	
	for (u in 1:length(uniq_nodenum_LGcpp))
		{
		# Find matches to Python nodenums
		matchTF = nodenum_LGcpp == uniq_nodenum_LGcpp[u]
		
		# Produce new nodenums
		nodenum_ORD1[matchTF] = nodenum_new[u]
		}
	
	splits = cbind(nodenum_ORD1, splits[,c("nodenum_LGcpp", "splits")], leftright, splits[,c("LnL", "relprob")])
	
	return(splits)
	}




#######################################################
# LGcpp_states_fn_to_table
#######################################################
#' Get the ML states per node, from C++ LAGRANGE output
#' 
#' C++ \code{LAGRANGE} outputs a list of states and state probabilities for each node. This function
#' converts them to a table.
#' 
#' \code{LAGRANGE} outputs just the states making up the top 95% of the probability, or 15 states, 
#' whichever comes first.
#' 
#' See \code{\link{LGcpp_MLstate_per_node}} for choosing the single ML state at each node, and 
#' see \code{\link{get_lagrange_nodenums}} for connecting these node numbers to APE node numbers.
#' 
#' @param states_fn The filename of a C++ \code{LAGRANGE} output file.
#' @return \code{states} A data.frame containing the node numbers, states, and state probabilities.
#' @export
#' @seealso \code{\link{get_lagrange_nodenums}}, \code{\link{LGcpp_MLstate_per_node}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @examples
#' test=1
#' 
#' # states_fn = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/
#' # inst/extdata/examples/Psychotria_M0/LAGRANGE_C++/
#' # Psychotria_M0_lgcpp_out_states00001.txt"
#' # LGcpp_states_fn_to_table(states_fn)
LGcpp_states_fn_to_table <- function(states_fn)
	{
	states = read.table(states_fn)
	names(states) = c("nodenum_LGcpp", "states", "relprob", "LnL")
	states = states[ , c("nodenum_LGcpp", "states", "LnL", "relprob")]
	states$LnL = -1 * states$LnL
	
	# state the states into left- and right- branch bottoms
	#leftright = t(sapply(X=states$states, FUN=strsplit2, split="\\|"))
	#row.names(leftright)=NULL
	#leftright = as.data.frame(leftright[,c(2,1)])
	
	# "BB" means "branch bottom"
	#names(leftright) = c("leftBB", "rightBB")
	#leftright
	
	
	# Re-do nodenums, from bottom up
	nodenum_LGcpp = states$nodenum_LGcpp
	uniq_nodenum_LGcpp = rev(sort(unique(nodenum_LGcpp)))
	nodenum_new = 1:length(uniq_nodenum_LGcpp)
	
	nodenum_ORD1 = rep(NA, length(nodenum_LGcpp))
	
	for (u in 1:length(uniq_nodenum_LGcpp))
		{
		# Find matches to Python nodenums
		matchTF = nodenum_LGcpp == uniq_nodenum_LGcpp[u]
		
		# Produce new nodenums
		nodenum_ORD1[matchTF] = nodenum_new[u]
		}
	
	states = cbind(nodenum_ORD1, states[,c("nodenum_LGcpp", "states")], states[,c("LnL", "relprob")])
	
	return(states)
	}








#######################################################
# Convert a prt_tree into a phylo4 object
#######################################################
#' prt_tree_to_phylo4
#' 
#' Converts a tree table (a prt_tree from the function \code{\link{prt}}, which prints trees to tables) to a
#' \code{\link[phylobase]{phylobase}} \code{\link[phylobase]{phylo4}} tree object.
#' 
#' @param prt_tr A prt_tree from the function \code{\link{prt}}.
#' @return \code{newtr} A \code{\link[phylobase]{phylobase}} \code{\link[phylobase]{phylo4}} tree object.
#' @export
#' @seealso \code{\link[phylobase]{phylo4}}, \code{\link{prt}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
#' 
prt_tree_to_phylo4 <- function(prt_tr)
	{
	newtr = cbind(prt_tr$label, prt_tr$node, prt_tr$ancestor, prt_tr$edge.length, prt_tr$node.type)
	newtr = adf(newtr)
	names(newtr) = c("label", "node", "ancestor", "edge.length", "node.type")
	
	return(newtr)
	}


#######################################################
# label_nodes_postorder_phylo3
#######################################################
#' Add postorder node number labels to a phylo3 tree object.
#' 
#' Adds \code{\link[phylobase]{phylobase}} \code{\link[phylobase]{phylo4}} postorder node number labels to a 
#' \code{\link[ape]{phylo}} tree object.
#' 
#' @param tr2 \code{\link[ape]{phylo}} tree object.
#' @return \code{tr2} A \code{\link[ape]{phylo}} tree object with node labels added.
#' @export
#' @seealso \code{\link[ape]{phylo}}, \code{\link[phylobase]{phylo4}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
#' 
label_nodes_postorder_phylo3 <- function(tr2)
	{
	require(phylobase)
	tr4 = as(tr2, "phylo4")
	tr5 = reorder(tr4, "postorder")
	
	tmp_edge = attr(tr5, "edge")
	tmp_edge2 = tmp_edge[tmp_edge[,2] > length(tr2$tip.label), ]
	tmp_edge2[,2]
	
	nodenums_to_change = tmp_edge2[,2]-length(tr2$tip.label)
	tr2$node.label[nodenums_to_change] = 1:tr2$Nnode
	#plot(tr2)
	#nodelabels(tr2$node.label)
	
	return(tr2)
	}

#######################################################
# postorder_nodes_phylo4_return_table
#######################################################
#' Get a table of node numbers, including DIVA node numbers
#' 
#' Various programs (annoyingly) label internal nodes in different ways.  This function
#' shows the corresponding node numbers for several different systems. This table can
#' then be used to translate, when the user wishes to plot the output from various 
#' programs on the nodes of a tree.  In particular, the last column contains the DIVA 
#' node-numbering scheme (\cite{Ronquist1996_DIVA}, \cite{Ronquist_1997_DIVA}).
#' 
#' There are many ways of numbering nodes in a tree.  This returns a matrix
#' containing (column 1) R's native internal numbering scheme, and (column 2)
#' the node numbers in the downpass numbering used by C++ \code{LAGRANGE}, in particular in 
#' their .bgkey output file.  Note that this is different
#' from \code{\link[ape]{ape}}'s \code{pruningwise} downpass ordering (see \code{\link{get_pruningwise_nodenums}}).
#' 
#' The python version of LAGRANGE labels internal nodes differently (sigh), but they are
#' in the same order at least, so can just be renumbered from 1 to \code{tr$Nnode} to get them
#' to match the C++ \code{LAGRANGE} node numbering.
#' 
#' DIVA has yet a different node numbering scheme; see \code{\link{postorder_nodes_phylo4_return_table}}
#'
#' @param tr4 A tree object in \code{\link[ape]{phylo}} or \code{\link[phylobase]{phylo4}} format.
#' @return \code{postorder_table} A data.frame showing the various corresponding node numbers.
#' @export
#' @seealso \code{\link{get_pruningwise_nodenums}}, \code{\link{get_lagrange_nodenums}}, \code{\link{prt}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#'	 @cite Ronquist1996_DIVA
#'	 @cite Ronquist_1997_DIVA
#' @examples
#' test=1
#' 
postorder_nodes_phylo4_return_table <- function(tr4)
	{
	require(phylobase)
	
	# Check if phylo instead of phylo4
	if (is(tr4, "phylo") == TRUE)
		{
		tr4 = as(tr4, "phylo4")
		
		# Force this, to get "descendants" to work properly;
		# see r-sig-phylo
		tr4@order = "unknown"
		}
	
	# Find root row
	rootnode = nodeId(tr4, type="root")

	tipnames = descendants(phy=tr4, node=rootnode, type="tips")
	
	tr5 = reorder(tr4, "postorder")
	tmp_edge = attr(tr5, "edge")
	tmp_edge2 = tmp_edge[tmp_edge[,2] > length(tipnames), ]
	tmp_edge2[,2]
	
	nodenums = tmp_edge2[,2]
	internal_nodenums = tmp_edge2[,2]-length(tipnames)
	postorder = 1:length(nodenums)
	
	
	# Get postorder numbering, DIVA-style
	ntips = summary(tr4)$nb.tips
	
	DIVA_postorder = ntips+postorder
	
	postorder_table = cbind(nodenums, internal_nodenums, postorder, DIVA_postorder)
	postorder_table = adf2(postorder_table)
	
	return(postorder_table)	
	}

#######################################################
# traverse_up
#######################################################
#' Traverse the tree from node up to the tips
#' 
#' This is a utility function for \code{\link{nodenums_bottom_up}}.
#'
#' @param tr4 A tree object in \code{\link[phylobase]{phylo4}} format.
#' @param startnode The node number to start the uppass at.
#' @param traverse_records A list of the nodes visited.
#' @return \code{traverse_records} 
#' @export
#' @seealso \code{\link[phylobase]{phylo4}}, 
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
#' 
traverse_up <- function(tr4, startnode, traverse_records)
	{
	tmp_children = children(tr4, startnode)
	traverse_records$childnum = traverse_records$childnum + 1
	traverse_records$newnodes[startnode] = traverse_records$childnum

	# print if desired
	#cat(paste(startnode, "	", traverse_records$childnum, "\n", sep=""))
	
	# if at a tip, return down the tree
	if (length(tmp_children) == 0)
		{
		return(traverse_records)
		}
	
	for (child in tmp_children)
		{	
		traverse_records = traverse_up(tr4, child, traverse_records)
		}
	return(traverse_records)
	}






#######################################################
# nodenums_bottom_up
#######################################################
#' Assign node labels in bottom-up, left-first format (as in e.g. r8s)
#' 
#' This function assigns node numbers by tracing up from the root.  This corresponds to the 
#' node numbers in e.g. r8s (\cite{Sanderson_2003_r8s}).
#'
#' @param tr A tree object in \code{\link[ape]{phylo}} format.
#' @return \code{traverse_records} 
#' @export
#' @seealso \code{\link[phylobase]{phylo4}}, 
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'   @cite Sanderson_2003_r8s
#'   @cite r8s
#'   @cite Marazzi_etal_Sanderson_2012_r8s_morph
#' @examples
#' test=1
#' 
nodenums_bottom_up <- function(tr)
	{
	tr4 = as(tr, "phylo4")
	tr4@order = "unknown"
	
	numnodes = nNodes(tr4) + nTips(tr4)
	newnodes = matrix(NA, nrow=numnodes, ncol=1)
	
	# Find root row
	rootnode = nodeId(tr4, type="root")
	#newnodes[rootnode] = 1
	
	startnode = rootnode
	childnum = 0
	traverse_records = c()
	traverse_records$childnum = 0
	traverse_records$newnodes = newnodes
	traverse_records = traverse_up(tr4, startnode, traverse_records)
	
	newnodes = traverse_records$newnodes
	old2new_nodes_table = cbind(1:numnodes, newnodes)
	old2new_nodes_table = adf(old2new_nodes_table)
	names(old2new_nodes_table) = c("orig_nodenums", "LR_nodenums")
	
	return(old2new_nodes_table)
	}



#######################################################
# get_nodenums
#######################################################
#' Get the unique node numbers in a tree
#' 
#' This is a utility function for \code{\link{get_nodenum_structural_root}}.
#'
#' @param t A tree object in \code{\link[ape]{phylo}} format.
#' @return \code{ordered_nodenames} The node numbers, in order.
#' @export
#' @seealso \code{\link[ape]{phylo}}, \code{\link{get_nodenum_structural_root}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' blah = 1
# this returns the NUMBERS identifying each node
get_nodenums <- function(t)
	{
	# get just the unique node numbers from the edge list (left column: start node; right column: end node):
	nodenames = unique(c(t$edge))
	ordered_nodenames = nodenames[order(nodenames)]
	return(ordered_nodenames)
	}

#######################################################
# get_nodenum_structural_root
#######################################################
#' Gets the root node 
#' 
#' This function gets the root node by finding the node not in the descendants list (edge[,2]). This
#' may be more reliable than e.g. assuming length(tr$tip.label)+1.
#'
#' @param t A tree object in \code{\link[ape]{phylo}} format.
#' @param print_nodenum Print the node numbers as you go through the list? Default FALSE.
#' @return \code{root_nodenums_list} 
#' @export
#' @seealso \code{\link[ape]{phylo}}, \code{\link{get_nodenums}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' blah=1
#' 
get_nodenum_structural_root <- function(t, print_nodenum=FALSE)
	{
	#numnodes = length(t$tip.label) + length(t$node.label)
	#ordered_nodes = 1:length(numnodes)
	
	ordered_nodes = get_nodenums(t)

	root_nodenums_list = c()
	for (n in 1:length(ordered_nodes))
		{
		tmpnode = ordered_nodes[n]
		if (tmpnode %in% t$edge[,2])
			{
			blah = TRUE
			}
		else
			{
			if (print_nodenum == TRUE)
				{
				cat("get_nodenum_structural_root(): Root nodenum = ", tmpnode, sep="")
				}
			root_nodenums_list = c(root_nodenums_list, tmpnode)
			}
		}
	return(root_nodenums_list)
	}





# print tree in hierarchical format
#######################################################
# prt
#######################################################
#' Print tree in table format
#' 
#' Learning and using APE's tree structure can be difficult and confusing because much of the information is
#' implicit.  This function prints the entire
#' tree to a table, and makes much of the implicit information explicit.  It is not particularly fast, but
#' it is useful.
#'
#' See \url{http://ape.mpl.ird.fr/ape_development.html} for the official documentation of R tree objects.
#' 
#' @param t A \code{\link[ape]{phylo}} tree object.
#' @param printflag Should the table be printed to screen?  Default TRUE.
#' @param relabel_nodes Manually renumber the internal nodes, if desired. Default FALSE.
#' @param time_bp_digits The number of digits to print in the time_bp (time before present) column. Default=7.
#' @param add_root_edge Should a root edge be added?  Default \code{TRUE}.
#' @param get_tipnames Should the list of tipnames descending from each node be printed as a string in another column?  
#' This is slow-ish, but useful for matching up nodes between differing trees. Default \code{FALSE}.
#' @param fossils_older_than Tips that are older than \code{fossils_older_than} will be marked as \code{TRUE} in a column called \code{fossil}.
#' This is not currently set to 0, because Newick files can have slight precision issues etc. that mean not all tips quite come to zero.  You 
#' can attempt to fix this with \code{\link{average_tr_tips}} (but make sure you do not inappropriately average in fossils!!).
#' @return \code{dtf} A \code{\link[base]{data.frame}} holding the table. (Similar to the printout of a \code{\link[phylobase]{phylo4}} object.)
#' @export
#' @seealso \code{\link[ape]{phylo}}, \code{\link{average_tr_tips}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://ape.mpl.ird.fr/ape_development.html}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
#' 
prt <- function(t, printflag=TRUE, relabel_nodes = FALSE, time_bp_digits=7, add_root_edge=TRUE, get_tipnames=FALSE, fossils_older_than=0.6)
	{
	# assemble beginning table
	
	# check if internal node labels exist
	if ("node.label" %in% attributes(t)$names == FALSE)
		{
		rootnum = get_nodenum_structural_root(t)
		
		new_node_labels = paste("inNode", rootnum:(rootnum+t$Nnode-1), sep="")
		t$node.label = new_node_labels
		}
	
	# or manually relabel the internal nodes, if desired
	if (relabel_nodes == TRUE)
		{
		rootnum = get_nodenum_structural_root(t)
		
		new_node_labels = paste("inNode", rootnum:(rootnum+t$Nnode-1), sep="")
		t$node.label = new_node_labels
		}
	
	labels = c(t$tip.label, t$node.label)
	ordered_nodenames = get_nodenums(t)
	#nodenums = 1:length(labels)
	node.types1 = rep("tip", length(t$tip.label))
	node.types2 = rep("internal", length(t$node.label))
	node.types2[1] = "root"
	node.types = c(node.types1, node.types2)
	
	# These are the index numbers of the edges below each node
	parent_branches = get_indices_where_list1_occurs_in_list2(ordered_nodenames, t$edge[,2])
	#parent_edges = parent_branches
	brlen_to_parent = t$edge.length[parent_branches]
	
	parent_nodes = t$edge[,1][parent_branches]
	daughter_nodes = lapply(ordered_nodenames, get_daughters, t)
	
	# print out the structural root, if desired
	root_nodenum = get_nodenum_structural_root(t)
	tmpstr = paste("prt(t): root=", root_nodenum, "\n", sep="")
	prflag(tmpstr, printflag=printflag)
	
	levels_for_nodes = unlist(lapply(ordered_nodenames, get_level, t))
	#tmplevel = get_level(23, t)
	#print(tmplevel)
	
	
	#height above root
	hts_at_end_of_branches_aka_at_nodes = t$edge.length
	hts_at_end_of_branches_aka_at_nodes = get_all_node_ages(t)
	h = hts_at_end_of_branches_aka_at_nodes

	# times before present, below (ultrametric!) tips
	# numbers are positive, i.e. in millions of years before present
	#                       i.e. mybp, Ma
	times_before_present = get_max_height_tree(t) - h

	
	# fill in the ages of each node for the edges
	edge_ages = t$edge
	edge_ages[,1] = h[t$edge[,1]]	# bottom of branch
	edge_ages[,2] = h[t$edge[,2]]	# top of branch


	# fill in the times before present of each node for the edges
	edge_times_bp = t$edge
	edge_times_bp[,1] = times_before_present[t$edge[,1]]	# bottom of branch
	edge_times_bp[,2] = times_before_present[t$edge[,2]]	# top of branch
	
	
	# If desired, get the list of all tipnames descended from a node, in alphabetical order
	if (get_tipnames == TRUE)
		{
		# Make the empty list
		list_of_clade_members_lists = rep(list(NA), length(ordered_nodenames))
		
		# Tips have only one descendant
		list_of_clade_members_lists[1:length(t$tip.label)] = t$tip.label
		list_of_clade_members_lists
		
		
		nontip_nodenums = (length(t$tip.label)+1) : length(ordered_nodenames)
		if (length(nontip_nodenums) > 1)
			{
			# More than 1 node
			nontip_nodenames = ordered_nodenames[nontip_nodenums]
			nontip_cladelists = sapply(X=nontip_nodenames, FUN=get_all_daughter_tips_of_a_node, t=t)
			nontip_cladelists
			
			nontip_cladelists_alphabetical = sapply(X=nontip_cladelists, FUN=sort)
			nontip_cladelists_alphabetical
			
			nontip_cladelists_alphabetical_str = sapply(X=nontip_cladelists_alphabetical, FUN=paste, collapse=",")
			nontip_cladelists_alphabetical_str
			
			# Store the results
			list_of_clade_members_lists[nontip_nodenums] = nontip_cladelists_alphabetical_str
			list_of_clade_members_lists
			} else {
			# Just one node
			nontip_nodenames = ordered_nodenames[nontip_nodenums]
			nontip_cladelists = sapply(X=nontip_nodenames, FUN=get_all_daughter_tips_of_a_node, t=t)
			nontip_cladewords = unlist(sapply(X=nontip_cladelists, FUN=strsplit, split=","))
			
			nontip_cladelists_alphabetical = sort(nontip_cladewords)
			nontip_cladelists_alphabetical
			
			nontip_cladelists_alphabetical_str = paste(nontip_cladelists_alphabetical, collapse=",", sep="")
			nontip_cladelists_alphabetical_str
			
			# Store the results
			list_of_clade_members_lists[nontip_nodenums] = nontip_cladelists_alphabetical_str
			list_of_clade_members_lists			
			}
			
		}

	
	# Add fossils TRUE/FALSE column.  You can turn this off with fossils_older_than=NULL.
	fossils = times_before_present > fossils_older_than

	# Obviously, internal nodes are irrelevant and should be NA
	tmpnodenums = (length(t$tip.label)+1) : ( length(t$tip.label) + t$Nnode )
	fossils[tmpnodenums] = NA
	
	if (get_tipnames == FALSE)
		{
		# Don't put in the list of clade names
		tmpdtf = cbind(1:length(ordered_nodenames), ordered_nodenames, levels_for_nodes, node.types, parent_branches, brlen_to_parent, parent_nodes, daughter_nodes, h, round(times_before_present, digits=time_bp_digits), fossils, labels)
		
		dtf = as.data.frame(tmpdtf, row.names=NULL)
		# nd = node
		
		# edge.length is the same as brlen_2_parent
		names(dtf) = c("node", "ord_ndname", "node_lvl", "node.type", "parent_br", "edge.length", "ancestor", "daughter_nds", "node_ht", "time_bp", "fossils", "label")
		
		# convert the cols from class "list" to some natural class
		dtf = unlist_dtf_cols(dtf, printflag=FALSE)
		} else {
		# Put in the list of clade names
		tmpdtf = cbind(1:length(ordered_nodenames), ordered_nodenames, levels_for_nodes, node.types, parent_branches, brlen_to_parent, parent_nodes, daughter_nodes, h, round(times_before_present, digits=time_bp_digits), fossils, labels, list_of_clade_members_lists)
		
		dtf = as.data.frame(tmpdtf, row.names=NULL)
		# nd = node
		
		# edge.length is the same as brlen_2_parent
		names(dtf) = c("node", "ord_ndname", "node_lvl", "node.type", "parent_br", "edge.length", "ancestor", "daughter_nds", "node_ht", "time_bp", "fossils", "label", "tipnames")
		
		# convert the cols from class "list" to some natural class
		dtf = unlist_dtf_cols(dtf, printflag=FALSE)		
		}
	

	
	
	
	
	# Add the root edge, if desired
	# (AND, only if t$root.edge exists)
	if ( (add_root_edge == TRUE) && (!is.null(t$root.edge)) )
		{
		root_row_TF = dtf$node.type == "root"
		root_edge_length = t$root.edge
		
		# Stick in this edge length
		dtf$edge.length[root_row_TF] = root_edge_length
		
		# Add the root edge length to all node heights
		dtf$node_ht = dtf$node_ht + root_edge_length
		}
	
	# print if desired
	prflag(dtf, printflag=printflag)
	
	#tree_strings = c()
	#root_str = get_node_info(root_nodenum, t)
	return(dtf)
	}

















##########################################
# LAGRANGE (Python version) utilities
##########################################


#######################################################
# parse_lagrange_output_old
#######################################################
#' Parse the output file from python \code{LAGRANGE} -- older version
#' 
#' Parse the output of a C++ \code{LAGRANGE} run.  
#' 
#' This function parses the output of \code{LAGRANGE}, obtained by a command such as the following, run at a UNIX/Mac
#' Terminal command line.  This is an older version useful for automating processing of
#' many files.
#' 
#' \code{cd /Users/nick/Desktop/__projects/_2011-07-15_Hannah_spider_fossils/_data/lagrange_for_nick}
#' 
#' \code{./lagrange_cpp palp_no_Lacun_v1_2nd387.lg > lagrange_results_v1_2nd387.txt}
#' 
#' C++ LAGRANGE can be obtained at \url{https://code.google.com/p/lagrange/}
#' 
#' @param outfn The C++ \code{LAGRANGE} output text file.
#' @param results_dir The directory \code{outfn} is in.
#' @param new_splits_fn Should a text file containing a table of the splits and their probabilities be output? Default \code{TRUE}.
#' @param new_states_fn Should a text file containing a table of the splits and their probabilities be output? Default \code{TRUE}, 
#' unlike python \code{LAGRANGE}, C++ \code{LAGRANGE} \emph{will} output the states at the nodes.
#' @param filecount The starting number for the filecount (relevant if one is processing many files).
#' @return sumstats A \code{\link[base]{data.frame}} containing the summary statistics (LnL, d and e rates, etc.)  The splits
#' filename is output to screen.
#' @export
#' @seealso \code{\link{get_lagrange_nodenums}}, \code{\link{LGpy_splits_fn_to_table}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @examples
#' test=1
#' 
parse_lagrange_output_old <- function(outfn, results_dir=getwd(), new_splits_fn = TRUE, new_states_fn = TRUE, filecount=0)
	{
	setup='
	outfn = "/Users/nick/Desktop/__projects/_2011-07-15_Hannah_spider_fossils/_data/lagrange_for_nick/_run_3rd226/lagrange_results_v1_3rd226.txt"
	results_dir = "/Users/nick/Desktop/__projects/_2011-07-15_Hannah_spider_fossils/_data/lagrange_for_nick/_results_3rd226/"
	'


	# Split the file into 1 file for each analysis
	tmplines = scan(outfn, what="character", sep="\n")
	
	results_dir = addslash(results_dir)
	
	# Data to store
	initial_log_likelihood = NULL
	dispersal_rates = NULL
	extinction_rates = NULL
	final_log_likelihood = NULL
	#ancestral_splits_array = NULL
	#ancestral_states_array = NULL
	
	# Go through the lines
	#filecount = 0
	for (i in 1:length(tmplines))
		{
		tmpline = tmplines[i]
		
		# Start a new output file
		if (grepl("starting likelihood calculations", tmpline) == TRUE)
			{
			filecount = filecount + 1
			cat("\nProcessing tree #", filecount, sep="")
			tmpi = sprintf("%05.0f", filecount)
	
			new_prefix = get_fn_prefix(fn=get_path_last(path=outfn))
			splits_table_fn = np(paste(addslash(results_dir), new_prefix, "_splits", tmpi, ".txt", sep=""))
			states_table_fn = np(paste(addslash(results_dir), new_prefix, "_states", tmpi, ".txt", sep=""))
	
			# you need to clear these files
			#new_splits_fn = TRUE
			#new_states_fn = TRUE
			}
		
		# parse initial ln likelihood
		if (grepl("initial -ln likelihood", tmpline) == TRUE)
			{
			words = strsplit_whitespace(tmpline)
			tmpnum = -1 * as.numeric(words[length(words)])
			initial_log_likelihood = c(initial_log_likelihood, tmpnum)
			}
	
		# parse final ln likelihood
		if (grepl("final -ln likelihood", tmpline) == TRUE)
			{
			words = strsplit_whitespace(tmpline)
			tmpnum = -1 * as.numeric(words[length(words)])
			final_log_likelihood = c(final_log_likelihood, tmpnum)
			}
	
		# parse dispersal/extinction lines
		if (grepl("dis: ", tmpline) == TRUE)
			{
			words = strsplit_whitespace(tmpline)
			tmpnum = as.numeric(words[2])
			dispersal_rates = c(dispersal_rates, tmpnum)
	
			tmpnum = as.numeric(words[4])
			extinction_rates = c(extinction_rates, tmpnum)
			}
		
		# Ancestral splits
		if (grepl("Ancestral splits for:", tmpline) == TRUE)
			{
			# first line
			words = strsplit_whitespace(tmpline)
			tmpnum1 = as.numeric(words[length(words)])
			
			# second line
			for (j in ((i+1):(i+20)))
				{
				# Test if the line exists
				if (is.na(tmplines[j]))
					{
					break()
					}
	
				#cat("\nstart", tmplines[j], "end", sep="")
				if (tmplines[j] == "" || grepl("Ancestral", tmplines[j]) == TRUE || grepl("initializing", tmplines[j]) == TRUE)
					{
					break
					} else {
					words = strsplit_whitespace(tmplines[j])
					words[3] = extract_numbers(words[3])
					outwords = c(tmpnum1, words)
					outstr = list2str(outwords, spacer="\t")
	
					#ancestral_splits_array = rbind(ancestral_splits_array, tmpline)
					# Append data to end of file
					if (new_splits_fn == TRUE)
						{
						write(outstr, file=splits_table_fn, append=FALSE, sep="\n")
						new_splits_fn = FALSE
						} else {
						write(outstr, file=splits_table_fn, append=TRUE, sep="\n")
						}
					}
				}
			}
	
		# Ancestral states
		if (grepl("Ancestral states for:", tmpline) == TRUE)
			{
			# first line
			words = strsplit_whitespace(tmpline)
			tmpnum1 = as.numeric(words[length(words)])
			
			# second line
			for (j in ((i+1):(i+20)))
				{
				# Test if the line exists
				if (is.na(tmplines[j]))
					{
					break()
					}
	
				#cat("\nstart", tmplines[j], "end", sep="")
				if (tmplines[j] == "" || grepl("Ancestral", tmplines[j]) == TRUE || grepl("initializing", tmplines[j]) == TRUE)
					{
					break
					} else {
					words = strsplit_whitespace(tmplines[j])
					words[3] = extract_numbers(words[3])
					outwords = c(tmpnum1, words)
					outstr = list2str(outwords, spacer="\t")
					
					#ancestral_states_array = rbind(ancestral_states_array, tmpline)
					# Append data to end of file
					if (new_states_fn == TRUE)
						{
						write(outstr, file=states_table_fn, append=FALSE, sep="\n")
						new_states_fn = FALSE
						} else {
						write(outstr, file=states_table_fn, append=TRUE, sep="\n")
						}
					}
				}
			}
		#print(outstr)
		}
	cat("\n")

	cat("\nStates output to: ", splits_table_fn, "\n", sep="")
	cat("\nSplits output to: ", splits_table_fn, "\nn", sep="")
	
	sumstats = cbind(initial_log_likelihood, dispersal_rates, extinction_rates, final_log_likelihood, splits_table_fn, states_table_fn)
	sumstats = adf2(sumstats)
	
	outfn = np(paste(results_dir, "sumstats.txt", sep=""))
	write.table(x=sumstats, file=outfn, append=FALSE, quote=FALSE, sep="	", row.names=FALSE, col.names=TRUE)
	# write_table_good(sumstats, paste(results_dir, "sumstats.txt", sep=""))

	return(sumstats)
	}



#######################################################
# parse_lagrange_python_output_old
#######################################################
#' Parse the output file from python \code{LAGRANGE} -- old version
#' 
#' Parse the output of a python \code{LAGRANGE} output file.  This is 
#' an older version useful
#' for automating the parsing of a large number of files.
#' 
#' Python LAGRANGE is run from a UNIX/Terminal command-line
#' with a command such as "\code{python lagrangefilename.py}".  You will need to have the "lagrange" python directory in 
#' your working directory.
#' 
#' The input file can be obtained from \url{http://www.reelab.net/lagrange/configurator/index} (\cite{Ree2009configurator}).
#' 
#' Python comes installed on many machines, or can be downloaded from the Enthought Python Distribution 
#' (\url{https://www.enthought.com/products/epd/}).
#' 
#' @param outfn The python \code{LAGRANGE} output text file.
#' @param results_dir The directory \code{outfn} is in.
#' @param new_splits_fn Should a text file containing a table of the splits and their probabilities be output? Default \code{TRUE}.
#' @param new_states_fn Should a text file containing a table of the splits and their probabilities be output? Default \code{FALSE}, 
#' as I don't believe python \code{LAGRANGE} will output the states at the nodes (C++ \code{LAGRANGE} will, however).
#' @param filecount The starting number for the filecount (relevant if one is processing many files).
#' @return sumstats A \code{\link[base]{data.frame}} containing the summary statistics (LnL, d and e rates, etc.)  The splits
#' filename is output to screen.
#' @export
#' @seealso \code{\link{get_lagrange_nodenums}}, \code{\link{LGpy_splits_fn_to_table}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{https://code.google.com/p/lagrange/}
#' \url{https://www.enthought.com/products/epd/}
#' \url{http://www.reelab.net/lagrange/configurator/index}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @examples
#' test=1
#' 
parse_lagrange_python_output_old <- function(outfn="output.results.txt", results_dir=getwd(), new_splits_fn = TRUE, new_states_fn = FALSE, filecount=0)
	{
	setup='
	wd = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/Psychotria_M0/LAGRANGE_python/"
	setwd(wd)
	outfn = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/Psychotria_M0/LGpy/Psychotria_5.2_demo.results.txt"
	results_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/Psychotria_M0/LGpy/"
	new_splits_fn = TRUE
	new_states_fn = FALSE
	filecount=0
	parse_lagrange_python_output(outfn, results_dir, new_splits_fn = TRUE, new_states_fn = FALSE, filecount=0)
	'

	# Split the file into 1 file for each analysis
	tmplines = scan(outfn, what="character", sep="\n")
	
	results_dir = addslash(results_dir)
	
	# Data to store
	initial_log_likelihood = NULL
	dispersal_rates = NULL
	extinction_rates = NULL
	final_log_likelihood = NULL
	#ancestral_splits_array = NULL
	#ancestral_states_array = NULL
	
	# Go through the lines
	#filecount = 0
	for (i in 1:length(tmplines))
		{
		tmpline = tmplines[i]
		
		# Start a new output file
		if (grepl("Global ML at root node:", tmpline) == TRUE)
			{
			filecount = filecount + 1
			cat("\nProcessing tree #", filecount, sep="")
			tmpi = sprintf("%05.0f", filecount)
	
			new_prefix = get_fn_prefix(fn=get_path_last(path=outfn))
			splits_table_fn = np(paste(addslash(results_dir), new_prefix, "_splits", tmpi, ".txt", sep=""))
			states_table_fn = np(paste(addslash(results_dir), new_prefix, "_states", tmpi, ".txt", sep=""))
	
			# you need to clear these files
			#new_splits_fn = TRUE
			#new_states_fn = TRUE
			}
		
		# parse initial ln likelihood
		if (grepl("initial -ln likelihood", tmpline) == TRUE)
			{
			words = strsplit_whitespace(tmpline)
			tmpnum = -1 * as.numeric(words[length(words)])
			initial_log_likelihood = c(initial_log_likelihood, tmpnum)
			}
	
		# parse final ln likelihood
		if (grepl("-lnL =", tmpline) == TRUE)
			{
			words = strsplit_whitespace(tmpline)
			tmpnum = -1 * as.numeric(words[length(words)])
			final_log_likelihood = c(final_log_likelihood, tmpnum)
			}
	
		# parse dispersal line
		if (grepl("dispersal = ", tmpline) == TRUE)
			{
			words = strsplit_whitespace(tmpline)
			tmpnum = as.numeric(words[3])
			dispersal_rates = c(dispersal_rates, tmpnum)
			}
		# parse extinction line
		if (grepl("extinction = ", tmpline) == TRUE)
			{
			words = strsplit_whitespace(tmpline)
			tmpnum = as.numeric(words[3])
			extinction_rates = c(extinction_rates, tmpnum)
			}
		
		# Ancestral splits
		if (grepl("At node N", tmpline) == TRUE)
			{
			# first line, e.g. "At node N36:"
			words = strsplit_whitespace(tmpline)
			tmpnum3 = words[length(words)]
			tmpnum2 = gsub(pattern="N", replacement="", x=tmpnum3)
			tmpnum2 = gsub(pattern="\\:", replacement="", x=tmpnum2)
			tmpnum1 = as.numeric(tmpnum2)

			# second line (there are only 15 values reported, plus header
			for (j in ((i+1):(i+20)))
				{
				# Test if the line exists
				if (is.na(tmplines[j]))
					{
					break()
					}

				# Test if the line is header:    split     lnL     Rel.Prob
				if ( (grepl("split", tmplines[j]) == TRUE) && (grepl("lnL", tmplines[j]) == TRUE) )
					{
					next()
					}

	
				#cat("\nstart", tmplines[j], "end", sep="")
				if (tmplines[j] == "" || grepl("At node N", tmplines[j]) == TRUE || grepl("initializing", tmplines[j]) == TRUE)
					{
					break()
					} else {
					words = strsplit_whitespace(tmplines[j])
					
					# Fix splits by removing []
					words[1] = gsub(pattern="\\[", replacement="", x=words[1])
					words[1] = gsub(pattern="\\]", replacement="", x=words[1])
					
					# -LnL
					words[2] = -1* as.numeric(extract_numbers(words[2]))
					
					# Rel. prob.
					words[3] = extract_numbers(words[3])
					outwords = c(tmpnum1, words)
					outstr = list2str(outwords, spacer="\t")
	
					#ancestral_splits_array = rbind(ancestral_splits_array, tmpline)
					# Append data to end of file
					if (new_splits_fn == TRUE)
						{
						write(outstr, file=splits_table_fn, append=FALSE, sep="\n")
						new_splits_fn = FALSE
						} else {
						write(outstr, file=splits_table_fn, append=TRUE, sep="\n")
						}
					}
				}
			}
	
		# Ancestral states
		if (grepl("Ancestral states for:", tmpline) == TRUE)
			{
			# first line
			words = strsplit_whitespace(tmpline)
			tmpnum1 = as.numeric(words[length(words)])
			
			# second line
			for (j in ((i+1):(i+20)))
				{
				# Test if the line exists
				if (is.na(tmplines[j]))
					{
					break()
					}
	
				#cat("\nstart", tmplines[j], "end", sep="")
				if (tmplines[j] == "" || grepl("Ancestral", tmplines[j]) == TRUE || grepl("initializing", tmplines[j]) == TRUE)
					{
					break
					} else {
					words = strsplit_whitespace(tmplines[j])
					words[3] = extract_numbers(words[3])
					outwords = c(tmpnum1, words)
					outstr = list2str(outwords, spacer="\t")
					
					#ancestral_states_array = rbind(ancestral_states_array, tmpline)
					# Append data to end of file
					if (new_states_fn == TRUE)
						{
						write(outstr, file=states_table_fn, append=FALSE, sep="\n")
						new_states_fn = FALSE
						} else {
						write(outstr, file=states_table_fn, append=TRUE, sep="\n")
						}
					}
				}
			}
		#print(outstr)
		}
	cat("\n")

	cat("\nSplits output to: ", splits_table_fn, "\n", sep="")
	#cat("\nStates output (if LAGRANGE Python does this, which I don't think it does) to: ", splits_table_fn, "\n", sep="")
	if (new_states_fn == TRUE)
		{
		cat("\nStates output (if LAGRANGE Python does this, which I don't think it does) to: ", states_table_fn, "\n", sep="")
		} else {
		states_table_fn = NA
		}
	
	
	sumstats = cbind(initial_log_likelihood, dispersal_rates, extinction_rates, final_log_likelihood, splits_table_fn, states_table_fn)
	sumstats = adf2(sumstats)
	outfn = np(paste(results_dir, "/sumstats.txt", sep=""))
	write.table(x=sumstats, file=outfn, append=FALSE, quote=FALSE, sep="	", row.names=FALSE, col.names=TRUE)
	
	return(sumstats)
	}






#######################################################
# parse_lagrange_python_output
#######################################################
#' Parse the output file from python \code{LAGRANGE}
#' 
#' Parse the output of a python \code{LAGRANGE}.  
#' 
#' Python LAGRANGE is run from a UNIX/Terminal command-line
#' with a command such as "\code{python lagrangefilename.py}".  You will need to have the "lagrange" python directory in 
#' your working directory.
#' 
#' The input file can be obtained from \url{http://www.reelab.net/lagrange/configurator/index} (\cite{Ree2009configurator}).
#' 
#' Python comes installed on many machines, or can be downloaded from the Enthought Python Distribution 
#' (\url{https://www.enthought.com/products/epd/}).
#' 
#' @param outfn The python \code{LAGRANGE} output text file.
#' @param outputfiles Should parsed output be written to files? Default FALSE.
#' @param results_dir The directory \code{outfn} is in.
#' @param new_splits_fn Should a text file containing a table of the splits and their probabilities be output? Default \code{TRUE}.
#' @param new_states_fn Should a text file containing a table of the states and their probabilities be output? Default \code{FALSE}, 
#' as I don't believe python \code{LAGRANGE} will output the states at the nodes (C++ \code{LAGRANGE} will, however).
#' @param filecount The starting number for the filecount (relevant if one is processing many files).
#' @param append Should results be appended to preexisting file? (default \code{FALSE})
#' @return sumstats A \code{\link[base]{data.frame}} containing the summary statistics (LnL, d and e rates, etc.)  The splits
#' filename is output to screen.
#' @export
#' @seealso \code{\link{get_lagrange_nodenums}}, \code{\link{LGpy_splits_fn_to_table}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{https://code.google.com/p/lagrange/}
#' \url{https://www.enthought.com/products/epd/}
#' \url{http://www.reelab.net/lagrange/configurator/index}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @examples
#' test=1
#' 
parse_lagrange_python_output <- function(outfn="output.results.txt", outputfiles=FALSE, results_dir=getwd(), new_splits_fn = TRUE, new_states_fn = FALSE, filecount=0, append=FALSE)
	{
	setup='
	wd = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/Psychotria_M0/LAGRANGE_python/"
	setwd(wd)
	outfn = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/Psychotria_M0/LGpy/Psychotria_5.2_demo.results.txt"
	results_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/Psychotria_M0/LGpy/"
	new_splits_fn = TRUE
	new_states_fn = FALSE
	filecount=0
	parse_lagrange_python_output(outfn, results_dir, new_splits_fn = TRUE, new_states_fn = FALSE, filecount=0)
	'

	# Split the file into 1 file for each analysis
	tmplines = scan(outfn, what="character", sep="\n")
	
	results_dir = addslash(results_dir)
	
	# Data to store
	initial_log_likelihood = NULL
	dispersal_rates = NULL
	extinction_rates = NULL
	final_log_likelihood = NULL
	#ancestral_splits_array = NULL
	#ancestral_states_array = NULL
	
	# Go through the lines
	#filecount = 0
	for (i in 1:length(tmplines))
		{
		tmpline = tmplines[i]
		
		# Start a new output file
		if (grepl("Global ML at root node:", tmpline) == TRUE)
			{
			filecount = filecount + 1
			cat("\nProcessing tree #", filecount, sep="")
			tmpi = sprintf("%05.0f", filecount)
	
			new_prefix = get_fn_prefix(fn=get_path_last(path=outfn))
			if (outputfiles == TRUE)
				{
				if (new_splits_fn == TRUE)
					{
					splits_table_fn = np(paste(addslash(results_dir), new_prefix, "_splits", tmpi, ".txt", sep=""))
					}
				if (new_states_fn == TRUE)
					{
					states_table_fn = np(paste(addslash(results_dir), new_prefix, "_states", tmpi, ".txt", sep=""))
					}
				}
	
			# you need to clear these files
			#new_splits_fn = TRUE
			#new_states_fn = TRUE
			}
		
		# parse initial ln likelihood
		if (grepl("initial -ln likelihood", tmpline) == TRUE)
			{
			words = strsplit_whitespace(tmpline)
			tmpnum = -1 * as.numeric(words[length(words)])
			initial_log_likelihood = c(initial_log_likelihood, tmpnum)
			}
	
		# parse final ln likelihood
		if (grepl("-lnL =", tmpline) == TRUE)
			{
			words = strsplit_whitespace(tmpline)
			tmpnum = -1 * as.numeric(words[length(words)])
			final_log_likelihood = c(final_log_likelihood, tmpnum)
			}
	
		# parse dispersal line
		if (grepl("dispersal = ", tmpline) == TRUE)
			{
			words = strsplit_whitespace(tmpline)
			tmpnum = as.numeric(words[3])
			dispersal_rates = c(dispersal_rates, tmpnum)
			}
		# parse extinction line
		if (grepl("extinction = ", tmpline) == TRUE)
			{
			words = strsplit_whitespace(tmpline)
			tmpnum = as.numeric(words[3])
			extinction_rates = c(extinction_rates, tmpnum)
			}
		
		# Ancestral splits
		if (grepl("At node N", tmpline) == TRUE)
			{
			# first line, e.g. "At node N36:"
			words = strsplit_whitespace(tmpline)
			tmpnum3 = words[length(words)]
			tmpnum2 = gsub(pattern="N", replacement="", x=tmpnum3)
			tmpnum2 = gsub(pattern="\\:", replacement="", x=tmpnum2)
			tmpnum1 = as.numeric(tmpnum2)

			# second line (there are only 15 values reported, plus header
			for (j in ((i+1):(i+20)))
				{
				# Test if the line exists
				if (is.na(tmplines[j]))
					{
					break()
					}

				# Test if the line is header:    split     lnL     Rel.Prob
				if ( (grepl("split", tmplines[j]) == TRUE) && (grepl("lnL", tmplines[j]) == TRUE) )
					{
					next()
					}

	
				#cat("\nstart", tmplines[j], "end", sep="")
				if (tmplines[j] == "" || grepl("At node N", tmplines[j]) == TRUE || grepl("initializing", tmplines[j]) == TRUE)
					{
					break()
					} else {
					words = strsplit_whitespace(tmplines[j])
					
					# Fix splits by removing []
					words[1] = gsub(pattern="\\[", replacement="", x=words[1])
					words[1] = gsub(pattern="\\]", replacement="", x=words[1])
					
					# -LnL
					words[2] = -1* as.numeric(extract_numbers(words[2]))
					
					# Rel. prob.
					words[3] = extract_numbers(words[3])
					outwords = c(tmpnum1, words)
					outstr = list2str(outwords, spacer="\t")
	
					#ancestral_splits_array = rbind(ancestral_splits_array, tmpline)
					# Append data to end of file
					# Write data to file, if outputfiles==TRUE
					if (outputfiles == TRUE)
						{
						if (new_splits_fn == TRUE)
							{
							write(outstr, file=splits_table_fn, append=append, sep="\n")
							#new_splits_fn = FALSE
							} else {
							pass=1
							}
						}
					}
				}
			}
	
		# Ancestral states
		if (grepl("Ancestral states for:", tmpline) == TRUE)
			{
			# first line
			words = strsplit_whitespace(tmpline)
			tmpnum1 = as.numeric(words[length(words)])
			
			# second line
			for (j in ((i+1):(i+20)))
				{
				# Test if the line exists
				if (is.na(tmplines[j]))
					{
					break()
					}
	
				#cat("\nstart", tmplines[j], "end", sep="")
				if (tmplines[j] == "" || grepl("Ancestral", tmplines[j]) == TRUE || grepl("initializing", tmplines[j]) == TRUE)
					{
					break
					} else {
					words = strsplit_whitespace(tmplines[j])
					words[3] = extract_numbers(words[3])
					outwords = c(tmpnum1, words)
					outstr = list2str(outwords, spacer="\t")
					
					#ancestral_states_array = rbind(ancestral_states_array, tmpline)
					# Append data to end of file
					# Write data to file, if outputfiles==TRUE
					if (outputfiles == TRUE)
						{
						if (new_states_fn == TRUE)
							{
							write(outstr, file=states_table_fn, append=append, sep="\n")
							#new_states_fn = FALSE
							} else {
							pass=1
							}
						}
					}
				}
			}
		#print(outstr)
		}
	cat("\n")

	# If no initial log-likelihood, but there are other values, put in NA
	if (!is.null(final_log_likelihood))
		{
		if (is.null(initial_log_likelihood))
			{
			initial_log_likelihood = NA
			}
		}


	# Write data to file, if outputfiles==TRUE
	if (outputfiles == TRUE)
		{

		if (new_splits_fn == TRUE)
			{
			cat("\nSplits output to: ", splits_table_fn, "\n", sep="")
			} else {
			splits_table_fn = ""
			}


		#cat("\nStates output (if LAGRANGE Python does this, which I don't think it does) to: ", splits_table_fn, "\n", sep="")
		if (new_states_fn == TRUE)
			{
			cat("\nStates output (if LAGRANGE Python does this, which I don't think it does) to: ", states_table_fn, "\n", sep="")
			} else {
			states_table_fn = ""
			}

		# Output to table
		sumstats = cbind(initial_log_likelihood, dispersal_rates, extinction_rates, final_log_likelihood, splits_table_fn, states_table_fn)
		sumstats = adf2(sumstats)


		# Output sumstats
		outfn = np(paste(results_dir, "/sumstats.txt", sep=""))
		write.table(x=sumstats, file=outfn, append=FALSE, quote=FALSE, sep="	", row.names=FALSE, col.names=TRUE)
		} else {
		# Output to table
		sumstats = cbind(initial_log_likelihood, dispersal_rates, extinction_rates, final_log_likelihood)
		sumstats = adf2(sumstats)
		}
	
	
	return(sumstats)
	}








##########################################
# LAGRANGE (C++ version) utilities
##########################################


#######################################################
# parse_lagrange_output
#######################################################
#' Parse the output file from python \code{LAGRANGE}
#' 
#' Parse the output of a C++ \code{LAGRANGE} run.  
#' 
#' This function parses the output of \code{LAGRANGE}, obtained by a command such as the following, run at a UNIX/Mac
#' Terminal command line.
#' 
#' \code{cd /Users/nick/Desktop/__projects/_2011-07-15_Hannah_spider_fossils/_data/lagrange_for_nick}
#' 
#' \code{./lagrange_cpp palp_no_Lacun_v1_2nd387.lg > lagrange_results_v1_2nd387.txt}
#' 
#' C++ LAGRANGE can be obtained at \url{https://code.google.com/p/lagrange/}
#' 
#' @param outfn The C++ \code{LAGRANGE} output text file.
#' @param outputfiles Should parsed output be written to files? Default FALSE.
#' @param results_dir The directory \code{outfn} is in.
#' @param new_splits_fn Should a text file containing a table of the splits and their probabilities be output? Default \code{FALSE}.
#' @param new_states_fn Should a text file containing a table of the states and their probabilities be output? Default \code{TRUE}, 
#' unlike python \code{LAGRANGE}, C++ \code{LAGRANGE} \emph{will} output the states at the nodes.
#' @param filecount The starting number for the filecount (relevant if one is processing many files).
#' @return sumstats A \code{\link[base]{data.frame}} containing the summary statistics (LnL, d and e rates, etc.)  The splits
#' filename is output to screen.
#' @export
#' @seealso \code{\link{get_lagrange_nodenums}}, \code{\link{LGpy_splits_fn_to_table}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @examples
#' test=1
#' 
parse_lagrange_output <- function(outfn, outputfiles=FALSE, results_dir=getwd(), new_splits_fn=FALSE, new_states_fn=TRUE, filecount=0)
	{
	setup='
	outfn = "/Users/nick/Desktop/__projects/_2011-07-15_Hannah_spider_fossils/_data/lagrange_for_nick/_run_3rd226/lagrange_results_v1_3rd226.txt"
	results_dir = "/Users/nick/Desktop/__projects/_2011-07-15_Hannah_spider_fossils/_data/lagrange_for_nick/_results_3rd226/"
	'


	# Split the file into 1 file for each analysis
	tmplines = scan(outfn, what="character", sep="\n")
	
	results_dir = addslash(results_dir)
	
	# Data to store
	initial_log_likelihood = NULL
	dispersal_rates = NULL
	extinction_rates = NULL
	final_log_likelihood = NULL
	#ancestral_splits_array = NULL
	#ancestral_states_array = NULL
	
	# Go through the lines
	#filecount = 0
	for (i in 1:length(tmplines))
		{
		tmpline = tmplines[i]
		
		# Start a new output file
		if (grepl("starting likelihood calculations", tmpline) == TRUE)
			{
			filecount = filecount + 1
			cat("\nProcessing tree #", filecount, sep="")
			tmpi = sprintf("%05.0f", filecount)
	
			new_prefix = get_fn_prefix(fn=get_path_last(path=outfn))
			if (outputfiles == TRUE)
				{
				if (new_splits_fn == TRUE)
					{
					splits_table_fn = np(paste(addslash(results_dir), new_prefix, "_splits", tmpi, ".txt", sep=""))
					}
				if (new_states_fn == TRUE)
					{
					states_table_fn = np(paste(addslash(results_dir), new_prefix, "_states", tmpi, ".txt", sep=""))
					}
				}
	
			# you need to clear these files
			#new_splits_fn = TRUE
			#new_states_fn = TRUE
			}
		
		# parse initial ln likelihood
		if (grepl("initial -ln likelihood", tmpline) == TRUE)
			{
			words = strsplit_whitespace(tmpline)
			tmpnum = -1 * as.numeric(words[length(words)])
			initial_log_likelihood = c(initial_log_likelihood, tmpnum)
			}
	
		# parse final ln likelihood
		if (grepl("final -ln likelihood", tmpline) == TRUE)
			{
			words = strsplit_whitespace(tmpline)
			tmpnum = -1 * as.numeric(words[length(words)])
			final_log_likelihood = c(final_log_likelihood, tmpnum)
			}
	
		# parse dispersal/extinction lines
		if (grepl("dis: ", tmpline) == TRUE)
			{
			words = strsplit_whitespace(tmpline)
			tmpnum = as.numeric(words[2])
			dispersal_rates = c(dispersal_rates, tmpnum)
	
			tmpnum = as.numeric(words[4])
			extinction_rates = c(extinction_rates, tmpnum)
			}
		
		# Ancestral splits
		if (grepl("Ancestral splits for:", tmpline) == TRUE)
			{
			# first line
			words = strsplit_whitespace(tmpline)
			tmpnum1 = as.numeric(words[length(words)])
			
			# second line
			for (j in ((i+1):(i+20)))
				{
				# Test if the line exists
				if (is.na(tmplines[j]))
					{
					break()
					}
	
				#cat("\nstart", tmplines[j], "end", sep="")
				if (tmplines[j] == "" || grepl("Ancestral", tmplines[j]) == TRUE || grepl("initializing", tmplines[j]) == TRUE)
					{
					break
					} else {
					words = strsplit_whitespace(tmplines[j])
					words[3] = extract_numbers(words[3])
					outwords = c(tmpnum1, words)
					outstr = list2str(outwords, spacer="\t")
	
					#ancestral_splits_array = rbind(ancestral_splits_array, tmpline)
					# Write data to file, if outputfiles==TRUE
					if (outputfiles == TRUE)
						{
						if (new_splits_fn == TRUE)
							{
							write(outstr, file=splits_table_fn, append=FALSE, sep="\n")
							new_splits_fn = FALSE
							} else {
							write(outstr, file=splits_table_fn, append=TRUE, sep="\n")
							}
						}
					}
				}
			}
	
		# Ancestral states
		if (grepl("Ancestral states for:", tmpline) == TRUE)
			{
			# first line
			words = strsplit_whitespace(tmpline)
			tmpnum1 = as.numeric(words[length(words)])
			
			# second line
			for (j in ((i+1):(i+20)))
				{
				# Test if the line exists
				if (is.na(tmplines[j]))
					{
					break()
					}
	
				#cat("\nstart", tmplines[j], "end", sep="")
				if (tmplines[j] == "" || grepl("Ancestral", tmplines[j]) == TRUE || grepl("initializing", tmplines[j]) == TRUE)
					{
					break
					} else {
					words = strsplit_whitespace(tmplines[j])
					words[3] = extract_numbers(words[3])
					outwords = c(tmpnum1, words)
					outstr = list2str(outwords, spacer="\t")
					
					#ancestral_states_array = rbind(ancestral_states_array, tmpline)
					# Append data to end of file
					# Write data to file, if outputfiles==TRUE
					if (outputfiles == TRUE)
						{
						if (new_states_fn == TRUE)
							{
							write(outstr, file=states_table_fn, append=FALSE, sep="\n")
							new_states_fn = FALSE
							} else {
							write(outstr, file=states_table_fn, append=TRUE, sep="\n")
							}
						}
					}
				}
			}
		#print(outstr)
		}
	cat("\n")



	# Write data to file, if outputfiles==TRUE
	if (outputfiles == TRUE)
		{
		# Make an output table
		sumstats = cbind(initial_log_likelihood, dispersal_rates, extinction_rates, final_log_likelihood, splits_table_fn, states_table_fn)
		sumstats = adf2(sumstats)

		cat("\nStates output to: ", splits_table_fn, "\n", sep="")
		cat("\nSplits output to: ", splits_table_fn, "\nn", sep="")
		
		outfn = np(paste(results_dir, "sumstats.txt", sep=""))
		write.table(x=sumstats, file=outfn, append=FALSE, quote=FALSE, sep="	", row.names=FALSE, col.names=TRUE)
		# write_table_good(sumstats, paste(results_dir, "sumstats.txt", sep=""))
		} else {
		# Output to table
		sumstats = cbind(initial_log_likelihood, dispersal_rates, extinction_rates, final_log_likelihood)
		sumstats = adf2(sumstats)
		}

	return(sumstats)
	}





