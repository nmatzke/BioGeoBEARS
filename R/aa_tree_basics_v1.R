

#' Get the node numbers of all tips descending from node 
#'
#' This function starts from the input node \code{nodenum},
#' and searches up the tree recursively for tip nodes.
#' It uses function \code{get_daughter_nodes} for the recursion.
#'
#' The recursion should make it faster than some other possible
#' methods that have to work over the whole tree.
#'
#' @param nodenum The node number to find tip descendants of.
#' @param tr An ape \code{phylo} object.
#' @return nodes The node numbers of the tips
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
#' # Load hard-coded tree
#' tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
#' get_daughter_tipnums(nodenum=1, tr)
#' get_daughter_tipnums(nodenum=2, tr)
#' get_daughter_tipnums(nodenum=3, tr)
#' get_daughter_tipnums(nodenum=4, tr)
#' get_daughter_tipnums(nodenum=5, tr)
#' 
get_daughter_tipnums <- function(nodenum, tr)
	{
	ex='
	get_daughter_tipnums(nodenum=1, tr)
	get_daughter_tipnums(nodenum=2, tr)
	get_daughter_tipnums(nodenum=3, tr)
	get_daughter_tipnums(nodenum=4, tr)
	get_daughter_tipnums(nodenum=5, tr)
	'
	nodes = get_daughter_nodes(nodenum, tr, nodes=NULL)
	tips_TF = nodes <= length(tr$tip.label)
	tipnums = nodes[tips_TF]
	return(tipnums)
	}


#' Get the tipnames of all tips descending from node 
#'
#' This function starts from the input node \code{nodenum},
#' and searches up the tree recursively for tip nodes (using
#' \code{get_daughter_tipnums}).  It uses function 
#' \code{get_daughter_nodes} for the recursion. Once tip
#' node numbers are found, the tip names are returned.
#'
#' The recursion should make it faster than some other possible
#' methods that have to work over the whole tree.
#'
#' @param nodenum The node number to find tip descendants of.
#' @param tr An ape \code{phylo} object.
#' @return nodtipnameses The names of the tips
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
#' # Load hard-coded tree
#' tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
#' get_daughter_tipnums(nodenum=1, tr)
#' get_daughter_tipnums(nodenum=2, tr)
#' get_daughter_tipnums(nodenum=3, tr)
#' get_daughter_tipnums(nodenum=4, tr)
#' get_daughter_tipnums(nodenum=5, tr)
#' 
get_daughter_tipnames <- function(nodenum, tr)
	{
	tipnums = get_daughter_tipnums(nodenum, tr)
	tipnames = tr$tip.label[tipnums]
	return(tipnames)
	}


#' Recursively get the nodes descending from a node 
#'
#' This function starts from the input node \code{nodenum},
#' and searches up the tree recursively for descendant nodes.
#' \code{get_daughter_nodes} is recursive, adding to the 
#' input/output list \code{nodes}.
#'
#' The recursion should make it faster than some other possible
#' methods that have to work over the whole tree.
#'
#' @param nodenum The node number to find tip descendants of.
#' @param tr An ape \code{phylo} object.
#' @param nodes The input/output list of nodes.
#' @return nodes The input/output list of nodes.
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
#' # Load hard-coded tree
#' tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
#' get_daughter_nodes(nodenum=1, tr, nodes=NULL)
#' get_daughter_nodes(nodenum=2, tr, nodes=NULL)
#' get_daughter_nodes(nodenum=3, tr, nodes=NULL)
#' get_daughter_nodes(nodenum=4, tr, nodes=NULL)
#' get_daughter_nodes(nodenum=5, tr, nodes=NULL)
#' get_daughter_nodes(nodenum=5, tr, nodes=c(0))
#' 
get_daughter_nodes <- function(nodenum, tr, nodes=NULL)
	{
	ex='
	get_daughter_nodes(nodenum=1, tr, nodes=NULL)
	get_daughter_nodes(nodenum=2, tr, nodes=NULL)
	get_daughter_nodes(nodenum=3, tr, nodes=NULL)
	get_daughter_nodes(nodenum=4, tr, nodes=NULL)
	get_daughter_nodes(nodenum=5, tr, nodes=NULL)
	get_daughter_nodes(nodenum=5, tr, nodes=c(0))
	'
	
	if(is.null(nodes))
		{
		nodes = vector()
		}
	daughter_nodes = tr$edge[which(tr$edge[,1]==nodenum),2]
	
	# Error check, in case the starting nodenum is a tip
	if ((length(daughter_nodes) == 0) && (length(nodes)==0))
		{
		nodes = c(nodes, nodenum)
		} else {
		nodes = c(nodes, daughter_nodes)
		}
		
	internal_nodes_indices = which(daughter_nodes > length(tr$tip.label))
	if(length(internal_nodes_indices) > 0)
		{
		for (i in 1:length(internal_nodes_indices))
			{
			nodes = get_daughter_nodes(nodenum=daughter_nodes[internal_nodes_indices[i]], tr=tr, nodes=nodes)
			}
		}
	return(nodes)
	}



#' Trace from a node up to its parents etc. a specified distance
#'
#' This function starts from the input node \code{nodenum},
#' and searches "down" the tree (rootwards) for parental nodes.
#' \code{trace_parents_up} is recursive, continuing to 
#' find parents-of-parent nodes until depthtime is reached,
#' at which point the node "owning" that branch is identified.
#' 
#' Note: in R, every node "owns" a branch below it (rootwards to it), except
#' (usually) the root node. Knowing the identity of this 
#' node is useful, e.g. if you want to find the appropriate 
#' location to insert a direct ancestor of a particular age.
#'
#' The recursion should make it faster than some other possible
#' methods that have to work over the whole tree.
#'
#' @param nodenum The node number to find tip descendants of.
#' @param t An ape \code{phylo} object.
#' @param depthtime The distance below the node \code{nodenum} at which 
#'        to stop search.
#' @return nodenum The parental node "owning" the branch where \code{depthtime} was reached.
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
#' # Load hard-coded tree
#' tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
#' trace_parents_up(nodenum=1, tr, depthtime=0.1)
#' trace_parents_up(nodenum=1, tr, depthtime=1.1)
#' trace_parents_up(nodenum=5, tr, depthtime=0.5)
#' 
#' # These return errors, the error messages explain the reason:
#' \dontrun{
#' trace_parents_up(nodenum=1, tr, depthtime=2.1)
#' trace_parents_up(nodenum=1, tr, depthtime=1.0)
#' trace_parents_up(nodenum=4, tr, depthtime=0.5)
#' }
#' 
trace_parents_up <- function(nodenum, t, depthtime)
	{
	ex='
	tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
	trace_parents_up(nodenum=1, tr, depthtime=0.1)
	trace_parents_up(nodenum=1, tr, depthtime=1.1)
	trace_parents_up(nodenum=1, tr, depthtime=2.1)
	trace_parents_up(nodenum=1, tr, depthtime=1.0)
	trace_parents_up(nodenum=4, tr, depthtime=0.5)
	trace_parents_up(nodenum=5, tr, depthtime=0.5)
	'
	# Trace from a node up to its parents etc., a specified distance
	parent_node = get_parent_for_trace_parents_up(nodenum, t)
	
	# print nodenum
	#print(nodenum)
	#print(parent_node)
	#print(depthtime)
	
	length_to_parent = t$edge.length[t$edge[,2] == nodenum]
	#cat("length_to_parent: ", length_to_parent, ", length=", length(length_to_parent), "\n", sep="")
	#cat("depthtime: ", depthtime, "\n", sep="")
	if (length(length_to_parent) == 0)
		{
		print("ERROR: trace_parents_up() -- no length_to_parent returned, probably overshot bottom of tree")
		return(NA)
		}
	if (length_to_parent == depthtime)
		{
		print("ERROR: trace_parents_up() doesn't want to find an EXACT match between depthtime and a node; this will lead to problems in hook addition!")
		}
	if (length_to_parent > depthtime)
		{
		# you're done!
		return(nodenum)
		}
	else
		{
		# burrow up to parents
		depthtime = depthtime - length_to_parent
		parent_node = trace_parents_up(parent_node, t, depthtime)
		return(parent_node)
		}
	return(nodenum)
	}


#' Get first parent for trace_parents_up
#'
#' Utility for use with \code{trace_parents_up}.
#'
#' @param nodenum The node number to find tip descendants of.
#' @param t An ape \code{phylo} object.
#' @param printflag If TRUE, print messages about success/failure of
#'                  finding the parent node. Default FALSE.
#' @return parent_nodenum The node directly parental to \code{nodenum}.
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
#' # Load hard-coded tree
#' tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
#' get_parent_for_trace_parents_up(nodenum=1, tr, printflag=FALSE)
#' get_parent_for_trace_parents_up(nodenum=1, tr, printflag=FALSE)
#' get_parent_for_trace_parents_up(nodenum=5, tr, printflag=FALSE)
#' get_parent_for_trace_parents_up(nodenum=1, tr, printflag=TRUE)
#' get_parent_for_trace_parents_up(nodenum=1, tr, printflag=TRUE)
#' get_parent_for_trace_parents_up(nodenum=4, tr, printflag=FALSE)
#' 
#' # This returns an error instead of NA:
#' \dontrun{
#' get_parent_for_trace_parents_up(nodenum=4, tr, printflag=TRUE)
#' }
#' 
get_parent_for_trace_parents_up <- function(nodenum, t, printflag=FALSE)
	{
	ex='
	# Load hard-coded tree
	tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
	get_parent_for_trace_parents_up(nodenum=1, tr, printflag=FALSE)
	get_parent_for_trace_parents_up(nodenum=1, tr, printflag=FALSE)
	get_parent_for_trace_parents_up(nodenum=5, tr, printflag=FALSE)
	get_parent_for_trace_parents_up(nodenum=1, tr, printflag=TRUE)
	get_parent_for_trace_parents_up(nodenum=1, tr, printflag=TRUE)
	get_parent_for_trace_parents_up(nodenum=4, tr, printflag=FALSE)

	# This returns an error instead of NA:
	get_parent_for_trace_parents_up(nodenum=4, tr, printflag=TRUE)
	'
	matching_edges = findall(nodenum, t$edge[,2])
	parent_nodenum = t$edge[,1][matching_edges][1]
	
	if (printflag)
		{
		print(paste("nodenum=", nodenum, " parent_nodenum=", parent_nodenum, sep=""))
		}
	if (is.na(parent_nodenum))
		{
		if (printflag)
			{
			print(paste("get_parent(): node ", nodenum, " has no parent, it's probably the root!\nAnd you missed whatever parent you were actually trying to find!", sep=""))
			}
		}
	return(parent_nodenum)
	}
	


#' Recursive algorithm to get distance between descendent and ancestor
#'
#' The problem with \code{ape}'s \code{\link{dist.nodes}} is that it
#' returns of matrix of size numnodes x numnodes. This gets very
#' large and efficient for large trees, especially if you just
#' want the distance between 2 particular nodes.
#' 
#' Here, \code{dist_between_direct_ancestors} uses
#' \code{\link{get_parent}} to recursively trace between
#' a descendant and ancestor node, calculating the total
#' distance (total branchlength).
#'
#' @param ancestor_node
#' @param descendant_node
#' @param t An ape \code{phylo} object.
#' @param totaldist The total distance (default starts at zero, as this is added
#'                  to recursively, as branches are iterated through).
#' @param printflag Prints warning if \code{ancestor_node == descendant_node}, then 
#"                  returns \code{totaldist} (usually 0). Default \code{FALSE}.
#' @return totaldist The total distance between the nodes
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
#' tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
#' dist_between_direct_ancestors(ancestor_node=4, descendant_node=1, t=tr, totaldist=0, printflag=FALSE)
#' dist_between_direct_ancestors(ancestor_node=5, descendant_node=1, t=tr, totaldist=0, printflag=FALSE)
#' dist_between_direct_ancestors(ancestor_node=5, descendant_node=2, t=tr, totaldist=0, printflag=FALSE)
#' dist_between_direct_ancestors(ancestor_node=4, descendant_node=3, t=tr, totaldist=0, printflag=FALSE)
#' dist_between_direct_ancestors(ancestor_node=3, descendant_node=3, t=tr, totaldist=0, printflag=FALSE)
#' 
#' \dontrun{
#' # Returns warning
#' dist_between_direct_ancestors(ancestor_node=3, descendant_node=3, t=tr, totaldist=0, printflag=TRUE)
#' 
#' # Returns error (as node 5 is not ancestral to tipnode 3)
#' dist_between_direct_ancestors(ancestor_node=5, descendant_node=3, t=tr, totaldist=0, printflag=FALSE)
#' }
dist_between_direct_ancestors <- function(ancestor_node, descendant_node, t, totaldist=0, printflag=FALSE)
	{
	ex='
	tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
	dist_between_direct_ancestors(ancestor_node=4, descendant_node=1, t=tr, totaldist=0, printflag=FALSE)
	dist_between_direct_ancestors(ancestor_node=5, descendant_node=1, t=tr, totaldist=0, printflag=FALSE)
	dist_between_direct_ancestors(ancestor_node=5, descendant_node=2, t=tr, totaldist=0, printflag=FALSE)
	dist_between_direct_ancestors(ancestor_node=4, descendant_node=3, t=tr, totaldist=0, printflag=FALSE)
	dist_between_direct_ancestors(ancestor_node=3, descendant_node=3, t=tr, totaldist=0, printflag=FALSE)
	
	# Returns warning
	dist_between_direct_ancestors(ancestor_node=3, descendant_node=3, t=tr, totaldist=0, printflag=TRUE)
	
	# Returns error (as node 5 is not ancestral to tipnode 3)
	dist_between_direct_ancestors(ancestor_node=5, descendant_node=3, t=tr, totaldist=0, printflag=FALSE)
	'


	# Recursive algorithm to get distance between descendent and ancestor
	
	# Error trap should operate before this
	
	if (ancestor_node == descendant_node)
		{
		if (printflag)
			{
			txt = "dist_between_direct_ancestors(): ancestor_node == descendant_node"
			warning(txt)
			print(txt)
			}
		return(totaldist)
		}
	
	parent_node = get_parent(descendant_node, t)
	dist_to_parent = t$edge.length[t$edge[,2] == descendant_node]
	
	totaldist = dist_to_parent + totaldist
	
	#print(paste(parent_node, ancestor_node, sep=""))
	if (parent_node == ancestor_node)
		{
		return(totaldist)
		}
	else
		{
		totaldist = dist_between_direct_ancestors(ancestor_node, parent_node, t, totaldist)
		return(totaldist)
		}
	}



#' Add a hook (a small side branch) to a phylogeny
#'
#' Given a tree, tipname, depth below the tip, and new tip name, adds the 
#' new tip as a "hook" (a very short side-branch).
#'
#' In \code{BioGeoBEARS}, direct ancestors are coded as very-very short side 
#' branches ("hooks") in the tree. This is an easy way to include direct ancestors, and 
#' save/load the resulting trees in standard formats like Newick.
#' 
#" It does, though, require that the program in question "knows" to read short 
#' branches this way. \code{BioGeoBEARS} does this, but other programs might not.
#'
#' @param t An ape \code{phylo} object.
#' @param tipname The name of the tip below which the hook will be added (the hook can
#'                be added many nodes below the specified tipname; the total distance 
#'                below the tip, "\code{depthtime}", is what determines this).
#' @param depthtime The length of time (total branchlength) below tip \code{tipname} 
#'                  where the hook will be added, with new direct ancestor tip 
#'                  \code{newtipname}.
#' @param plottree Should the updated tree be plotted? Default \code{FALSE}.
#' @param printflag If 2 or higher, print all messages about success/failure of
#'                  steps of function. 1 prints less in some cases. Default 0.
#' @param newtipname Name of new direct ancestor tip added via a hook. If 
#' \code{newtipname="default"}, new tips get "hook" and then the node number, 
#' e.g. hook414. In \code{add_hooks} (but not directly in \code{add_hook}), 
#' if \code{newtipnames="tipnames"}, new tips get either 
#' (1) tipname_age (for hooks on tip branches), or (2) hook_tipname_age (for hooks 
#' on internal branches; just one of the tips is listed)
#' @return newtree2 The revised tree, with the hook added. Default name is (cleverly) "default".
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
#' tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
#' add_hook(t=tr, tipname="gorilla", depthtime=0.5, brlen_of_side_branch=0.0000001, plottree=FALSE, printflag=0, newtipname="default")
#' add_hook(t=tr, tipname="gorilla", depthtime=1.5, brlen_of_side_branch=0.0000001, plottree=FALSE, printflag=0, newtipname="default")
#' add_hook(t=tr, tipname="chimp", depthtime=1.5, brlen_of_side_branch=0.0000001, plottree=FALSE, printflag=0, newtipname="default")
#' add_hook(t=tr, tipname="human", depthtime=1.5, brlen_of_side_branch=0.0000001, plottree=FALSE, printflag=0, newtipname="default")
#' add_hook(t=tr, tipname="chimp", depthtime=1.5, brlen_of_side_branch=0.0000001, plottree=FALSE, printflag=0, newtipname="fossil_tip1")
#' add_hook(t=tr, tipname="human", depthtime=1.5, brlen_of_side_branch=0.0000001, plottree=FALSE, printflag=0, newtipname="fossil_tip2")
#' 
#' # Produces error
#' add_hook(t=tr, tipname="human", depthtime=2.5, brlen_of_side_branch=0.0000001, plottree=FALSE, printflag=0, newtipname="default")
#' 
add_hook <- function(t, tipname, depthtime, brlen_of_side_branch=0.0000001, plottree=FALSE, printflag=0, newtipname="default")
	{
	ex='
	tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
	add_hook(t=tr, tipname="gorilla", depthtime=0.5, brlen_of_side_branch=0.0000001, plottree=FALSE, printflag=0, newtipname="default")
	add_hook(t=tr, tipname="gorilla", depthtime=1.5, brlen_of_side_branch=0.0000001, plottree=FALSE, printflag=0, newtipname="default")
	add_hook(t=tr, tipname="chimp", depthtime=1.5, brlen_of_side_branch=0.0000001, plottree=FALSE, printflag=0, newtipname="default")
	add_hook(t=tr, tipname="human", depthtime=1.5, brlen_of_side_branch=0.0000001, plottree=FALSE, printflag=0, newtipname="default")
	add_hook(t=tr, tipname="chimp", depthtime=1.5, brlen_of_side_branch=0.0000001, plottree=FALSE, printflag=0, newtipname="fossil_tip1")
	add_hook(t=tr, tipname="human", depthtime=1.5, brlen_of_side_branch=0.0000001, plottree=FALSE, printflag=0, newtipname="fossil_tip2")
	
	# Produces error
	add_hook(t=tr, tipname="human", depthtime=2.5, brlen_of_side_branch=0.0000001, plottree=FALSE, printflag=0, newtipname="default")
	
	'
	
	# Add a hook (a small side tip) to a phylogeny
	#
	# e.g.:
	# Do spatial variogram by doing points from many different species
	# add tips to tree
	#cat("owls(((Strix_aluco:4.2,Asio_otus:4.2):3.1,Athene_noctua:7.3):6.3,Tyto_alba:13.5);", file = "ex.tre", sep = "\n")
	#t <- read.tree("ex.tre")
	#prt(t)
	#newtree = add_hook(t, brlen_of_side_branch=0.0000001, plottree = TRUE)

	if (printflag >= 2)
		{
		cat("add_hook(): Adding below tipname ",  tipname, "\n", sep="")
		}

	newtree = t
	#daughter_nodenum_to_add_hook_below = 4
	#height_below_daughter_at_which_to_add_it = 1.0
	#brlen_of_side_branch=0.0000001
	
	# find the node, below this you will add the 
	tip_nodenum = which(t$tip.label == tipname)
	
	if (printflag >= 2)
		{
		print(paste("addhook(): tip_nodenum = ", tip_nodenum, sep=""))
		}

	daughter_nodenum_to_add_hook_below = trace_parents_up(tip_nodenum, t, depthtime)
	
	# Error trap in case of e.g. overrun of bottom of tree
	if (is.na(daughter_nodenum_to_add_hook_below))
		{
		txt = "add_hook() error: daughter_nodenum_to_add_hook_below is NA, probably overshot bottom of tree (?); not adding hook. Returning input tree.\n"
		cat(txt)
		warning(txt)
		return(t)
		}
	
	if (printflag >= 2)
		{
		print(paste("addhook(): internal nodenum to add below = ", daughter_nodenum_to_add_hook_below, sep=""))
		}
	
	tip_to_ancestor_node_dist = dist_between_direct_ancestors(daughter_nodenum_to_add_hook_below, tip_nodenum, t, totaldist=0)
	
	
	
	height_below_daughter_at_which_to_add_it = depthtime - tip_to_ancestor_node_dist
	
	
	# add a new tip to the list of tips (this is the hook)
	new_tip_nodenum = get_nodenum_structural_root(t)
	
	# bump up all of the nodenums above the node tip by 1
	newtree$edge[t$edge >= new_tip_nodenum] = t$edge[t$edge >= new_tip_nodenum] + 1
	
	# add a new internal node at the end
	new_inNode = max(newtree$edge) + 1
	
	# add two new edges, and replace the old edge
	#print(t$edge[,2])
	#print(daughter_nodenum_to_add_hook_below)
	#print(t$edge[,2] == daughter_nodenum_to_add_hook_below)
	old_edge_num = which(t$edge[,2] == daughter_nodenum_to_add_hook_below)
	
	# extract the edgenums before and after this insertion (exceptions for in case
	# if the modified row is the first or last row)
	if (old_edge_num == 1)
		{
		first_old_edges_rownums = NULL
		} else {
		first_old_edges_rownums = 1:(old_edge_num-1) #newtree$edge[1:(old_edge_num-1), ]
		}
	if (old_edge_num == nrow(t$edge))
		{
		second_old_edges_rownums = NULL
		} else {
		second_old_edges_rownums = (old_edge_num+1):nrow(t$edge) # newtree$edge[, ]
		}
	
	
	
	# replace the edge, keeping the old parent (which may have increased by 1! use newtree!!), put the new internal node as the daughter)
	replacement_edge_row = newtree$edge[old_edge_num, ]
	replacement_edge_row[2] = new_inNode
	
	# subtract the distance below the daughter, from the top
	replacement_edge_length = t$edge.length[old_edge_num] - height_below_daughter_at_which_to_add_it
	
	
	# make the new edge, which goes below the old daughter node
	# you have to bump the daughter_nodenum_to_add_hook_below if it is
	# >= to the new_tip_nodenum
	if (daughter_nodenum_to_add_hook_below >= new_tip_nodenum)
		{
		daughter_nodenum_to_add_hook_below = daughter_nodenum_to_add_hook_below + 1
		}
	new_edge_below_old_daughter_node = c(new_inNode, daughter_nodenum_to_add_hook_below)
	new_edge_below_old_daughter_node_edge_length = height_below_daughter_at_which_to_add_it
	
	# make the new edge, which goes below the new tip: c(parent, daughter)
	new_edge_below_new_tip = c(new_inNode, new_tip_nodenum)
	new_edge_below_new_tip_edge_length = brlen_of_side_branch
	
	
	# add the edge rows before the one that is replaced, then the replaced edge, then the other old edges, then the 2 new edges
	new_edge_table = rbind(newtree$edge[first_old_edges_rownums, ], replacement_edge_row, newtree$edge[second_old_edges_rownums, ], new_edge_below_old_daughter_node, new_edge_below_new_tip)
	
	new_edgelength_list = c(t$edge.length[first_old_edges_rownums], replacement_edge_length, t$edge.length[second_old_edges_rownums], new_edge_below_old_daughter_node_edge_length, new_edge_below_new_tip_edge_length)
	
	# it MAY be important that the node numbers be INTEGER, not NUMERIC
	newtree$edge = matrix(as.integer(new_edge_table), ncol=2, byrow=FALSE)
	
	#row.names(newtree$edge) = NULL
	#newtree$edge[,1] = as.integer(newtree$edge[,1])
	#newtree$edge[,2] = as.integer(newtree$edge[,2])
	
	newtree$edge.length = new_edgelength_list
	#row.names(newtree$edge.length) = NULL
	
	# update number of internal nodes
	newtree$Nnode = t$Nnode + 1
	
	# add the new tip to the end of the list of tips
	if (newtipname == "default")
		{
		newtipname = paste("hook", new_tip_nodenum, sep="")
		} else {
		newtipname = newtipname
		}
	newtree$tip.label = c(t$tip.label, newtipname)
	
	if (printflag >= 2)
		{
		cat("\nAdding tip: ",  newtipname, sep="")
		}	
	
	# some crap to fix the tree formatting somehow
	# I mean, really, the tree was friggin logically correct, but
	# hanging plot and dist.nodes and probably anything
	# using reorder(phy, "pruningwise"), but I couldn't figure out why
	# I guess the order of the tips is important for some reason?
	# like maybe leftmost tip is leftmost in branching?
	# wtf kind of data architecture is this?
	# anyway, screw it, writing to Newick and reading back fixes it.
	newtree = reorder(newtree)
	#tmpfn = "tmp_junktree.tree"
	#write.tree(newtree, tmpfn)
	newtree2 = read.tree(file="", text=write.tree(newtree, file=""))

	
	#cat("\nDone adding tips.\n")
	
	# plot, if desired:
	if (plottree == TRUE)
		{
		cat("add_hook(): plotting/printing the resulting tree...\n", sep="")
		prt(newtree)
		plot(newtree)
		}
	
	return(newtree2)
	} # END add_hook <- function(t, tipname, depthtime, brlen_of_side_branch=0.0000001, plottree = FALSE, printflag=0, newtipname="default")





################################################################################
# TREE MODIFICATION FUNCTIONS (e.g. adding hooks, choosing certain branches)
################################################################################
# newtipname="default" means new tips get "hook" and then the node number, e.g. hook414
# newtipname="tipnames" means new tips get either 
#    tipname_age      (for hooks on tip branches), or 
#    hook_tipname_age (for hooks on internal branches; just one of the tips is listed)


#' Add hooks (small side branches) to a phylogeny
#'
#' Given a tree, tipname, depth below the tip, and new tip name, adds the 
#' new tip as a "hook" (a very short side-branch).
#'
#' In \code{BioGeoBEARS}, direct ancestors are coded as very-very short side 
#' branches ("hooks") in the tree. This is an easy way to include direct ancestors, and 
#' save/load the resulting trees in standard formats like Newick.
#' 
#" It does, though, require that the program in question "knows" to read short 
#' branches this way. \code{BioGeoBEARS} does this, but other programs might not.
#'
#' @param tr An ape \code{phylo} object.
#' @param list_of_times_before_present A vector of the lengths of time before the present
#'                (highest tip). Hooks will be added at each point where a branch 
#'                crosses the timepoint.
#' @param brlen_of_side_branch The branch length of the hook branch. Default 0.0000001.
#' @param plottree Should the updated tree be plotted? Default \code{FALSE}.
#' @param printflag If 2 or higher, print all messages about success/failure of
#'                  steps of function. 1 prints less in some cases. Default 0.
#' @param newtipnames If \code{newtipnames="default"}, new tips get "hook" and then the node number, 
#' e.g. hook414. If \code{newtipnames="tipnames"}, new tips get either 
#' (1) tipname_age (for hooks on tip branches), or (2) hook_tipname_age (for hooks 
#' on internal branches; just one of the tips is listed)
#' @return newtree2 The revised tree, with the hook added. Default name is (cleverly) "default".
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
#' tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
#' tr2 = add_hooks(tr, list_of_times_before_present=c(0.5,1.5), brlen_of_side_branch=0.0000001, plottree=FALSE, printflag=0, newtipnames="default")
#' write.tree(tr2, file="")
#' tr3 = add_hooks(tr, list_of_times_before_present=c(1.5), brlen_of_side_branch=0.0000001, plottree=FALSE, printflag=0, newtipnames="default")
#' write.tree(tr3, file="")
#' 
add_hooks <- function(tr, list_of_times_before_present, brlen_of_side_branch=0.0000001, plottree=FALSE, printflag=0, newtipnames="default")
	{
	ex='
	tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
	tr2 = add_hooks(tr, list_of_times_before_present=c(0.5,1.5), brlen_of_side_branch=0.0000001, plottree=FALSE, printflag=0, newtipnames="default")
	write.tree(tr2, file="")
	tr3 = add_hooks(tr, list_of_times_before_present=c(0.4,0.5,1.5), brlen_of_side_branch=0.0000001, plottree=FALSE, printflag=0, newtipnames="default")
	write.tree(tr3, file="")
	tr4 = add_hooks(tr, list_of_times_before_present=c(1.5), brlen_of_side_branch=0.0000001, plottree=FALSE, printflag=0, newtipnames="default")
	write.tree(tr4, file="")
	'
	# Take a list of ages, add hooks to any branch existing at that age

	# OK, RE-FREAKING DO
	# Gather a list of tip labels, and then record the distance below those tips that
	# you would go down to attach a hook
	hooktree = tr
	list_of_daughter_tipnames_to_add_hooks_below = c()
	list_of_ages_below_daughter = c()   # assumes ultrametric
	list_of_ages_absolute = c()
	list_of_new_tipnames = c()
	ntips = length(hooktree$tip.label)
	
	# Gather the list of tip labels, distance below each one, and new tip labels
	for (i in 1:length(list_of_times_before_present))
		{
		
		# Get the edges that exist at the time_slice in question
		time_slice = as.numeric(list_of_times_before_present[i])
		edge_times_bp = get_edge_times_before_present(hooktree)
		edges_that_exist_in_the_right_time = edges_existing_at_correct_time_bp_TF(time_slice, edge_times_bp)
		
		# get the nodes daughter to the branches that match
		nodenums_to_add_hooks_to = hooktree$edge[,2][edges_that_exist_in_the_right_time]

		if (printflag >= 1.5)
			{
			txt = paste0("Adding ", length(nodenums_to_add_hooks_to), " hooks for time #", i, "/", length(list_of_times_before_present), ": ", time_slice, " m.y.a.")
			cat("\n")
			cat(txt)
			}


		# calculate the times parent to these daughters at which to insert the hooks		
		#times_before_daughter_nodes = time_slice - edge_times_bp[edges_that_exist_in_the_right_time, 2]
		
		# trace these nodes to their tips in the (UNALTERED ORIGINAL) tree
# 		if (printflag >= 1.5)
# 			{
# 			cat("\n\tj(nodenum): ")
# 			}
		for (j in 1:length(nodenums_to_add_hooks_to))
			{
			nodenum = nodenums_to_add_hooks_to[j]

			if (printflag >= 2)
				{
				txt = paste0(j, ":", nodenum, " ")
				cat(txt)
				}

			
			# If the node is a tip
			if (nodenum <= ntips)
				{
				daughter_tipname_to_add_hooks_below = hooktree$tip.label[nodenum]
				
				edgenums = 1:nrow(edge_times_bp)
				edgenum = edgenums[tr$edge[,2] == nodenum]
				node_time_of_daughter = edge_times_bp[edgenum, 2]
				
				time_before_daughter_nodes = time_slice - node_time_of_daughter
				
				list_of_daughter_tipnames_to_add_hooks_below = c(list_of_daughter_tipnames_to_add_hooks_below, daughter_tipname_to_add_hooks_below)
				list_of_ages_below_daughter = c(list_of_ages_below_daughter, time_before_daughter_nodes)
				list_of_ages_absolute = c(list_of_ages_absolute, time_slice)
				list_of_new_tipnames = c(list_of_new_tipnames, daughter_tipname_to_add_hooks_below)
				} else {
				# Get *a* tip to compare time slice to
				temp_tipnames = get_all_daughter_tips_of_a_node(nodenum, hooktree)
				temp_tipname = temp_tipnames[1]
				namenums = 1:length(hooktree$tip.label)
				temp_tipnum = namenums[hooktree$tip.label == temp_tipname]
				list_of_daughter_tipnames_to_add_hooks_below = c(list_of_daughter_tipnames_to_add_hooks_below, hooktree$tip.label[temp_tipnum])
				
				#if (newtipnames == "tipnames")
				#	{
				hook_name_base = paste0("nodeBelow_", hooktree$tip.label[temp_tipnum])
				list_of_new_tipnames = c(list_of_new_tipnames, hook_name_base)
				#	list_of_daughter_tipnames_to_add_hooks_below = c(list_of_daughter_tipnames_to_add_hooks_below, hook_name_base)
				#	} # END if (newtipnames == "tipnames")


				#
				edgenums = 1:nrow(edge_times_bp)
				edgenum = edgenums[tr$edge[,2] == temp_tipnum]
				node_time_of_daughter = edge_times_bp[edgenum, 2]
				
				time_before_daughter_nodes = time_slice - node_time_of_daughter

				list_of_ages_below_daughter = c(list_of_ages_below_daughter, time_before_daughter_nodes)
				list_of_ages_absolute = c(list_of_ages_absolute, time_slice)
				} # END if (nodenum <= ntips)
			
			if (printflag >= 3)
				{
				print("print(list_of_ages_below_daughter):")
				print(list_of_ages_below_daughter)
				}
			}
		}
		

	if (printflag >= 1.5)
		{
		txt = paste0("\nAdding ", length(list_of_daughter_tipnames_to_add_hooks_below), " hooks to tree: ")
		cat(txt)
		}

	
	# Now, attach the hooks
	for (i in 1:length(list_of_daughter_tipnames_to_add_hooks_below))
		{
		#print(paste("i=", i, sep=""))
		tipname = list_of_daughter_tipnames_to_add_hooks_below[i]
		depthtime = as.numeric(list_of_ages_below_daughter[i])
		
		if (printflag >= 3)
			{
			print("print(depthtime):")
			print(depthtime)
			print("print(list_of_ages_below_daughter):")
			print(list_of_ages_below_daughter)
			}

		if (printflag >= 1.5)
			{
			txt = paste0(i, " ")
			cat(txt)
			}
			
		if (newtipnames == "default")
			{
			hooktree = add_hook(hooktree, tipname, depthtime, plottree=plottree, printflag=printflag, newtipname=newtipnames)
			}
		if (newtipnames == "tipnames")
			{
			new_name = paste0(list_of_new_tipnames[i], "_", list_of_ages_absolute[i])
			hooktree = add_hook(hooktree, tipname, depthtime, plottree=plottree, printflag=printflag, newtipname=new_name)
			}
		}
	
	return(hooktree)
	}



#' Find the branches that cross a particular timepoint
#'
#' This function finds the branches (edges) that cross a particular
#' timepoint. The input is a numedges x 2 matrix, with the left
#' column specifying the age of bottom (oldest end) of the branch,
#' and the right column specifying the age of the top (youngest end)
#' of the branch.  The input (edge_times_bp) can be obtained, 
#' with a little effort, from the \code{prt} function (prints tree to a table).
#'
#' @param time_slice The timepoint of interest (a single number).
#' @param edge_times_bp A numedges x 2 matrix, with the left
#' column specifying the age of bottom (oldest end) of the branch,
#' and the right column specifying the age of the top (youngest end)
#' of the branch.
#' @param The number of digits to round the branch start/end times to. The default
#'        (roundto=5) equates to 10 calendar years, if the unit of the branchlengths
#'        is millions of years, as is typical. The rounding can be helpful if 
#'        the living tips do not up to exact 0 my before present, due to e.g. 
#'        rounding errors in calculating the summary tree, and if it important 
#'        whether or not lineages exist at timeslice=0.0 before present.
#' @return edges_that_exist_in_the_right_time A TRUE/FALSE vector.
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
#' # Load tree and set up edge times matrix
#' tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
#' trtable = prt(tr, printflag=FALSE)
#' branch_tops = trtable$time_bp
#' branch_bots = trtable$time_bp + trtable$edge.length
#' branchnum = trtable$parent_br
#' edge_times_bp = cbind(branch_bots, branch_tops)
#' edge_times_bp
#' 
#' # Remove root node, which has no branch
#' edge_times_bp = edge_times_bp[!is.na(branchnum), ]
#' branchnum = branchnum[!is.na(branchnum)]
#' # Sort by branch number
#' reorder_vals = order(branchnum)
#' edge_times_bp = edge_times_bp[reorder_vals,]
#' row.names(edge_times_bp) = branchnum[reorder_vals]
#' edge_times_bp
#' 
#' edges_existing_at_correct_time_bp_TF(time_slice=0.5, edge_times_bp, roundto=5)
#' 
#' edges_existing_at_correct_time_bp_TF(time_slice=1.5, edge_times_bp, roundto=5)
#' 
#' edges_existing_at_correct_time_bp_TF(time_slice=0.0, edge_times_bp, roundto=5)
#' 
edges_existing_at_correct_time_bp_TF <- function(time_slice, edge_times_bp, roundto=5)
	{
	ex='
	# Load tree and set up edge times matrix
	tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
	trtable = prt(tr, printflag=FALSE)
	branch_tops = trtable$time_bp
	branch_bots = trtable$time_bp + trtable$edge.length
	branchnum = trtable$parent_br
	edge_times_bp = cbind(branch_bots, branch_tops)
	edge_times_bp

	# Remove root node, which has no branch
	edge_times_bp = edge_times_bp[!is.na(branchnum), ]
	branchnum = branchnum[!is.na(branchnum)]
	# Sort by branch number
	reorder_vals = order(branchnum)
	edge_times_bp = edge_times_bp[reorder_vals,]
	row.names(edge_times_bp) = branchnum[reorder_vals]
	edge_times_bp

	edges_existing_at_correct_time_bp_TF(time_slice=0.5, edge_times_bp, roundto=5)

	edges_existing_at_correct_time_bp_TF(time_slice=1.5, edge_times_bp, roundto=5)

	edges_existing_at_correct_time_bp_TF(time_slice=0.0, edge_times_bp, roundto=5)
	'
	# Needed for add_hooks:
	# find the edges that exist in the right time
	
	# (note: round to a default of 5 digits (a single year) with roundto; this is 
	#  important for whether or not lineages exist at time=0 before present)

	# timepoint is younger or equal to the oldest end of the branch
	edges_that_start_below_time = round(edge_times_bp[, 1], digits=roundto) > time_slice
	
	# timepoint is older than the youngest end of the branch	
	edges_that_end_after_time = round(edge_times_bp[, 2], digits=roundto) <= time_slice
	edges_that_exist_in_the_right_time = edges_that_start_below_time + edges_that_end_after_time == 2
	return(edges_that_exist_in_the_right_time)
	}


