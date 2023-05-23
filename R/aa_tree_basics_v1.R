

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













motA_tree_example <- function()
	{
	trstr = "(((((((((((((((((((((((('AAC73831.1_TolPal_system_protein_TolQ_Escherichi':2e-06,'QKY92913.1_TolPal_system_protein_TolQ_Shigella_s':2e-06):0.069602,'QXB11428.1_TolPal_system_protein_TolQ_Klebsiella':0.015345):0.035819,('AJJ55217.1_protein_TolQ_Yersinia_pseudotuberculo':2e-06,'AJJ31712.1_protein_TolQ_Yersinia_pestis_':2e-06):0.138185):0.398956,('BAC59320.1_TolQ_protein_Vibrio_parahaemolyticus_':0.344151,'QEO41016.1_protein_TolQ_Vibrio_cholerae_':0.120534):0.174246):0.320143,'ADY82611.1_tolerance_to_group_A_colicins_singles':1.14206):0.340639,'AAG04358.1_TolQ_protein_Pseudomonas_aeruginosa_P':0.471411):0.357661,((('APW43236.1_protein_TolQ_Rhodoferax_saidenbachens':0.432405,'ATU66653.1_protein_TolQ_Rhizobacter_gummiphilus_':0.450307):0.166467,('QKS27990.1_MAG_protein_TolQ_Candidatus_Accumuli':0.448172,'UBQ04634.1_protein_TolQ_Tepidimonas_taiwanensis_':0.34701):0.266215):0.42606,(('MBU2984437.1_protein_TolQ_Saccharophagus_degrada':2e-06,'ABD81792.1_MotA/TolQ/ExbB_proton_channel_Sacchar':2e-06):0.988781,'QEI18804.1_protein_TolQ_Cellvibrio_japonicus_':0.656449):0.270271):0.018258):0.364169,((((('AIQ88683.1_TolQ_protein_Methylobacterium_oryzae_':0.354101,'ABA05973.1_Cell_division_and_transportassociated_':0.28909):0.109819,'ABQ61581.1_protein_TolQ_Brucella_ovis_ATCC_25840':0.524003):0.209494,(('AEO47594.1_MotA/TolQ/ExbB_proton_channel_Rhodosp':0.432166,'CDK97994.1_TolQ_protein_Magnetospirillum_gryphis':0.199811):0.082287,'AIK95558.1_biopolymer_transporter_ExbB_Candidatu':0.801192):0.247759):0.412765,('ADE86923.1_TolQ_protein_Rhodobacter_capsulatus_S':0.223393,'WP_002720852.1_protein_TolQ_Cereibacter_sphaeroi':0.366171):0.412355):0.119263,(('AEI88825.1_Tolq_transport_protein_Candidatus_Mid':0.939789,'BDB95989.1_hypothetical_protein_HYD_1220_Candida':1.55138):0.042938,('AIQ91305.1_TolQ_protein_Methylobacterium_oryzae_':0.522723,'AEO48685.1_MotA/TolQ/ExbB_proton_channel_Rhodosp':0.859677):0.584851):0.073635):0.31233):0.066741,(((((((((((('AAC76042.1_Ton_complex_subunit_ExbB_Escherichia_':2e-06,'QKY94846.1_tolpal_systemassociated_acylCoA_thioes':2e-06):0.258223,'QXB09426.1_tolpal_systemassociated_acylCoA_thioes':0.431278):0.336062,('AJJ56099.1_tonBsystem_energizer_ExbB_Yersinia_ps':3e-06,'AJJ31614.1_tonBsystem_energizer_ExbB_Yersinia_pe':0.004788):0.46789):0.89846,'ADE86105.1_biopolymer_transport_protein_ExbB_Rho':1.06068):0.039342,(('ABA04528.1_outer_membrane_transport_energization_':2e-06,'ABA05559.1_outer_membrane_transport_energization_':0.014491):0.780174,('ABA05299.1_outer_membrane_transport_energization_':2e-06,'ABA03969.1_outer_membrane_transport_energization_':2e-06):0.65744):0.149302):2e-06,(('WP_015921325.1_tonBsystem_energizer_ExbB_Cereiba':1.19569,'ABQ60495.1_ExbB_Brucella_ovis_ATCC_25840_':0.368682):0.183052,'CDK97637.1_membrane_spanning_protein_in_TonBExbBE':0.831166):0.066221):0.645544,('ACL95886.1_TonB_accessory_protein_exbB_Caulobact':0.147727,'ADL01024.1_MotA/TolQ/ExbB_proton_channel_Brevund':0.30299):0.850101):0.262202,('ADL00174.1_protein_TolQ_Brevundimonas_subvibrioi':0.630686,'ACL96805.1_TolQ_protein_Caulobacter_vibrioides_N':0.431295):0.83407):0.251321,'AHZ85825.1_adventurous_gliding_motility_protein_X':2.64646):3e-06,'BAG40198.1_TolQ_protein_Orientia_tsutsugamushi_s':1.47687):0.027664,'QNQ08778.1_protein_TolQ_Sphingomonas_alpina_':0.848102):0.338104,'AUX20279.1_transport_protein_tolQ_Sorangium_cel':1.89005):0.266819):0.28453,(((((((((('QEI19895.1_MotA/TolQ/ExbB_proton_channel_family_p':0.50147,'AAG03587.1_transport_protein_ExbB_Pseudomonas_ae':0.493563):0.12735,('ADY81763.1_putative_biopolymer_transport_protein_':0.17744,'ADY81773.1_putative_biopolymer_transport_protein_':0.229745):0.557315):0.736719,('APW44370.1_flagellar_motor_protein_MotA_Rhodofer':0.59549,'ATU67891.1_MotA/TolQ/ExbB_proton_channel_family_p':0.463267):0.684262):0.310969,'ATU63282.1_biopolymer_transporter_ExbB_Rhizobact':1.25497):0.100491,(('ADY81560.1_putative_biopolymer_transport_protein_':0.205257,'ATU64277.1_MotA/TolQ/ExbB_proton_channel_family_p':0.43972):0.242199,'QNQ10683.1_MotA/TolQ/ExbB_proton_channel_family_p':1.27365):0.631488):0.036614,(('QEI19039.1_MotA/TolQ/ExbB_proton_channel_family_p':0.490647,'APW42557.1_biopolymer_transporter_Rhodoferax_sai':0.690305):0.089273,'UBQ06174.1_MotA/TolQ/ExbB_proton_channel_family_p':0.824145):0.504873):2e-06,'SNW07143.1_biopolymer_transport_protein_ExbB_Eik':1.70777):0.046461,'QKS29082.1_MAG_MotA/TolQ/ExbB_proton_channel_fam':1.66271):0.156351,'ADY83951.1_MotA/TolQ/ExbB_proton_channel_Acineto':1.65909):0.652219,(((((((((('AOO66396.1_TolQ_proton_channel_family_protein_Su':0.629765,'ADG93510.1_tonBsystem_energizer_ExbB_Arcobacter_':0.231833):0.814538,'ADN08701.1_outer_membrane_transport_energization_':0.82404):0.07116,'AKF25943.1_biopolymer_transporter_ExbB_Sulfurovu':0.894162):0.070563,(('ALF47140.1_TonB_system_transport_protein_ExbB_Ca':0.214746,'CAL35725.1_putative_exbB/tolQ_family_transport_pr':0.662943):0.585602,'AOO66288.1_ferric_siderophore_transport_system_b':0.634973):0.550846):0.201748,('ADO46046.1_tonBsystem_energizer_ExbB_Hydrogenoba':0.34993,'ABC78648.1_tolQ_protein_Syntrophus_aciditrophicu':0.509458):0.058795):0.123019,'AHE95642.1_biopolymer_transporter_ExbB_Thermocri':0.372539):0.071807,'AAC07592.1_biopolymer_transport_exbB_Aquifex_aeo':0.510379):0.448808,'ADD67695.1_tonBsystem_energizer_ExbB_Denitrovibr':0.610274):0.235899,('SLM48220.1_Biopolymer_transport_protein_ExbB_Nit':2e-06,'ALA57346.1_Biopolymer_transport_protein_exbB_Nit':0.060497):0.448261):1.03357,'AOO66426.1_TolQ_proton_channel_family_protein_Su':0.451963):1.04702):3e-06):0.233863,((((((((((((('ABF85904.1_MotA/TolQ/ExbB_proton_channel_family_p':2e-06,'UEO04322.1_MotA/TolQ/ExbB_proton_channel_family_p':2e-06):0.085479,('ABF90492.1_MotA/TolQ/ExbB_proton_channel_family_p':2e-06,'UEO05991.1_MotA/TolQ/ExbB_proton_channel_family_p':2e-06):0.096728):0.089314,('ADO75002.1_MotA/TolQ/ExbB_proton_channel_family_p':0.142425,'ADO76074.1_MotA/TolQ/ExbB_proton_channel_family_p':0.154506):0.03206):0.067138,'AFE08685.1_MotA/TolQ/ExbB_proton_channel_family_p':0.207663):0.655108,('AUX23532.1_biopolymer_transporter_Sorangium_cell':0.877703,'AWV89332.1_MotA/TolQ/ExbB_proton_channel_family_p':0.942067):0.210157):2e-06,'ACL64950.1_MotA/TolQ/ExbB_proton_channel_Anaerom':1.02824):0.294827,'UPT75425.1_MAG_MotA/TolQ/ExbB_proton_channel_fam':0.92655):0.140998,('UPT76060.1_MAG_MotA/TolQ/ExbB_proton_channel_fam':0.895438,'EKT87439.2_TolQ_transporter_Leptospira_santarosa':1.14831):0.404557):2e-06,'UPT74212.1_MAG_MotA/TolQ/ExbB_proton_channel_fam':1.38813):0.198775,'AHZ86092.1_tolQ_protein_Bdellovibrio_bacteriovor':1.44852):3e-06,(('QQS07012.1_MAG_MotA/TolQ/ExbB_proton_channel_fam':0.129346,'QQS07007.1_MAG_MotA/TolQ/ExbB_proton_channel_fam':0.116091):1.29881,'QQS07765.1_MAG_MotA/TolQ/ExbB_proton_channel_fam':0.955966):0.95838):0.183146,((((('UEO07578.1_MotA/TolQ/ExbB_proton_channel_family_p':2e-06,'ABF91342.1_MotA/TolQ/ExbB_proton_channel_family_p':2e-06):0.194194,'AFE10516.1_MotA/TolQ/ExbB_proton_channel_family_p':0.149439):0.05875,'ADO73220.1_MotA/TolQ/ExbB_proton_channel_family_p':0.200083):0.501726,'ACL67494.1_MotA/TolQ/ExbB_proton_channel_Anaerom':0.807212):0.554986,'UPT72845.1_MAG_MotA/TolQ/ExbB_proton_channel_fam':1.83281):0.520921):0.176672,('AEO48101.1_chemotaxis_sensory_transducer_Rhodosp':0.063406,'AEO48100.1_chemotaxis_sensory_transducer_Rhodosp':0.002513):3.453):0.260945):0.047767,(((((((((('BAH38987.1_putative_biopolymer_transport_protein_':0.018302,'AMW06740.1_hypothetical_protein_GEMMAAP_10105_Ge':0.327583):0.063863,'QJR34845.1_hypothetical_protein_HKW67_04610_Gemm':0.713004):1.29787,('ATX81939.1_Cell_division_and_transportassociated_':0.071255,'ATX79746.1_Cell_division_and_transportassociated_':0.186058):0.518032):0.402454,((('ACI22061.1_TolQ_protein_Thermodesulfovibrio_yell':0.605273,'BCB95125.1_protein_TolQ_Dissulfurispira_thermoph':0.54296):0.332587,'ABC76751.1_tolQ_protein_Syntrophus_aciditrophicu':0.759051):0.209644,'ACX73814.1_MotA/TolQ/ExbB_proton_channel_Fibroba':2.08833):0.002675):0.153132,'AWV88602.1_flagellar_motor_protein_MotA_Bradymon':0.8054):0.204279,(((('AFL88265.1_biopolymer_transport_protein_Terriglo':2e-06,'AFL88606.1_biopolymer_transport_protein_Terriglo':2e-06):0.617671,'ACO32834.1_transporter_MotA/TolQ/ExbB_proton_cha':0.536337):0.506285,'QOY87643.1_MotA/TolQ/ExbB_proton_channel_family_p':1.0962):0.566472,'QTD52661.1_MotA/TolQ/ExbB_proton_channel_family_p':1.18207):0.258228):0.078483,((((('UEO06636.1_MotA/TolQ/ExbB_proton_channel_family_p':2e-06,'ABF91581.1_transporter_MotA/TolQ/ExbB_proton_cha':2e-06):0.034464,'AFE07999.1_MotA/TolQ/ExbB_proton_channel_family_p':0.109091):0.412265,'ADO74176.1_Adventurous_gliding_motility_protein_X':0.109336):0.501078,'ACL64052.1_protein_TolQ_Anaeromyxobacter_dehalog':0.301497):0.45502,('CAL34348.1_biopolymer_transport_protein_Campylob':1.8191,'AAR33363.1_biopolymer_transport_membrane_proton_c':0.602386):0.155984):0.143531):0.060876,((((('ALA60092.1_Protein_TolQ_Nitrospira_moscoviensis_':0.38678,'SLM46359.1_Protein_TolQ_Nitrospira_japonica_':0.483768):0.661906,'BAM06957.1_putative_MotA/TolQ/ExbB_proton_channel':1.29099):0.159585,'QOY89541.1_MotA/TolQ/ExbB_proton_channel_family_p':3.81491):0.251212,('AQV01680.1_protein_TolQ_Desulfococcus_multivoran':0.648969,'CAG37138.1_related_to_biopolymer_transport_protei':0.82466):0.383878):0.027504,(('ACF44367.1_protein_TolQ_Pelodictyon_phaeoclathra':0.56662,'ACF11057.1_protein_TolQ_Chlorobaculum_parvum_NCI':0.433568):0.183337,('ABL64907.1_Cell_division_and_transportassociated_':0.26398,'ACD89825.1_protein_TolQ_Chlorobium_limicola_DSM_':0.25662):0.491208):0.668154):0.08459):0.375425,((((((((('AGW38400.1_putative_transport_protein_Chlamydia_':0.160061,'AAP98742.1_ExbB_Chlamydia_pneumoniae_TW183_':0.192181):0.094634,('AVM88605.1_MotA/TolQ/ExbB_proton_channel_family_p':0.070429,'AAC68198.1_polysaccharide_transporter_Chlamydia_':0.025028):0.284105):1.5303,'CCB87689.1_putative_uncharacterized_protein_Para':1.01923):1.41155,((('QJE97904.1_DUF2341_domaincontaining_protein_Lute':0.512494,'QJE95399.1_hypothetical_protein_HHL09_06260_Lute':0.198523):0.707168,'AAG04082.1_transport_protein_ExbB2_Pseudomonas_a':1.39727):0.042735,'QNQ08565.1_DUF2341_domaincontaining_protein_Sphi':0.468215):1.24757):2e-06,(('ACB76604.1_MotA/TolQ/ExbB_proton_channel_Opitutu':0.473718,'ATC65746.1_flagellar_motor_protein_MotA_Nibricoc':0.110352):0.243871,'AWI09757.1_flagellar_motor_protein_MotA_Ereboglo':0.146099):1.00583):0.250037,('AKJ63346.1_Biopolymer_transport_protein_ExbB_Kir':0.946912,'AVM45096.1_hypothetical_protein_C5Q97_10465_Vict':1.20075):0.249423):2e-06,((('ALA61164.1_MotA/TolQ/ExbB_proton_channel_Nitrosp':1.19337,DGBBEKCF_00127_TolPal_system_protein_TolQ_Verruc:0.684914):0.084122,'ALJ56928.1_Biopolymer_transport_protein_ExbB_Can':0.901312):3e-06,'QSR84055.1_MotA/TolQ/ExbB_proton_channel_family_p':0.455245):0.505717):0.661057,'AGR80198.1_TolQ_protein_Anaplasma_phagocytophilu':2.15842):0.018765,('EKT87417.1_flagellar_motor_protein_MotA_Leptospi':3.57174,'QQS05379.1_MAG_MotA/TolQ/ExbB_proton_channel_fam':1.67268):0.464054):0.321069):3e-06,(((((('AOO64467.1_putative_biopolymer_transport_protein_':0.697758,'ALF47162.1_TolPal_system_subunit_TolQ_Campylobac':0.240941):0.210073,'CAL34280.1_putative_MotA/TolQ/ExbB_proton_channel':0.514531):0.488635,'AKF24510.1_biopolymer_transporter_Sulfurovum_lit':0.525099):0.147449,'ADG93991.1_MotA/TolQ/ExbB_proton_channel_Arcobac':0.619506):0.054094,'ADN08827.1_MotA/TolQ/ExbB_proton_channel_Sulfuri':0.770756):0.961298,'QTD51577.1_MotA/TolQ/ExbB_proton_channel_family_p':5.81166):1.4159):0.149761):0.345047,(((((((('ABF92363.1_transporter_MotA/TolQ/ExbB_proton_cha':2e-06,'UEO03725.1_MotA/TolQ/ExbB_proton_channel_family_p':2e-06):0.075847,'AFE10515.1_MotA/TolQ/ExbB_proton_channel_family_p':0.04689):0.03988,'ADO70054.1_TolQ_protein_Stigmatella_aurantiaca_D':0.041059):0.249848,'AUX23472.1_flagellar_motor_protein_MotA_Sorangiu':0.381529):0.650582,((('ABF92421.1_tonB_system_transport_protein_ExbB/Tol':2e-06,'UEO04809.1_MotA/TolQ/ExbB_proton_channel_family_p':2e-06):0.190041,'AFE08686.1_TonB_system_transport_protein_ExbB/Tol':0.154783):0.729317,'AUX24210.1_flagellar_motor_protein_MotA_Sorangiu':0.837449):0.434782):0.072534,'ACL63649.1_MotA/TolQ/ExbB_proton_channel_Anaerom':0.741689):0.286499,((('AFL86415.1_biopolymer_transport_protein_Terriglo':0.145127,'ACO31817.1_transporter_MotA/TolQ/ExbB_proton_cha':0.118209):0.129677,'QOY88067.1_MotA/TolQ/ExbB_proton_channel_family_p':0.290722):0.134529,'QUV79767.1_MotA/TolQ/ExbB_proton_channel_family_p':0.312933):0.242717):0.072065,('BAH39052.1_putative_biopolymer_transport_protein_':0.058868,'QJR34785.1_flagellar_motor_protein_MotA_Gemmatim':0.104988):0.938522):0.609225):0.143988,'UPT72618.1_MAG_MotA/TolQ/ExbB_proton_channel_fam':1.56685):0.178704,((('BAC61768.1_ExbBlike_protein_Vibrio_parahaemolyti':0.379994,'QEO43084.1_MotA/TolQ/ExbB_proton_channel_family_p':0.148335):2.77307,'BAC89082.1_gll1141_Gloeobacter_violaceus_PCC_742':1.37768):0.396201,'AKL98121.1_putative_Transporter_MotA/TolQ/ExbB_p':1.3107):0.159921):0.266638,'BAS27557.1_MotA/TolQ/ExbB_proton_channel_family_p':1.5998):0.057561,(((((((((('ABL64659.1_outer_membrane_transport_energization_':0.087037,'ACD90995.1_MotA/TolQ/ExbB_proton_channel_Chlorob':0.109938):0.064872,'ACF43024.1_MotA/TolQ/ExbB_proton_channel_Pelodic':0.055476):0.078225,'ACF11913.1_MotA/TolQ/ExbB_proton_channel_Chlorob':0.043534):0.363544,'ACF14782.1_MotA/TolQ/ExbB_proton_channel_Chloroh':0.311314):0.168301,'ACF14374.1_MotA/TolQ/ExbB_proton_channel_Chloroh':0.882939):0.338866,((('ABG59284.1_outer_membrane_transport_energization_':0.583005,'AKD02510.1_flagellar_motor_protein_MotA_Pontibac':0.387608):0.040708,'AEW03299.1_outer_membrane_transport_energization_':0.555337):0.084073,('BCI63438.1_flagellar_motor_protein_MotA_Coprobac':0.200661,'UAK41566.1_MotA/TolQ/ExbB_proton_channel_family_p':0.314469):0.400741):0.21486):0.880384,'QVL32859.1_MotA/TolQ/ExbB_proton_channel_family_p':2.20341):0.267963,((((('CAD76944.1_probable_TolQ_Rhodopirellula_baltica_':0.672355,'QDU94732.1_colicin_uptake_protein_TolQ_Lignipire':0.416108):0.178091,'ADB17385.1_MotA/TolQ/ExbB_proton_channel_Pirellu':0.42968):0.357844,'QDU27583.1_colicin_uptake_protein_TolQ_Anatilimn':1.2654):0.988197,'ADY58078.1_MotA/TolQ/ExbB_proton_channel_Rubinis':0.984922):0.43534,'ADV63710.1_MotA/TolQ/ExbB_proton_channel_Isospha':1.48435):1.20805):0.067652,'QTD53280.1_MotA/TolQ/ExbB_proton_channel_family_p':0.97588):0.01644,(((('ABG60812.1_outer_membrane_transport_energization_':0.802211,'AEW00830.1_outer_membrane_transport_energization_':0.737497):0.22059,'AKD04823.1_biopolymer_transporter_ExbB_Pontibact':0.606008):0.091494,('BCI63565.1_biopolymer_transporter_ExbB_Coprobact':0.286298,'UAK41330.1_MotA/TolQ/ExbB_proton_channel_family_p':0.250399):0.424589):0.629313,(('ADG92738.1_MotA/TolQ/ExbB_proton_channel_Arcobac':0.840677,'ADG92304.1_MotA/TolQ/ExbB_proton_channel_Arcobac':0.447437):0.514819,'ACX75330.1_MotA/TolQ/ExbB_proton_channel_Fibroba':1.36836):1.12859):0.190464):0.114134):0.040425,(((((((((((('SEH78546.1_mota/tolq/exbb_proton_channel_family_':1.28181,'QJE97269.1_MotA/TolQ/ExbB_proton_channel_family_p':1.11377):0.426751,('QSR84464.1_MotA/TolQ/ExbB_proton_channel_family_p':0.460866,DGBBEKCF_01299_TolPal_system_protein_TolQ_Verruc:0.538879):0.670377):0.428623,('ADV62354.1_hypothetical_protein_Isop_1771_Isosph':2.18463,'ADV60824.1_MotA/TolQ/ExbB_proton_channel_Isospha':1.19737):0.097779):0.060076,((('AWI09431.1_flagellar_motor_protein_MotA_Ereboglo':0.134429,'ATC62595.1_flagellar_motor_protein_MotA_Nibricoc':0.507827):0.832644,'ACB77180.1_MotA/TolQ/ExbB_proton_channel_Opitutu':0.043195):1.62577,'AQU99647.1_MotA_Desulfococcus_multivorans_':2.53836):0.403566):0.429735,((((('ADB17134.1_MotA/TolQ/ExbB_proton_channel_Pirellu':0.127099,'QDU25761.1_Biopolymer_transport_protein_ExbB_Ana':0.361853):0.300156,'QDU93364.1_Biopolymer_transport_protein_ExbB_Lig':0.368495):0.161857,'ADY58376.1_MotA/TolQ/ExbB_proton_channel_Rubinis':0.662233):0.077234,'CAD71937.1_probable_biopolymer_transport_ExbBrela':0.773593):0.645433,'AWV88047.1_MotA/TolQ/ExbB_proton_channel_family_p':1.03469):0.169576):0.131007,(((((('MBU2985725.1_MotA/TolQ/ExbB_proton_channel_family':0.001762,'ABD81321.1_MotA/TolQ/ExbB_proton_channel_Sacchar':0.0224):0.508966,'QEI19442.1_MotA/TolQ/ExbB_proton_channel_family_p':0.287323):0.292793,'AAG06371.1_probable_tolQtype_transport_protein_P':0.597557):0.279103,'ADY81529.1_hypothetical_protein_BDGL_000943_Acin':1.39259):0.044342,((('UBQ05954.1_MotA/TolQ/ExbB_proton_channel_family_p':0.405313,'ATU65616.1_MotA/TolQ/ExbB_proton_channel_family_p':0.363274):0.04759,'APW42522.1_flagellar_motor_protein_MotA_Rhodofer':0.350807):0.527396,'QKS29102.1_MAG_MotA/TolQ/ExbB_proton_channel_fam':0.383985):0.557229):0.274961,('ATX80498.1_biopolymer_transport_protein_ExbB_Mar':0.059595,'ATX81215.1_biopolymer_transport_protein_ExbB_Mar':3e-06):0.795094):0.279843):0.146973,(((('ADD67874.1_MotA/TolQ/ExbB_proton_channel_Denitro':0.420642,'QAR32927.1_MotA/TolQ/ExbB_proton_channel_family_p':0.383857):0.107429,'AEI15149.1_MotA/TolQ/ExbB_proton_channel_Flexist':0.327872):0.186835,'USF25003.1_Biopolymer_transport_protein_ExbB_Muc':0.527987):0.54538,'QTD50766.1_MotA/TolQ/ExbB_proton_channel_family_p':1.13664):0.298264):0.083828,((('ADD67609.1_MotA/TolQ/ExbB_proton_channel_Denitro':0.547118,'QAR33430.1_MotA/TolQ/ExbB_proton_channel_family_p':0.723745):1.13749,'BCB95942.1_biopolymer_transporter_ExbB_Dissulfur':1.38896):0.358829,'SEI00245.1_mota/tolq/exbb_proton_channel_family_':1.71674):0.08372):0.163862,(((((((('ACB75629.1_MotA/TolQ/ExbB_proton_channel_Opitutu':0.408411,'AWI08882.1_flagellar_motor_protein_MotA_Ereboglo':0.813975):3e-06,'ATC64798.1_flagellar_motor_protein_MotA_Nibricoc':0.450255):0.062704,'ATC62712.1_flagellar_motor_protein_MotA_Nibricoc':0.263953):0.629389,(('AKJ65464.1_colicin_uptake_protein_TolQ_Kiritimat':0.49366,'ARN57240.1_Biopolymer_transport_protein_ExbB_Sed':0.780033):2.08989,'ATU63381.1_hypothetical_protein_CPZ87_01790_Rhiz':3.19187):2e-06):0.260216,((('QSR85862.1_MotA/TolQ/ExbB_proton_channel_family_p':0.765489,PKPEBJJI_00344_hypothetical_protein_CAIZXV01_Ver:0.521677):0.346807,('QSR84237.1_MotA/TolQ/ExbB_proton_channel_family_p':0.630927,'QSR83870.1_MotA/TolQ/ExbB_proton_channel_family_p':0.27706):0.71318):0.13978,('AKJ64362.1_Biopolymer_transport_protein_ExbB_Kir':0.99965,'QJE97567.1_MotA/TolQ/ExbB_proton_channel_family_p':1.72602):3e-06):0.010961):0.114118,'BAM06053.1_biopolymer_transport_protein_Leptospi':1.94991):0.209361,((((('QZZ23337.1_MotA/TolQ/ExbB_proton_channel_family_p':0.549063,'QIZ70555.1_MotA/TolQ/ExbB_proton_channel_family_p':0.407):0.076934,'UNU17462.1_MotA/TolQ/ExbB_proton_channel_family_p':0.349549):0.084183,'QZZ22141.1_MotA/TolQ/ExbB_proton_channel_family_p':0.384591):0.521419,('BAC90343.1_MotA/TolQ/ExbB_family_proton_channel_p':2e-06,'BAC89328.1_glr1387_Gloeobacter_violaceus_PCC_742':2e-06):0.787608):0.239824,((('QZZ19520.1_MotA/TolQ/ExbB_proton_channel_family_p':0.644929,'QIZ71789.1_MotA/TolQ/ExbB_proton_channel_family_p':0.263144):0.249408,'UNU18280.1_MotA/TolQ/ExbB_proton_channel_family_p':0.184539):0.418415,'QZZ21680.1_MotA/TolQ/ExbB_proton_channel_family_p':0.918197):0.500123):0.688398):0.596947,(('QVL33094.1_MotA/TolQ/ExbB_proton_channel_family_p':6.51399,'CAG37702.1_related_to_biopolymer_transport_protei':1.45167):2e-06,('ATX79248.1_outer_membrane_transport_energization_':0.227101,'ATX82358.1_outer_membrane_transport_energization_':0.160104):1.36086):0.753911):0.06271):0.049451,(((((DGBBEKCF_01621_TolPal_system_protein_TolQ_Verruc:0.508706,'QSR84174.1_MotA/TolQ/ExbB_proton_channel_family_p':0.60196):0.492843,'ALJ56692.1_colicin_uptake_protein_TolQ_Candidatu':0.828523):2e-06,('SEH96646.1_mota/tolq/exbb_proton_channel_family_':0.829807,'QJE97381.1_MotA/TolQ/ExbB_proton_channel_family_p':0.817261):0.791752):0.308196,('CCB87487.1_putative_uncharacterized_protein_Para':1.15512,'AKJ65567.1_Biopolymer_transport_protein_ExbB_Kir':0.885876):0.162557):0.020789,((('ATC63405.1_flagellar_motor_protein_MotA_Nibricoc':0.18928,'AWI08418.1_flagellar_motor_protein_MotA_Ereboglo':0.402942):0.109031,'ACB73533.1_MotA/TolQ/ExbB_proton_channel_Opitutu':0.280297):0.994296,('QAT17178.1_hypothetical_protein_BU251_05250_Cand':2.14349,'AVM45415.1_biopolymer_transporter_ExbB_Victivall':0.547125):0.430903):3e-06):0.520388):1.999999999e-06,((((((((((('BAC61496.1_TonB_system_transport_protein_ExbB2_V':0.485211,'QEO41300.1_MotA/TolQ/ExbB_proton_channel_family_p':0.325922):0.777833,'BAC58428.1_TonB_system_transport_protein_ExbB2_V':0.976754):0.147129,('MBU2986357.1_MotA/TolQ/ExbB_proton_channel_family':2e-06,'ABD79618.1_MotA/TolQ/ExbB_proton_channel_Sacchar':2e-06):0.613762):0.120488,'QEI18325.1_MotA/TolQ/ExbB_proton_channel_family_p':0.86031):0.23083,'ADG92230.1_MotA/TolQ/ExbB_proton_channel_Arcobac':0.950659):0.193054,'QEI18408.1_MotA/TolQ/ExbB_proton_channel_family_p':1.31441):1.54412,'ACX74512.1_MotA/TolQ/ExbB_proton_channel_Fibroba':1.77891):2e-06,'AKJ64799.1_ExbB_proton_channel_family_protein_Ki':2.42966):1.999999999e-06,'AEV99079.1_multisensor_hybrid_histidine_kinase_N':6.89379):0.747705,'EKT88184.1_biopolymer_transporter_Leptospira_san':1.3177):0.197423,'QAT17590.1_ferric_siderophore_transport_system_b':1.40156):0.275101):0.059778,(((((((((('QAR33783.1_MotA/TolQ/ExbB_proton_channel_family_p':1.16463,'ADD67745.1_MotA/TolQ/ExbB_proton_channel_Denitro':1.33992):1.63982,('ADD67744.1_MotA/TolQ/ExbB_proton_channel_Denitro':0.733009,'QAR33784.1_hypothetical_protein_EP073_10320_Geov':0.57442):3e-06):0.809866,('ADD68377.1_MotA/TolQ/ExbB_proton_channel_Denitro':0.791083,'AQU99646.1_MotA_Desulfococcus_multivorans_':0.614617):0.40194):0.06954,((('AKJ65465.1_Biopolymer_transport_protein_ExbB_Kir':0.405745,'ARN57239.1_Biopolymer_transport_protein_ExbB_Sed':0.425053):0.448894,'AKJ64800.1_Biopolymer_transport_protein_ExbB_Kir':0.912755):0.117603,('ABD82072.1_MotA/TolQ/ExbB_proton_channel_Sacchar':2e-06,'MBU2984729.1_MotA/TolQ/ExbB_proton_channel_family':0.005014):1.69697):0.544445):0.125259,(((('ABD79617.1_MotA/TolQ/ExbB_proton_channel_Sacchar':2e-06,'MBU2986356.1_MotA/TolQ/ExbB_proton_channel_family':2e-06):0.456921,'QEI18324.1_energy_transducer_TonB_Cellvibrio_jap':0.174728):0.474641,(('BAC58429.1_putative_TolR_Vibrio_parahaemolyticus':0.78751,'ADG92231.1_MotA/TolQ/ExbB_proton_channel_Arcobac':0.617082):0.19559,'QEI18407.1_MotA/TolQ/ExbB_proton_channel_family_p':0.8254):0.169289):1.999999999e-06,('BAC61495.1_biopolymer_transport_protein_ExbBrelat':0.330976,'QEO41299.1_MotA/TolQ/ExbB_proton_channel_family_p':0.188381):0.740235):0.62494):0.337896,(((('AVQ29260.1_MotA/TolQ/ExbB_proton_channel_family_p':0.850968,'QQS88496.1_MotA/TolQ/ExbB_proton_channel_family_p':0.707407):0.366874,'AVQ29523.1_MotA/TolQ/ExbB_proton_channel_family_p':0.778225):0.705983,'ACX74513.1_MotA/TolQ/ExbB_proton_channel_Fibroba':1.17645):1.999999999e-06,'ACX73657.1_MotA/TolQ/ExbB_proton_channel_Fibroba':1.70775):0.178501):0.18957,(((((((((('ABF90777.1_MotA/TolQ/ExbB_proton_channel_family_p':2e-06,'UEO05650.1_MotA/TolQ/ExbB_proton_channel_family_p':2e-06):0.13871,('AFE07995.1_MotA/TolQ/ExbB_proton_channel_family_p':0.097934,'ADO68798.1_MotA/TolQ/ExbB_proton_channel_family_p':0.227476):0.009103):0.285359,'ACL67782.1_MotA/TolQ/ExbB_proton_channel_Anaerom':0.282328):0.278602,'AUX22205.1_flagellar_motor_protein_MotA_Sorangiu':0.742597):0.164597,'AHZ86892.1_adventurous_gliding_motility_protein_R':1.01577):0.293861,('AHZ84904.1_gliding_motility_protein_Bdellovibrio':0.661783,'AHZ86042.1_gliding_motility_protein_Bdellovibrio':0.466693):0.635105):0.395406,(((('UEO02151.1_MotA/TolQ/ExbB_proton_channel_family_p':2e-06,'ABF87376.1_MotA/TolQ/ExbB_proton_channel_family_p':2e-06):0.097046,'ADO72879.1_MotA/TolQ/ExbB_proton_channel_family_p':0.241671):0.169098,'AFE07994.1_MotA/TolQ/ExbB_proton_channel_family_p':0.223134):0.712921,'AUX22204.1_flagellar_motor_protein_MotA_Sorangiu':0.754537):0.384318):2e-06,'AWV89791.1_MotA/TolQ/ExbB_proton_channel_family_p':0.802262):0.229713,(('ACX75484.1_MotA/TolQ/ExbB_proton_channel_Fibroba':1.73928,'AHZ84037.1_adventurous_gliding_motility_protein_R':1.151):0.204484,'AWV89461.1_MotA/TolQ/ExbB_proton_channel_family_p':1.24468):0.147047):0.072348,((('ABD82818.1_MotA/TolQ/ExbB_proton_channel_Sacchar':2e-06,'MBU2986144.1_MotA/TolQ/ExbB_proton_channel_family':2e-06):0.258385,'QEI20720.1_MotA/TolQ/ExbB_proton_channel_family_p':0.291143):0.214921,('ABD82214.1_MotA/TolQ/ExbB_proton_channel_Sacchar':0.004561,'MBU2985168.1_MotA/TolQ/ExbB_proton_channel_family':3e-06):0.619933):0.746528):0.412177):0.20147,(((('BBM36826.1_MotA/TolQ/ExbB_proton_channel_Pseudol':0.214906,'ACV38437.1_MotA/TolQ/ExbB_proton_channel_Leptotr':0.389921):0.95025,('QQS88170.1_MotA/TolQ/ExbB_proton_channel_family_p':0.730275,'AVQ28424.1_MotA/TolQ/ExbB_proton_channel_family_p':0.263072):0.439691):0.259989,(('CCG56832.1_biopolymer_transport_protein_ExbB_Bra':0.281014,'AEM22241.1_biopolymer_transport_protein_ExbB_Bra':0.231878):0.902093,'EKT86397.1_biopolymer_transporter_ExbB_Leptospir':1.40291):0.434102):0.214922,'UAK42981.1_MotA/TolQ/ExbB_proton_channel_family_p':1.96804):0.022134):3e-06,'QQS05296.1_MAG_MotA/TolQ/ExbB_proton_channel_fam':2.16259):0.153272,(('ADB17458.1_MotA/TolQ/ExbB_proton_channel_Pirellu':0.229529,'QDU31292.1_Biopolymer_transport_protein_ExbB_Ana':0.324079):0.42311,'QDU92956.1_Biopolymer_transport_protein_ExbB_Lig':0.325832):1.14652):0.092498):1.999999999e-06):0.167015,(('AEB12216.1_MotA/TolQ/ExbB_proton_channel_Marinit':0.686212,'ADI15492.1_MotA/TolQ/ExbB_proton_channel_Trueper':0.610874):0.315758,('ALW88052.1_biopolymer_transporter_Deinococcus_ac':0.257805,'ADV67905.1_MotA/TolQ/ExbB_proton_channel_Deinoco':0.645301):1.17542):0.479744):3e-06,(((((((((('ADB18835.1_MotA/TolQ/ExbB_proton_channel_Pirellu':0.269711,'QDU25322.1_Biopolymer_transport_protein_ExbB_Ana':0.279851):0.175588,'ADB15644.1_hypothetical_protein_Psta_0959_Pirell':2.70738):0.126662,'QDU98188.1_Biopolymer_transport_protein_ExbB_Lig':0.223822):0.141409,'CAD73493.1_probable_tolQ_protein_Rhodopirellula_':0.318097):0.496007,'ADY58081.1_MotA/TolQ/ExbB_proton_channel_Rubinis':0.370769):0.269227,(('QVL32579.1_MotA/TolQ/ExbB_proton_channel_family_p':1.07219,'ADV63187.1_MotA/TolQ/ExbB_proton_channel_Isospha':0.893668):0.156438,'ADV63892.1_MotA/TolQ/ExbB_proton_channel_Isospha':0.786668):0.416206):0.599325,(((('ADB15652.1_MotA/TolQ/ExbB_proton_channel_Pirellu':0.260609,'QDU98468.1_Biopolymer_transport_protein_ExbB_Lig':0.462796):0.129747,'QDU30851.1_Biopolymer_transport_protein_ExbB_Ana':0.499234):0.060761,'CAD74695.1_probable_TolQtype_transport_protein_R':0.953448):0.447448,'AKJ63368.1_Biopolymer_transport_protein_ExbB_Kir':1.66947):0.427361):0.022313,(((('ADB17547.1_MotA/TolQ/ExbB_proton_channel_Pirellu':0.256055,'QDU29130.1_colicin_uptake_protein_TolQ_Anatilimn':0.414425):0.173349,'QDU99088.1_colicin_uptake_protein_TolQ_Lignipire':0.427615):0.128746,'CAD77369.1_probable_tolQtype_transport_protein_R':0.908057):1.2292,(('QDU94893.1_Biopolymer_transport_protein_ExbB_Lig':0.527887,'ADB16663.1_MotA/TolQ/ExbB_proton_channel_Pirellu':0.450664):0.692154,'ADY61978.1_MotA/TolQ/ExbB_proton_channel_Rubinis':1.08251):0.207087):0.064061):0.251502,(((('QSR84626.1_MotA/TolQ/ExbB_proton_channel_family_p':0.69498,DGBBEKCF_00463_TolPal_system_protein_TolQ_Verruc:0.422673):0.599138,'QJE95191.1_MotA/TolQ/ExbB_proton_channel_family_p':0.818213):0.746812,('QJE97562.1_MotA/TolQ/ExbB_proton_channel_family_p':1.26025,PKPEBJJI_00297_TolPal_system_protein_TolQ_CAIZXV:1.11703):0.440007):2e-06,'AKJ64099.1_Biopolymer_transport_protein_ExbB_Kir':1.32081):0.182503):0.457852,((((('AHE96346.1_TolQ_Thermocrinis_ruber_':0.391164,'ADC89509.1_MotA/TolQ/ExbB_proton_channel_Thermoc':0.328045):0.158293,'ADO44512.1_MotA/TolQ/ExbB_proton_channel_Hydroge':0.329873):0.462357,'AAC07752.1_TolQlike_protein_Aquifex_aeolicus_VF5':0.547087):1.25946,('ARN57520.1_Biopolymer_transport_protein_ExbB_Sed':2.07262,'SEH93354.1_mota/tolq/exbb_proton_channel_family_':1.97651):0.71283):0.08752,('CCW36281.1_Biopolymer_transport_proteins_Chthono':0.451958,'AIE86461.1_MotA/TolQ/ExbB_proton_channel_Fimbrii':0.259927):1.17398):0.07224):0.034896):0.556216,'CAG37703.1_related_to_transport_protein_TolQ_Des':1.73869):0.007995,'ACX73656.1_MotA/TolQ/ExbB_proton_channel_Fibroba':1.94911):0.119976,((((((((((((((((((((((((((((('ACL94252.1_flagellar_motor_stator_protein_MotA_C':0.169924,'ADL02462.1_MotA/TolQ/ExbB_proton_channel_Brevund':0.181287):0.246726,('WP_002721884.1_flagellar_motor_stator_protein_Mot':0.160958,'ADE87210.1_chemotaxis_protein_MotA_Rhodobacter_c':0.21375):0.319058):0.060501,'AEI88772.1_flagellar_motor_protein_MotA_duplicate':0.75628):0.113335,'QNQ09341.1_flagellar_motor_stator_protein_MotA_S':0.739784):0.26845,('AIQ89863.1_Flagellar_motor_rotation_protein_MotA_':0.480121,'ABQ62189.1_chemotaxis_motA_protein_Brucella_ovis':0.312526):0.803508):0.089174,'AEI89251.1_flagellar_motor_protein_MotA_Candidat':1.43634):0.025678,((('CDK98316.1_putative_Motility_protein_A_Magnetosp':0.343079,'AEO48360.1_motA_chemotaxis_motility_protein_A_':0.25254):0.203904,'AIQ92722.1_Flagellar_motor_rotation_protein_MotA_':0.524532):0.070523,'AIK95837.1_flagellar_motor_protein_MotA_Candidat':0.472026):0.271682):0.170305,((((('APW43194.1_flagellar_motor_stator_protein_MotA_R':0.047468,'APW42908.1_flagellar_motor_stator_protein_MotA_R':0.159638):0.0579,'UBQ05764.1_flagellar_motor_stator_protein_MotA_T':0.170245):0.093718,'ATU66302.1_flagellar_motor_stator_protein_MotA_R':0.170562):0.190293,'QKS29226.1_MAG_flagellar_motor_stator_protein_Mo':0.452908):0.047405,((('AJJ32465.1_flagellar_motor_stator_protein_MotA_Y':0.003481,'AJJ54987.1_flagellar_motor_stator_protein_MotA_Y':2e-06):0.038983,'QXB11670.1_flagellar_motor_stator_protein_MotA_K':0.073836):0.037413,('QKY95842.1_flagellar_motor_stator_protein_MotA_S':0.006937,'AAC74960.1_motility_protein_A_Escherichia_coli_s':2e-06):0.173392):0.509113):0.183117):0.118009,((((('AJJ56783.1_flagellar_motor_stator_protein_MotA_Y':0.00363,'AJJ32757.1_flagellar_motor_stator_protein_MotA_Y':2e-06):0.482612,'ATU65955.1_flagellar_motor_stator_protein_MotA_R':0.507491):0.129971,'BAC62899.1_chemotaxis_LafT_protein_Vibrio_paraha':0.350955):0.732453,'AAG08339.1_chemotaxis_protein_MotA_Pseudomonas_a':0.506879):0.190922,'QEI18926.1_flagellar_motor_stator_protein_MotA_C':0.476319):0.109341):0.190383,(('ABC77078.1_chemotaxis_protein_Syntrophus_aciditr':0.396103,'BCB96094.1_flagellar_motor_stator_protein_MotA_D':0.18619):0.19866,'ACI21445.1_lateral_flagellar_motor_protein_MotA_':0.350845):0.297901):0.09125,(((('QOY88835.1_flagellar_motor_stator_protein_MotA_P':0.349304,'QOY86330.1_flagellar_motor_stator_protein_MotA_P':0.367043):0.159876,'ACO33337.1_putative_chemotaxis_MotA_protein_Acid':0.388978):0.135559,'QOY85080.1_flagellar_motor_stator_protein_MotA_P':0.545117):0.150379,'QQS05761.1_MAG_flagellar_motor_stator_protein_Mo':1.02595):0.125916):0.081958,(((('BAH37338.1_chemotaxis_MotA_protein_Gemmatimonas_':0.040541,'AMW06484.1_flagellar_motor_protein_MotA_Gemmatim':0.061673):0.035848,'QJR36484.1_flagellar_motor_stator_protein_MotA_G':0.078263):0.372576,'ARU42506.1_flagellar_motor_stator_protein_MotA_A':0.622046):0.031441,'CCW34551.1_flagellar_motor_stator_protein_MotA_C':0.497923):0.237245):0.096495,((('ATC63650.1_flagellar_motor_protein_MotA_Nibricoc':0.189675,'ACB73691.1_MotA/TolQ/ExbB_proton_channel_Opitutu':0.382897):0.586976,(PKPEBJJI_00564_Motility_protein_A_CAIZXV01_Verru:0.010191,DGBBEKCF_02484_Motility_protein_A_Verrucomicrobi:3e-06):1.01779):0.187351,'QDU25289.1_Chemotaxis_protein_LafT_Anatilimnocol':0.782278):0.323313):0.094938,'AHZ83632.1_flagellar_motor_protein_MotA_Bdellovi':0.813146):1.63367,(((('ACO34195.1_signal_recognition_particledocking_pro':1.86277,'CAB12348.1_twocomponent_sensor_histidine_kinase_':1.4607):0.975632,'QDU31046.1_hypothetical_protein_ETAA8_61990_Anat':1.92904):0.771159,('ADC89566.1_MotA/TolQ/ExbB_proton_channel_Thermoc':0.237579,'AHE95900.1_flagellar_motor_protein_MotA_Thermocr':0.105614):1.02851):3e-06,('AOO65115.1_flagellar_motor_rotation_protein_MotA_':0.424117,'CAL34488.1_putative_flagellar_motor_proton_channe':0.518654):0.673594):0.373441):0.224826,(((('CAB13242.1_motility_protein_A_MotA_component_of_':0.291524,'AJH79515.1_motA/TolQ/ExbB_proton_channel_family_p':0.307546):0.206139,'ATO49967.1_flagellar_motor_protein_MotA_Brevibac':0.369735):0.189505,'ALU34477.1_MotA/TolQ/ExbB_proton_channel_Clostri':0.636898):0.050893,'UPH48256.1_flagellar_motor_stator_protein_MotA_L':0.863926):0.303751):0.138273,(((((((('AFL88814.1_flagellar_motor_component_Terriglobus':2e-06,'AFL88473.1_flagellar_motor_component_Terriglobus':2e-06):0.543577,'ACO33277.1_chemotaxis_MotA_protein_Acidobacteriu':0.313056):0.2926,('QOY85398.1_flagellar_motor_protein_Paludibaculum':0.729573,'QOY88866.1_flagellar_motor_protein_Paludibaculum':0.717278):0.290422):0.343157,('AHZ85791.1_flagellar_motor_protein_MotA_Bdellovi':0.779096,'ACL65905.1_MotA/TolQ/ExbB_proton_channel_Anaerom':0.606688):0.31446):0.089279,('ALA58549.1_Flagellar_motor_protein_MotA_Nitrospi':0.483308,'BAM06016.1_flagellar_motor_component_Leptospiril':0.791432):0.298055):0.030024,(((('MBU2985819.1_flagellar_motor_protein_Saccharopha':0.005378,'ABD81419.1_MotA/TolQ/ExbB_proton_channel_Sacchar':1.999999999e-06):0.252391,'QEI19744.1_flagellar_motor_protein_Cellvibrio_ja':0.268776):0.197435,'AAG04849.1_MotC_Pseudomonas_aeruginosa_PAO1_':0.362147):0.43912,'QKS28544.1_MAG_flagellar_motor_protein_Candidat':0.605061):0.287984):0.0651,(('BCB96091.1_chemotaxis_protein_MotA_Dissulfurispi':0.321201,'ACI21187.1_chemotaxis_MotA_protein_Thermodesulfo':0.558445):0.543872,'AAR36419.1_flagellar_basal_body_stator_protein_Mo':0.513238):0.15004):0.207645,'QFG02534.1_flagellar_motor_protein_Tepidiforma_b':1.02362):0.094179):0.040343,(((((('QTD49844.1_hypothetical_protein_J3U87_30040_Sulf':4.03983,'CAB14951.1_sodium_channel_statorforce_generator_s':0.609304):0.405478,'QAT63225.1_motility_protein_A_Tissierella_sp._JN':1.14903):0.114918,'ATO49114.1_chemotaxis_protein_Brevibacillus_late':2.52648):0.325145,'EKT88755.1_endoflagellar_motor_protein_Leptospir':1.45764):0.070431,(('BAD41956.1_flagellar_motor_protein_MotA_Symbioba':0.978354,'APC08935.1_chemotaxis_protein_PomA_Moorella_ther':0.827816):0.050204,'ABY95376.1_MotA/TolQ/ExbB_proton_channel_Thermoa':0.770162):0.204934):3.000000001e-06,((('APC07539.1_chemotaxis_protein_PomA_Moorella_ther':0.684114,'ADL42743.1_MotA/TolQ/ExbB_proton_channel_Caldice':0.706771):0.075976,'ALU36624.1_Chemotaxis_protein_MotA_Clostridium_a':0.911061):0.094372,'ACL22206.1_MotA/TolQ/ExbB_proton_channel_Desulfi':0.612543):0.129638):0.052941):0.086448,'CCW35905.1_Flagellar_motor_component_Chthonomona':0.740042):0.270277,(((((((((((((((('QDU95463.1_MotA/TolQ/ExbB_proton_channel_family_p':0.592944,'QDU27457.1_MotA/TolQ/ExbB_proton_channel_family_p':0.562006):0.059306,'UAK41219.1_inorganic_phosphate_transporter_Bacte':3.20834):3e-06,'ADB18690.1_hypothetical_protein_Psta_4037_Pirell':0.473317):0.385573,'CAD71788.1_probable_biopolymer_transport_protein_':1.59075):0.487311,'ADY61197.1_MotA/TolQ/ExbB_proton_channel_Rubinis':1.34739):0.758582,'ADV63425.1_hypothetical_protein_Isop_2860_Isosph':2.09765):3e-06,(('QVL33268.1_MotA/TolQ/ExbB_proton_channel_family_p':1.61783,'ADB15801.1_hypothetical_protein_Psta_1118_Pirell':1.39619):0.094759,('QIZ73618.1_flagellar_motor_protein_MotA_Oxynema_':1.87797,'QDU28596.1_hypothetical_protein_ETAA8_36990_Anat':0.820655):0.351074):0.438565):0.661109,'QTD50492.1_MotA/TolQ/ExbB_proton_channel_family_p':1.9837):2e-06,'QOY89286.1_MotA/TolQ/ExbB_proton_channel_family_p':2.70424):2e-06,'QTD52802.1_MotA/TolQ/ExbB_proton_channel_family_p':1.07165):1.22058,(((((('ABA05996.1_conserved_hypothetical_protein_Nitrob':0.253511,'CDK97459.1_putative_MotA/TolQ/ExbB_proton_channel':0.245405):0.333144,'AEO47577.1_MotA/TolQ/ExbB_proton_channel_family_p':0.743491):0.223027,'ADE86713.1_motA/TolQ/ExbB_proton_channel_family_p':0.839052):0.165261,('AIQ93688.1_protein_of_unassigned_function_Methyl':0.878118,'BDB96353.1_flagellar_motor_protein_MotA_Candidat':2.67635):2e-06):2.84679,'QIZ70514.1_hypothetical_protein_HCG48_07910_Oxyn':3.11038):3e-06,'BAM05883.1_hypothetical_protein_LFE_0155_Leptosp':4.87833):1.5408):0.67488,(((((('AAG05312.1_hypothetical_protein_PA1924_Pseudomon':1.1829,'MBU2987451.1_MotA/TolQ/ExbB_proton_channel_family':0.455714):0.474875,'APW42573.1_biopolymer_transporter_ExbD_Rhodofera':0.501918):1.73314,('BCB95963.1_flagellar_motor_protein_MotA_Dissulfu':2.19326,'AVM44306.1_hypothetical_protein_C5Q97_06090_Vict':1.18785):0.111679):0.180194,'QVL34843.1_MotA/TolQ/ExbB_proton_channel_family_p':1.37357):0.267301,'ALU35076.1_MotA/TolQ/ExbB_proton_channel_Clostri':3.10351):0.018531,'ADY61698.1_MotA/TolQ/ExbB_proton_channel_Rubinis':1.82321):1.6867):1.57051,'AEO47025.1_MotA/TolQ/ExbB_proton_channel_Rhodosp':1.39368):0.125334,((('ABD82476.1_MotA/TolQ/ExbB_proton_channel_Sacchar':1.999999999e-06,'MBU2985330.1_flagellar_motor_protein_PomA_Saccha':1.999999999e-06):0.48309,('BAC58952.1_sodiumdriven_polar_flagellar_protein_M':0.030229,'QEO41881.1_flagellar_motor_protein_PomA_Vibrio_c':0.105294):0.476309):0.279058,('ABA04498.1_MotA/TolQ/ExbB_proton_channel_Nitroba':0.351606,'AIQ92767.1_MotA/TolQ/ExbB_proton_channel_Methylo':0.440989):0.938354):0.162831):0.035996,'WP_015920816.1_MotA/TolQ/ExbB_proton_channel_fami':1.30598):0.118059,(((('ADD67090.1_MotA/TolQ/ExbB_proton_channel_Denitro':0.299587,'QAR32385.1_motility_protein_A_Geovibrio_thiophil':0.267998):0.083063,'USF23676.1_Chemotaxis_protein_PomA_Mucispirillum':0.474213):0.166358,(('ALF47869.1_flagellar_motor_protein_Campylobacter':0.362931,'ADN09321.1_MotA/TolQ/ExbB_proton_channel_Sulfuri':0.291574):0.21575,'ADG94140.1_MotA/TolQ/ExbB_proton_channel_Arcobac':1.37653):0.27581):0.092589,'QTD53836.1_MotA/TolQ/ExbB_proton_channel_family_p':0.843872):0.187921):0.09542):0.050166,(('UOY13875.1_motility_protein_A_Treponema_pallidum':2e-06,'AWG41586.1_Mot_family_proton_H+_or_sodium_Na+_':2e-06):1.12573,'SCA58681.1_Motility_protein_A_Chlamydiales_bacte':0.979948):0.378638):0.029533,((('QDU26963.1_Chemotaxis_protein_PomA_Anatilimnocol':0.607933,'ADB18997.1_MotA/TolQ/ExbB_proton_channel_Pirellu':0.280091):0.23233,('QDU97593.1_Chemotaxis_protein_PomA_Lignipirellul':0.247831,'CAD76169.1_chemotaxis_pomA_protein_Rhodopirellul':0.432513):0.213067):0.141555,'ADY60068.1_MotA/TolQ/ExbB_proton_channel_Rubinis':0.540423):0.376295):0.076929,((('AEM21069.1_Flagellar_Motor_Protein_Brachyspira_i':0.044131,'CCG57546.1_flagellar_Motor_Protein_MotA_Brachysp':0.027625):0.804229,'EKT86693.1_flagellar_motor_protein_MotP_Leptospi':0.472122):0.312596,'APC07901.1_chemotaxis_protein_PomA_Moorella_ther':0.886524):0.15405):0.080529,((((('ATX78565.1_chemotaxis_protein_MotA_Mariprofundus':0.110438,'ATX82900.1_chemotaxis_protein_MotA_Mariprofundus':2e-06):0.582748,'CAG37394.1_related_to_flagellar_motor_apparatus_':0.603945):0.091094,'RLG45601.1_MAG_motility_protein_A_Candidatus_Ko':0.397182):0.23591,'QAT61614.1_motility_protein_A_Tissierella_sp._JN':0.882476):2e-06,('UPT73159.1_MAG_MotA/TolQ/ExbB_proton_channel_fam':1.16638,'RMG20999.1_MAG_motility_protein_A_Euryarchaeota':0.528733):2e-06):0.140441):3e-06,'AAC07083.1_flagellar_motor_protein_MotA_Aquifex_':1.27153):0.118342,'BAS27507.1_flagellar_motor_protein_MotP_Limnocho':0.781146):0.37118,(('ATX78592.1_chemotaxis_protein_MotA_Mariprofundus':0.037117,'ATX82928.1_chemotaxis_protein_MotA_Mariprofundus':2e-06):0.538677,'UBQ06690.1_MotA/TolQ/ExbB_proton_channel_family_p':0.909575):0.69947):0.046146,((('AEO48324.1_MotA/TolQ/ExbB_proton_channel_Rhodosp':0.355415,'CDK98228.1_putative_flagellar_motor_component_Ma':0.387703):0.499021,'QQS07318.1_MAG_MotA/TolQ/ExbB_proton_channel_fam':1.32647):0.863427,('AHZ83417.1_motility_protein_A_Bdellovibrio_bacte':1.6966,'QNQ12173.1_MotA/TolQ/ExbB_proton_channel_family_p':1.83174):0.106009):0.035266):2.507401,((('ADY59818.1_hypothetical_protein_Plabr_2216_Rubin':0.961368,'CAD77532.1_hypothetical_proteintransmembrane_pred':0.396126):4.2646,('ABD82073.1_MotA/TolQ/ExbB_proton_channel_Sacchar':2e-06,'MBU2984728.1_MotA/TolQ/ExbB_proton_channel_family':0.048154):2.32667):1.10582,'ADD68376.1_MotA/TolQ/ExbB_proton_channel_Denitro':2.62558):3e-06):0.133345):0.2420115,((((((((('ABL66229.1_hypothetical_protein_Cpha266_2232_Chl':2.42883,'BAH37227.1_hypothetical_membrane_protein_Gemmati':3.53482):0.543033,'QTD54135.1_hypothetical_protein_J3U87_16945_Sulf':1.19201):0.911492,'AEV97736.1_OmpA/MotB_domain_protein_Niastella_ko':0.539167):3e-06,'AHZ86390.1_hypothetical_protein_EP01_15820_Bdell':5.77972):3e-06,(((('ALF47932.1_putative_membrane_protein_Campylobact':1.33806,'CAL34744.1_putative_membrane_protein_Campylobact':0.182275):0.373583,'AOO65010.1_membrane_protein_Sulfurospirillum_hal':0.319174):0.182659,'AKF24681.1_hypothetical_protein_YH65_04225_Sulfu':1.33223):2.39547,'ACX76639.1_hypothetical_protein_Fisuc_3059_Fibro':1.97132):0.032292):2.166,'AHE96301.1_biopolymer_transporter_ExbB_Thermocri':2.52289):3e-06,'AIZ45138.1_hypothetical_protein_QR90_08530_Deino':5.49063):0.335235,('ACX76601.1_hypothetical_protein_Fisuc_3021_Fibro':5.17726,'AOO66357.1_hypothetical_protein_SHALO_2598_Sulfu':3.03838):1.67244):1.56302,'AUX27689.1_ABC_transporter_permease_Sorangium_ce':6.17168):0.2420115);"
	tr = ape::read.tree(file="", text=trstr)
	return(tr)
	}

# Check a tree for all common issues
checktree <- function(tr, with_explanations=FALSE)
	{
	junk='
	tr = motA_tree_example()
	
	checkdf = checktree(tr, with_explanations=TRUE)
	checkdf

	checkdf = checktree(tr, with_explanations=FALSE)
	checkdf
	
	' # END junk
	
	

	# Load the tree
	# 2017-08-04: Added check for NEXUS input tree
	# tr = read.tree(inputs$trfn
	#tr = check_trfn(trfn=inputs$trfn)
	
	good_bgb    = c(T, T, T, T, F, T, F, F, T, T, T, T)
	good_bioinf = c(T, T, T, T, F, T, F, F, T, T, T, T, T, T, T, T, T)
	
	checknames = c("Object exists?", 
	"Is it a 'phylo' object?",
	"Is it a rooted tree?",
	"Is it a binary (bifurcating) tree?",
	"Are there singleton nodes?",
	"Are the tip.labels unique?", 
	"Are there any negative branchlengths?",
	"Are there any zero (0.0) branchlengths?",
	"Are the tipnames free of spaces?", # end of BGB check
	"Are the tipnames free of single quotes (')?",
	"Are the tipnames free of single quotes (`)?",
	'Are the tipnames free of double quotes (")?',
	"Are the tipnames free of equals (=)?",
	"Are the tipnames free of parentheses ('(')?",
	"Are the tipnames free of parentheses (')')?",
	"Are the tipnames free of brackets ('[')?",
	"Are the tipnames free of brackets (']')?",
	"Passes all checks for BioGeoBEARS?",
	"Passes all checks for easy bioinformatics?")
	
	nodenums_to_check = rep("", times=length(checknames))
	
	checkTFs = rep(FALSE, times=length(checknames))
	checktxt = rep("", times=length(checknames))
	
	i = 0
	
	# Make sure it exists
	if (exists("tr") == FALSE)
		{
		checkTFs[(i=i+1)] = FALSE
		checktxt[i] = "Object 'tr' does not exist, according to exists(tr)."
		} else {
		checkTFs[(i=i+1)] = TRUE
		}
	
	if (("phylo" %in% class(tr)) == FALSE)
		{
		checkTFs[(i=i+1)] = FALSE
		checktxt[i] = "Object 'tr' does not have 'phylo' in class(tr)."

		} else {
		checkTFs[(i=i+1)] = TRUE
		}

	# Check for rooted tree
	if (is.rooted(tr) == FALSE)
		{
		stoptxt = paste("\nchecktree() says: Your tree is not rooted, i.e. is.rooted(tr) returns FALSE.\n", 
		"\nYou must fix the Newick file. APE's root() function is an option, or e.g. re-rooting by hand in FigTree.  However, you will want to make sure that all your tips still come up to the present (assuming you have a typical molecular tree, i.e. no fossils). \n", sep="")
		checkTFs[(i=i+1)] = FALSE
		checktxt[i] = stoptxt
		} else {
		checkTFs[(i=i+1)] = TRUE
		}



	trtable = prt(tr, printflag=FALSE)

	# Check for polytomies
	if (is.binary(tr) == FALSE)
		{
		stoptxt = paste("\nchecktree() says: Your tree not bifurcating, i.e. is.binary(tr) returns FALSE.\n", 
		"\nYou must fix the Newick file. APE's multi2di() function is an option.  See ?checktree for comments.\n", sep="")
		checkTFs[(i=i+1)] = FALSE
		checktxt[i] = stoptxt
		
		# Find the polytomy nodes
		tmp = sapply(X=trtable$daughter_nds, FUN=length)
		polytomy_TF = tmp > 2
		polytomy_nodes = (1:nrow(trtable))[polytomy_TF]
		polytomy_nodes_txt = paste0(polytomy_nodes, collapse=",", sep="")
		nodenums_to_check[i] = polytomy_nodes_txt
		} else {
		checkTFs[(i=i+1)] = TRUE
		}

	# Check for singletons
	if (has.singles(tr) == FALSE)
		{
		stoptxt = paste("\nchecktree() says: Your tree not bifurcating, because it as singletons (direct ancestor) nodes. I.e. has.singles(tr) returns FALSE.\n", 
		"\nYou must fix the Newick file. APE's collapse.singles() function is an option.  See ?checktree for comments.\n", sep="")
		checkTFs[(i=i+1)] = FALSE
		checktxt[i] = stoptxt


		# Find the singleton nodes
		tmp = sapply(X=trtable$daughter_nds, FUN=length)
		singleton_TF = tmp == 1
		singleton_nodes = (1:nrow(trtable))[singleton_TF]
		singleton_nodes_txt = paste0(singleton_nodes, collapse=",", sep="")
		nodenums_to_check[i] = singleton_nodes_txt


		} else {
		checkTFs[(i=i+1)] = TRUE
		}



	
	# Check that all tipnames are unique
	tipnames = tr$tip.label
	uniq_tipnames = unique(tipnames)
	if (length(uniq_tipnames) != length(tipnames))
		{
		stoptxt = paste("\n\nchecktree() says: your tree has non-unique tipnames. Make all tipnames unique in both your tree file and geography file.  Current tipnames are listed below:\n\n", sep="")
		
		checkTFs[(i=i+1)] = FALSE
		checktxt[i] = stoptxt
		
		namecounts = c(table(tipnames)) 
		non_uniq_TF = namecounts > 1
		non_uniq_names = names(namecounts)[non_uniq_TF]
		
		rownums_in_table = non_uniq_names %in% trtable$label
		nodes_txt = paste0(rownums_in_table, collapse=",", sep="")
		nodenums_to_check[i] = nodes_txt
		
		} else {
		checkTFs[(i=i+1)] = TRUE
		}
nodenums_to_check


	# Check for negative branchlengths
	brlen_equal_below_0_TF = tr$edge.length <= 0
	if (sum(brlen_equal_below_0_TF) > 0)
		{
		tr_table = prt(tr, printflag=FALSE)
		rows_w_BL0_TF = tr_table$edge.length <= 0
		rows_w_BL0_TF[is.na(rows_w_BL0_TF)] = FALSE
		
		nodenums = tr_table$node[rows_w_BL0_TF]
		branchlengths = tr_table$edge.length[rows_w_BL0_TF]
		edge_nums = tr_table$parent_br[rows_w_BL0_TF]
		tmptable = cbind(nodenums, branchlengths, edge_nums)
		tmptable = as.data.frame(tmptable, stringsAsFactors=FALSE)
		
		tmptxt = paste0(nodenums, collapse=",")
		tmptxt2 = paste0(edge_nums, collapse=",")
		
		stoptxt = paste("\nchecktree() says: FATAL ERROR in checktree(): the input tree has branchlengths <= 0, at some nodes edge numbers:\n\nThis can sometimes happen in e.g. MCC (majority clade consensus) trees output by BEAST's TreeAnnotator.\nYou must fix the Newick file. See ?checktree, and PhyloWiki, for comments. One option is to try impose_min_brlen() and then write the tree to a file.\n", sep="")
		checkTFs[(i=i+1)] = FALSE
		checktxt[i] = stoptxt

		nodes_txt = paste0(nodenums, collapse=",", sep="")
		nodenums_to_check[i] = nodes_txt

		} else {
		checkTFs[(i=i+1)] = TRUE
		}

	# Check for zero branchlengths
	brlen_equal_0_TF = tr$edge.length == 0
	if (sum(brlen_equal_below_0_TF) > 0)
		{
		tr_table = prt(tr, printflag=FALSE)
		rows_w_BL0_TF = tr_table$edge.length <= 0
		rows_w_BL0_TF[is.na(rows_w_BL0_TF)] = FALSE
		
		nodenums = tr_table$node[rows_w_BL0_TF]
		branchlengths = tr_table$edge.length[rows_w_BL0_TF]
		edge_nums = tr_table$parent_br[rows_w_BL0_TF]
		tmptable = cbind(nodenums, branchlengths, edge_nums)
		tmptable = as.data.frame(tmptable, stringsAsFactors=FALSE)
		
		tmptxt = paste0(nodenums, collapse=",")
		tmptxt2 = paste0(edge_nums, collapse=",")
		
		stoptxt = paste("\nchecktree() says: FATAL ERROR in checktree(): the input tree has branchlengths == 0, at some nodes edge numbers:\n\nThis can sometimes happen in e.g. MCC (majority clade consensus) trees output by BEAST's TreeAnnotator.\nYou must fix the Newick file. See ?checktree, and PhyloWiki, for comments. One option is to try impose_min_brlen() and then write the tree to a file.\n", sep="")
		checkTFs[(i=i+1)] = FALSE
		checktxt[i] = stoptxt

		nodes_txt = paste0(nodenums, collapse=",", sep="")
		nodenums_to_check[i] = nodes_txt

		} else {
		checkTFs[(i=i+1)] = TRUE
		}
	




	# Check that all tipnames have no spaces
	TF = grepl(" ", tipnames)
	if (sum(TF) > 0)
		{
		stoptxt = paste("\n\nchecktree() says: your tree has tipnames with spaces. Take these out of your Newick file and re-run. Tipnums are listed.\n\n", sep="")

		checkTFs[(i=i+1)] = FALSE
		checktxt[i] = stoptxt
		
		nodenums = (1:length(tipnames))[TF]
		nodes_txt = paste0(nodenums, collapse=",", sep="")
		nodenums_to_check[i] = nodes_txt
		} else {
		checkTFs[(i=i+1)] = TRUE
		}


	# Check that all tipnames have no single-quotes (')
	TF = grepl("'", tipnames)
	if (sum(TF) > 0)
		{
		stoptxt = paste("\n\nchecktree() says: your tree has tipnames with apostrophes ('). Take these out of your Newick file and re-run. Tipnums are listed.\n\n", sep="")
		
		checkTFs[(i=i+1)] = FALSE
		checktxt[i] = stoptxt

		nodenums = (1:length(tipnames))[TF]
		nodes_txt = paste0(nodenums, collapse=",", sep="")
		nodenums_to_check[i] = nodes_txt

		} else {
		checkTFs[(i=i+1)] = TRUE
		}

	# Check that all tipnames have no single-quotes (`)
	TF = grepl("`", tipnames)
	if (sum(TF) > 0)
		{
		stoptxt = paste("\n\nchecktree() says: your tree has tipnames with apostrophes (`). Take these out of your Newick file and re-run. Tipnums are listed.\n\n", sep="")
		
		checkTFs[(i=i+1)] = FALSE
		checktxt[i] = stoptxt

		nodenums = (1:length(tipnames))[TF]
		nodes_txt = paste0(nodenums, collapse=",", sep="")
		nodenums_to_check[i] = nodes_txt

		} else {
		checkTFs[(i=i+1)] = TRUE
		}

	# Check that all tipnames have no double-quotes ("")
	TF = grepl('"', tipnames)
	if (sum(TF) > 0)
		{
		stoptxt = paste('\n\nchecktree() says: your tree has tipnames with apostrophes ("). Take these out of your Newick file and re-run. Tipnums are listed.\n\n', sep="")
		
		checkTFs[(i=i+1)] = FALSE
		checktxt[i] = stoptxt

		nodenums = (1:length(tipnames))[TF]
		nodes_txt = paste0(nodenums, collapse=",", sep="")
		nodenums_to_check[i] = nodes_txt

		} else {
		checkTFs[(i=i+1)] = TRUE
		}

	# Check that all tipnames have no equals (=)
	TF = grepl("\\=", tipnames)
	if (sum(TF) > 0)
		{
		stoptxt = paste("\n\nchecktree() says: your tree has tipnames with equals signs ('='). Take these out of your Newick file and re-run. Tipnums are listed.\n\n", sep="")
		
		checkTFs[(i=i+1)] = FALSE
		checktxt[i] = stoptxt

		nodenums = (1:length(tipnames))[TF]
		nodes_txt = paste0(nodenums, collapse=",", sep="")
		nodenums_to_check[i] = nodes_txt

		} else {
		checkTFs[(i=i+1)] = TRUE
		}

	# Check that all tipnames have no brackets ("(")
	TF = grepl("\\(", tipnames)
	if (sum(TF) > 0)
		{
		stoptxt = paste("\n\nchecktree() says: your tree has tipnames with parentheses ('('). Take these out of your Newick file and re-run. Tipnums are listed.\n\n", sep="")
		
		checkTFs[(i=i+1)] = FALSE
		checktxt[i] = stoptxt

		nodenums = (1:length(tipnames))[TF]
		nodes_txt = paste0(nodenums, collapse=",", sep="")
		nodenums_to_check[i] = nodes_txt

		} else {
		checkTFs[(i=i+1)] = TRUE
		}

	TF = grepl("\\)", tipnames)
	if (sum(TF) > 0)
		{
		stoptxt = paste("\n\nchecktree() says: your tree has tipnames with parentheses (')'). Take these out of your Newick file and re-run. Tipnums are listed.\n\n", sep="")
		
		checkTFs[(i=i+1)] = FALSE
		checktxt[i] = stoptxt

		nodenums = (1:length(tipnames))[TF]
		nodes_txt = paste0(nodenums, collapse=",", sep="")
		nodenums_to_check[i] = nodes_txt

		} else {
		checkTFs[(i=i+1)] = TRUE
		}

	TF = grepl("\\[", tipnames)
	if (sum(TF) > 0)
		{
		stoptxt = paste("\n\nchecktree() says: your tree has tipnames with square brackets ('['). Take these out of your Newick file and re-run. Tipnums are listed.\n\n", sep="")
		
		checkTFs[(i=i+1)] = FALSE
		checktxt[i] = stoptxt

		nodenums = (1:length(tipnames))[TF]
		nodes_txt = paste0(nodenums, collapse=",", sep="")
		nodenums_to_check[i] = nodes_txt

		} else {
		checkTFs[(i=i+1)] = TRUE
		}


	TF = grepl("\\]", tipnames)
	if (sum(TF) > 0)
		{
		stoptxt = paste("\n\nchecktree() says: your tree has tipnames with square brackets (']'). Take these out of your Newick file and re-run. Tipnums are listed.\n\n", sep="")
		
		checkTFs[(i=i+1)] = FALSE
		checktxt[i] = stoptxt

		nodenums = (1:length(tipnames))[TF]
		nodes_txt = paste0(nodenums, collapse=",", sep="")
		nodenums_to_check[i] = nodes_txt

		} else {
		checkTFs[(i=i+1)] = TRUE
		}

	
	
	# BioGeoBEARS check
	TF1 = checkTFs[1:length(good_bgb)]
	if (all(TF1 == good_bgb) == FALSE)
		{
		stoptxt = paste("\n\nchecktree() says: Not all BioGeoBEARS tree checks passed (checkTFs 1-", length(good_bgb), ").", sep="")
		
		checkTFs[(i=i+1)] = FALSE
		checktxt[i] = stoptxt

		} else {
		checkTFs[(i=i+1)] = TRUE
		}
	
	TF1 = checkTFs[1:length(good_bioinf)]
	if (all(TF1 == good_bioinf) == FALSE)
		{
		stoptxt = paste("\n\nchecktree() says: Not all 'good checks for bioinformatics labels' passed (checkTFs 1-", length(good_bioinf), ").", sep="")
		
		checkTFs[(i=i+1)] = FALSE
		checktxt[i] = stoptxt

		} else {
		checkTFs[(i=i+1)] = TRUE
		}
	
	
	
	
	if (with_explanations == TRUE)
		{
		tmpmat = cbind(checknames, checkTFs, checktxt, nodenums_to_check)
		} else {
		tmpmat = cbind(checknames, checkTFs)
		}
	checkdf = as.data.frame(tmpmat, stringsAsFactors=FALSE)
	row.names=NULL
	checkdf
	
	return(checkdf)
	} # END checktree
	

