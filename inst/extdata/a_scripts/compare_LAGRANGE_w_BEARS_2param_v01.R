#######################################################
# Compare LAGRANGE-Python, LAGRANGE-C++, and BioGeoBEARS runs
# for a particular directory
#######################################################


#######################################################
# Loading the packages
#######################################################
# After installing these packages, library() them
library(ape)		# R phylogenetics package
library(rexpokit)	# Matrix exponentiation package (by Matzke) using the FORTRAN library EXPOKIT
					# (by Robert Sidje) to do rapid matrix exponentiation; may require that you
					# have the FORTRAN compiler, gfortran, installed.
library(cladoRcpp)	# Package for rapid calculation of cladogenesis models (by Matzke) in C++ and Rcpp. 
library(BioGeoBEARS)# Package using rexpokit & cladoRcpp to do model testing and ancestral range estimation
					# for historical biogeography (by Matzke)

#######################################################
# Source R scripts (until they go into BioGeoBEARS)
#######################################################
sourcedir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/"
sourceall(path=sourcedir)


#######################################################
# Set the directory -- the directories/files should follow the standard naming
#######################################################
# Normally, with the installed package
extdata_dir = system.file("extdata", package="BioGeoBEARS")
# tmp hard code NICK FIX: 
extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"
setwd(extdata_dir)


#######################################################
# Load tree and geog files
#######################################################
datadir = "examples/Psychotria_M0/LGpy/"
wd = slashslash(paste(extdata_dir, datadir, sep="/"))
setwd(wd)
trfn = "treefile.newick"
geogfn = "treegeog.data"
tr = read.tree(file=trfn)
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
areanames = getareas_from_tipranges_object(tipranges=tipranges)
tipranges = order_tipranges_by_tree_tips(tipranges, tr)
tiprange_names = tipranges_to_area_strings(tipranges=tipranges)
tiprange_names


#######################################################
# Python LAGRANGE -- parse results
#######################################################
datadir = "examples/Psychotria_M0/LGpy/"
wd = slashslash(paste(extdata_dir, datadir, sep="/"))
setwd(wd)
sumstats = parse_lagrange_python_output(outfn="LGpy_results.txt")
splits_table_fn = sumstats$splits_table_fn

splits = LGpy_splits_fn_to_table(splits_fn=splits_table_fn)
splits

MLsplits_LGpy = LGpy_MLsplit_per_node(splits=splits)
MLsplits_LGpy

MLsplits_LGpy = map_LG_MLsplits_to_tree(MLsplits_LGpy, tr=tr, tiprange_names=tiprange_names, type="python")
MLsplits_LGpy

# Plot on the corners, with nice auto-colors
MLsplits = map_LG_MLsplits_to_tree_corners(MLsplits=MLsplits_LGpy, tr, tipranges, type="python")
title("ML splits, local optimization\n(Python LAGRANGE)")

#######################################################
# C++ LAGRANGE -- parse results
#######################################################
datadir = "examples/Psychotria_M0/LGcpp/"
wd = slashslash(paste(extdata_dir, datadir, sep="/"))
setwd(wd)
sumstats = parse_lagrange_output(outfn="LGcpp_results.txt")
splits_table_fn = sumstats$splits_table_fn

splits = LGcpp_splits_fn_to_table(splits_fn=splits_table_fn)
splits

MLsplits_LGcpp = LGpy_MLsplit_per_node(splits=splits)
MLsplits_LGcpp

MLsplits_LGcpp = map_LG_MLsplits_to_tree(MLsplits_LGcpp, tr=tr, tiprange_names=tiprange_names, type="C++")
MLsplits_LGcpp

# Plot on the corners, with nice auto-colors
MLsplits = map_LG_MLsplits_to_tree_corners(MLsplits=MLsplits_LGcpp, tr, tipranges, removechar="_", type="C++")
title("ML splits, local optimization\n(C++ LAGRANGE)")


#######################################################
# Add the states
#######################################################
states_table_fn = sumstats$states_table_fn

states = LGcpp_states_fn_to_table(states_fn=states_table_fn)
states

MLstates_LGcpp = LGcpp_MLstate_per_node(states=states)
MLstates_LGcpp

# Plot on the corners, with nice auto-colors
MLstates = map_LG_MLstates_to_tree(MLstates=MLstates_LGcpp, tr, tipranges, removechar="_", type="C++", newplot=FALSE)
#title("ML states, local optimization\n(C++ LAGRANGE)")


