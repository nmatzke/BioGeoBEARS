#######################################################
# Check Hannah Wood, Matzke et al. (2012) dated consensus tree
#######################################################

library(ape)
library(BioGeoBEARS)

wd = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples_new/Palp_M0_unconstrained/"
setwd(wd)

trfn = "tree.newick"
geogfn = "palp_no_Lacun.data"



# Read tree
tr = read.tree(trfn)
plot(tr)
axisPhylo()
title("original tree before pruning")

tipnames = tr$tip.label
length(tipnames)


# Read tipranges
tipranges = getranges_from_LagrangePHYLIP(geogfn)
tipranges_names = rownames(tipranges@df)


# Cut tips
tipnames_to_cut_TF = tipnames %in% tipranges_names == FALSE
tipnames_to_cut = tipnames[tipnames_to_cut_TF]
tipnames_to_cut


# Prune tree
tr2 = drop.tip(tr, tip=tipnames_to_cut)

tr2
plot(tr2)
axisPhylo()
title("tree after pruning")

outfn = "tree_pruned.newick"
write.tree(tr2, file=outfn)

