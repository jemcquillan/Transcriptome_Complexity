library("phytools")
library(dplyr)

setwd("/Volumes/GoogleDrive/My Drive/OrthoDB/Ortholog_Final_Files/Fungi_Complexity/Fungi_tree")

# Full dataset for fungi, the tree and the traits
fungi_full_traits <- read.csv("Fungi_gtf_classes_full_list_Evo_Rates.csv")
row.names(fungi_full_traits) <- fungi_full_traits$group
# ^^^^ This need to be done since when reading row names, the default will be numeric for the row.
# So we rename them to the names that need to match up to the tip.labels. Without
#  R will not mach to row names to the colmn that we designate, here fungi_full_traits$group.

#Read trees and check tip labels
full_fungi_tree = ape::read.tree(file="fungi_class_pruned_no_polytomy.tre")
plot(full_fungi_tree)
full_fungi_tree$tip.label
fungi_full_traits$group

#Sanity check the names are the same between the traits and the tree
match(full_fungi_tree$tip.label,fungi_full_traits$group, nomatch = "not a match")
length(match(full_fungi_tree$tip.label,fungi_full_traits$group, nomatch = 0))

# Match tree tip labels to the trait lables, should be the same species names (if not edit tree to match gtf)
print(fungi_full_traits$WT_groupEpT)
print(rownames(fungi_full_traits))
q = setNames(fungi_full_traits$WT_groupEpT, rownames(fungi_full_traits))
length(q)#Check the no tips dropped between the trait and tree files

# Likelihood test for rate variation in a continuous trait
brownian_chisqt <- brownie.lite(full_fungi_tree, q, maxit=1000, test="chisq", se=NULL)

brownian_sim <- brownie.lite(full_fungi_tree, q, maxit=1000, test="simulation", nsim=100, se=NULL)

# Whole-transcriptome trait metric objects
fungus_WT_TpG = setNames(fungi_full_traits$WT_groupTpG, rownames(fungi_full_traits))
fungus_WT_EpT = setNames(fungi_full_traits$WT_groupEpT, rownames(fungi_full_traits))
fungus_WT_EpG = setNames(fungi_full_traits$WT_groupEpG, rownames(fungi_full_traits))

# Ortholog trait metric objects
fungus_ortho_TpG = setNames(fungi_full_traits$ortho_groupTpG, rownames(fungi_full_traits))
fungus_ortho_EpT = setNames(fungi_full_traits$ortho_groupEpT, rownames(fungi_full_traits))
fungus_ortho_EpG = setNames(fungi_full_traits$ortho_groupEpG, rownames(fungi_full_traits))

# Create the multitree ofr analysis between ortho vs WT
fungus_multiTree = rep(full_fungi_tree,2)

# Input for rate test along tree, furthermore post-hoc tests
fungus_TpG_inpt = list(fungus_WT_TpG, fungus_ortho_TpG)
fungus_EpT_inpt = list(fungus_WT_EpT, fungus_ortho_EpT)
fungus_EpG_inpt = list(fungus_WT_EpG, fungus_ortho_EpG)

# Stored object for post-hoc tests
fungus_TpG_outpt=ratebytree(fungus_multiTree, fungus_TpG_inpt)
fungus_EpT_outpt=ratebytree(fungus_multiTree, fungus_EpT_inpt)
fungus_EpG_outpt=ratebytree(fungus_multiTree, fungus_EpG_inpt)

# Post-hoc tests
posthoc(fungus_TpG_outpt)
posthoc(fungus_EpT_outpt)
posthoc(fungus_EpG_outpt)

# Example visualizations of traits on the phylogenies
fungus_Trait_Tree_ortho <- contMap(full_fungi_tree,fungus_ortho_EpT,fsize=c(0.6,1),outline=FALSE)
fungus_Trait_Tree_WT <- contMap(full_fungi_tree,fungus_WT_EpT,fsize=c(0.6,1),outline=FALSE)

plotTree.wBars(fungus_Trait_Tree_WT$tree,fungus_WT_EpT,method="plotSimmap",
               tip.labels=TRUE,fsize=0.7,colors=fungus_Trait_Tree_WT$cols)
add.color.bar(1.0,fungus_Trait_Tree_WT$cols,title="trait value",lims=fungus_Trait_Tree_WT$lims,prompt=FALSE,
              x=0.9*par()$usr[1],y=0.9*par()$usr[3])

