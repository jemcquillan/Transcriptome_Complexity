library("phytools")

setwd("/Volumes/GoogleDrive/My\ Drive/OrthoDB/Ortholog_Final_Files/Drosophila_Complexity/Drosophila_tree")

# Full dataset for Deuterostomes, the tree and the traits
drosophila_traits_check = read.csv("Drosophila_WG_gtf_species_full_list.csv", sep = ",")
drosophila_traits_WG = read.table("Drosophila_WG_gtf_species_full_list.csv", sep = ",", row.names=1, header = TRUE)
drosophila_traits_ortho = read.table("Drosophila_Ortholog_gtf_species_full_list.csv", sep = ",", row.names=1, header = TRUE)
drosophila_traits <- paste0("ortho_", colnames(drosophila_traits_ortho)), paste0("WG_", colnames(drosophila_traits_WG))
drosophila_traits_check$species

#Read trees and check tip labels
full_dros_tree = ape::read.tree(file="Drosophila_Species_Tree_List.tre")
plot(full_dros_tree)
full_dros_tree$tip.label

#Sanity check the names are the same between the traits and the tree
match(full_dros_tree$tip.label,drosophila_traits_check$species, nomatch = "not a match")
length(match(full_dros_tree$tip.label,drosophila_traits$species, nomatch = 0))

# Match tree tip labels to the trait lables, should be the same species names (if not edit tree to match gtf)
q = setNames(drosophila_traits$WG_mean_exonPerTranscript, rownames(drosophila_traits))
length(q)#Check the no tips dropped between the trait and tree files

# Likelihood test for rate variation in a continuous trait
brownian_chisqt <- brownie.lite(full_deut_tree, q, maxit=999999999, test="chisq", se=NULL)

brownian_sim <- brownie.lite(full_deut_tree, q, maxit=1000000, test="simulation", nsim=100, se=NULL)

v = setNames(drosophila_traits$ortho_mean_exonPerTranscript, rownames(drosophila_traits))
multiTree = rep(full_deut_tree,2)
inpt=list(q,v)
outpt=ratebytree(multiTree, inpt)
posthoc(outpt, p.adjs="BH")
