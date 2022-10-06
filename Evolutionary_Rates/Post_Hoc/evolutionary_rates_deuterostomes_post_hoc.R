library("phytools")
library(dplyr)

setwd("/Volumes/GoogleDrive/My\ Drive/OrthoDB/Ortholog_Final_Files/Deuterostome_Complexity/Deuterostome_tree")
# setwd("set your working directory here")

# Full dataset for Deuterostomes, the tree and the traits
deuterostome_traits_WT = read.csv("Deuterostome_WT_gtf_species_full_list.csv", sep = ",", row.names=1, header = TRUE)
deuterostome_traits_ortho = read.csv("Deuterostome_Ortholog_gtf_species_full_list.csv", sep = ",", row.names=1, header = TRUE)
colnames(deuterostome_traits_WT) <- paste0("WT_", colnames(deuterostome_traits_WT))
colnames(deuterostome_traits_ortho) <- paste0("ortho_", colnames(deuterostome_traits_ortho))

#Prepare the dataframes
deuterostome_traits_WT_prepared <- deuterostome_traits_WT %>% tibble::rownames_to_column('species') %>% 
  select(-c('WT_commonName', 'WT_group', 'WT_database', 'WT_taxid', 'WT_link', 'WT_filePrefix'))
deuterostome_traits_ortho_prepared <- deuterostome_traits_ortho %>% tibble::rownames_to_column('species') %>%
  select(-c('ortho_commonName', 'ortho_group', 'ortho_database', 'ortho_taxid', 'ortho_link', 'ortho_filePrefix'))

## Merge the orthologs and WT traits
deuterostome_full_traits <- merge(deuterostome_traits_WT_prepared,deuterostome_traits_ortho_prepared,by = 'species', all = T)
row.names(deuterostome_full_traits) <- deuterostome_full_traits$species
write.csv(deuterostome_full_traits, 'deuterostome_full_traits.csv', row.names = FALSE)
# ^^^ Write out the traits into a file for future use as needed.

#Read trees and check tip labels
full_deut_tree = ape::read.tree(file="Deuterostome_Species_TimeTree.tre")
plot(full_deut_tree)
full_deut_tree$tip.label
deuterostome_full_traits$species

#Sanity check the names are the same between the traits and the tree
match(full_deut_tree$tip.label,deuterostome_full_traits$species, nomatch = "not a match")
length(match(full_deut_tree$tip.label,deuterostome_full_traits$species, nomatch = 0))

# Match tree tip labels to the trait lables, should be the same species names (if not edit tree to match gtf)
print(deuterostome_full_traits$WT_mean_exonPerTranscript)
print(rownames(deuterostome_full_traits))
q <- setNames(deuterostome_full_traits$WT_mean_exonPerTranscript, rownames(deuterostome_full_traits))
length(q)#Check the no tips dropped between the trait and tree files

# Likelihood test for rate variation in a continuous trait
brownian_chisqt <- brownie.lite(full_deut_tree, q, maxit=1000, test="chisq", se=NULL)

brownian_sim <- brownie.lite(full_deut_tree, q, maxit=1000, test="simulation", nsim=100, se=NULL)

# Whole-transcriptome trait metric objects
deut_WT_TpG = setNames(deuterostome_full_traits$WT_mean_transcriptPerGene, rownames(deuterostome_full_traits))
deut_WT_EpT = setNames(deuterostome_full_traits$WT_mean_exonPerTranscript, rownames(deuterostome_full_traits))
deut_WT_EpG = setNames(deuterostome_full_traits$WT_mean_exonPerGene, rownames(deuterostome_full_traits))

# Ortholog trait metric objects
deut_ortho_TpG = setNames(deuterostome_full_traits$ortho_mean_transcriptPerGene, rownames(deuterostome_full_traits))
deut_ortho_EpT = setNames(deuterostome_full_traits$ortho_mean_exonPerTranscript, rownames(deuterostome_full_traits))
deut_ortho_EpG = setNames(deuterostome_full_traits$ortho_mean_exonPerGene, rownames(deuterostome_full_traits))

# Create the multitree ofr analysis between ortho vs WT
deut_multiTree = rep(full_deut_tree,2)

# Input for rate test along tree, furthermore post-hoc tests
deut_TpG_inpt = list(deut_WT_TpG, deut_ortho_TpG)
deut_EpT_inpt = list(deut_WT_EpT, deut_ortho_EpT)
deut_EpG_inpt = list(deut_WT_EpG, deut_ortho_EpG)

# Stored object for post-hoc tests
deut_TpG_outpt=ratebytree(deut_multiTree, deut_TpG_inpt)
deut_EpT_outpt=ratebytree(deut_multiTree, deut_EpT_inpt)
deut_EpG_outpt=ratebytree(deut_multiTree, deut_EpG_inpt)

# Post-hoc tests
posthoc(deut_TpG_outpt)
posthoc(deut_EpT_outpt)
posthoc(deut_EpG_outpt)

# Example visualizations of traits on the phylogenies
deut_Trait_Tree_ortho <- contMap(full_deut_tree,deut_ortho_EpT,fsize=c(0.6,1),outline=FALSE)
deut_Trait_Tree_WT <- contMap(full_deut_tree,deut_WT_EpT,fsize=c(0.6,1),outline=FALSE)
contMap(full_deut_tree,deut_ortho_EpT,fsize=c(0.6,1),outline=FALSE)

plotTree.wBars(full_deut_tree, deut_WT_EpT, fsize=0.7,type="fan")

plotTree.wBars(deut_Trait_Tree_WT$tree,deut_WT_EpT,method="plotSimmap",
               tip.labels=TRUE,fsize=0.7,colors=deut_Trait_Tree_WT$cols)
add.color.bar(1.0,deut_Trait_Tree_WT$cols,title="trait value",lims=deut_Trait_Tree_WT$lims,prompt=FALSE,
              x=0.9*par()$usr[1],y=0.9*par()$usr[3])
