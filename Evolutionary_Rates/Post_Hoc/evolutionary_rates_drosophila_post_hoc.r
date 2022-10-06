library("phytools")
library(dplyr)

setwd("/Volumes/GoogleDrive/My\ Drive/OrthoDB/Ortholog_Final_Files/Drosophila_Complexity/Drosophila_tree")
# setwd("set your working directory here")

# Full dataset files for Drosophila, the tree and the traits, to prepare
drosophila_traits_WT = read.table("Drosophila_WT_gtf_species_full_list.csv", sep = ",", row.names=1, header = TRUE)
drosophila_traits_ortho = read.table("Drosophila_Ortholog_gtf_species_full_list.csv", sep = ",", row.names=1, header = TRUE)
colnames(drosophila_traits_WT) <- paste0("WT_", colnames(drosophila_traits_WT))
colnames(drosophila_traits_ortho) <- paste0("ortho_", colnames(drosophila_traits_ortho))

#Prepare the dataframes
drosophila_traits_WT_prepared <- drosophila_traits_WT %>% tibble::rownames_to_column('species') %>% 
  select(-c('WT_commonName', 'WT_group', 'WT_database', 'WT_taxid', 'WT_link', 'WT_filePrefix'))
drosophila_traits_ortho_prepared <- drosophila_traits_ortho %>% tibble::rownames_to_column('species') %>%
  select(-c('ortho_commonName', 'ortho_group', 'ortho_database', 'ortho_taxid', 'ortho_link', 'ortho_filePrefix'))

## Merge the orthologs and WT traits
drosophila_full_traits <- merge(drosophila_traits_WT_prepared,drosophila_traits_ortho_prepared,by = 'species',all = T)
row.names(drosophila_full_traits) <- drosophila_full_traits$species
write.csv(drosophila_full_traits, 'drosophila_full_traits.csv', row.names = FALSE)
# ^^^ Write out the traits into a file for future use as needed.

#Read trees and check tip labels
full_dros_tree = ape::read.tree(file="Drosophila_Species_Tree_List.tre")
plot(full_dros_tree)
full_dros_tree$tip.label
drosophila_full_traits$species

#Sanity check the names are the same between the traits and the tree
match(full_dros_tree$tip.label,drosophila_full_traits$species, nomatch = "not a match")
length(match(full_dros_tree$tip.label,drosophila_full_traits$species, nomatch = 0))

# Match tree tip labels to the trait labels, should be the same species names (if not edit gtf to match tree)
print(drosophila_full_traits$WT_mean_exonPerTranscript)
print(rownames(drosophila_full_traits))
q = setNames(drosophila_full_traits$WT_mean_exonPerTranscript, rownames(drosophila_full_traits))
length(q)#Check the no tips dropped between the trait and tree files

# Likelihood test for rate variation in a continuous trait
brownian_chisqt <- brownie.lite(full_dros_tree, q, maxit=1000, test="chisq", se=NULL)

brownian_sim <- brownie.lite(full_dros_tree, q, maxit=1000, test="simulation", nsim=100, se=NULL)

# Whole-transcriptome trait metric objects
dros_WT_TpG = setNames(drosophila_full_traits$WT_mean_transcriptPerGene, rownames(drosophila_full_traits))
dros_WT_EpT = setNames(drosophila_full_traits$WT_mean_exonPerTranscript, rownames(drosophila_full_traits))
dros_WT_EpG = setNames(drosophila_full_traits$WT_mean_exonPerGene, rownames(drosophila_full_traits))

# Ortholog trait metric objects
dros_ortho_TpG = setNames(drosophila_full_traits$ortho_mean_transcriptPerGene, rownames(drosophila_full_traits))
dros_ortho_EpT = setNames(drosophila_full_traits$ortho_mean_exonPerTranscript, rownames(drosophila_full_traits))
dros_ortho_EpG = setNames(drosophila_full_traits$ortho_mean_exonPerGene, rownames(drosophila_full_traits))

# Create the multitree ofr analysis between ortho vs WT
dros_multiTree = rep(full_dros_tree,2)

# Input for rate test along tree, furthermore post-hoc tests
dros_TpG_inpt = list(dros_WT_TpG,dros_ortho_TpG)
dros_EpT_inpt = list(dros_WT_EpT,dros_ortho_EpT)
dros_EpG_inpt = list(dros_WT_EpG,dros_ortho_EpG)

# Stored object for post-hoc tests
dros_TpG_outpt=ratebytree(dros_multiTree, dros_TpG_inpt)
dros_EpT_outpt=ratebytree(dros_multiTree, dros_EpT_inpt)
dros_EpG_outpt=ratebytree(dros_multiTree, dros_EpG_inpt)

# Post-hoc tests
posthoc(dros_TpG_outpt)
posthoc(dros_EpT_outpt)
posthoc(dros_EpG_outpt)

# Example visualizations of traits on the phylogenies
dros_Trait_Tree_ortho <- contMap(full_dros_tree,dros_ortho_EpT,fsize=c(0.6,1),outline=FALSE)
dros_Trait_Tree_WT <- contMap(full_dros_tree,dros_WT_EpT,fsize=c(0.6,1),outline=FALSE)
contMap(full_dros_tree,dros_ortho_EpT,fsize=c(0.6,1),outline=FALSE)

plotTree.wBars(full_dros_tree, dros_WT_EpT, fsize=0.7,type="fan")

plotTree.wBars(dros_Trait_Tree_WT$tree,dros_WT_EpT,method="plotSimmap",
               tip.labels=TRUE,fsize=0.7,colors=dros_Trait_Tree_WT$cols)
add.color.bar(1.0,dros_Trait_Tree_WT$cols,title="trait value",lims=dros_Trait_Tree_WT$lims,prompt=FALSE,
              x=0.9*par()$usr[1],y=0.9*par()$usr[3])
