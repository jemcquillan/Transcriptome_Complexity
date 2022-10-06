library("phytools")
library(dplyr)

setwd("/Volumes/GoogleDrive/My\ Drive/OrthoDB/Ortholog_Final_Files/Plantae_Complexity/Plantae_tree")
# setwd("set your working directory here")

# Full dataset for plantaes, the tree and the traits
plantae_traits_WT = read.csv("plantae_WT_gtf_species_full_list.csv", sep = ",", row.names=1, header = TRUE)
plantae_traits_ortho = read.csv("plantae_Ortholog_gtf_species_full_list.csv", sep = ",", row.names=1, header = TRUE)
colnames(plantae_traits_WT) <- paste0("WT_", colnames(plantae_traits_WT))
colnames(plantae_traits_ortho) <- paste0("ortho_", colnames(plantae_traits_ortho))

#Prepare the dataframes
plantae_traits_WT_prepared <- plantae_traits_WT %>% tibble::rownames_to_column('species') %>% 
  select(-c('WT_commonName', 'WT_group', 'WT_database', 'WT_taxid', 'WT_link', 'WT_filePrefix'))
plantae_traits_ortho_prepared <- plantae_traits_ortho %>% tibble::rownames_to_column('species') %>%
  select(-c('ortho_commonName', 'ortho_group', 'ortho_database', 'ortho_taxid', 'ortho_link', 'ortho_filePrefix'))

## Merge the orthologs and WT traits
plantae_full_traits <- merge(plantae_traits_WT_prepared,plantae_traits_ortho_prepared,by = 'species',all = T)
row.names(plantae_full_traits) <- plantae_full_traits$species
write.csv(plantae_full_traits, 'plantae_full_traits.csv', row.names = FALSE)
# ^^^ Write out the traits into a file for future use as needed.

#Read trees and check tip labels
full_plant_tree = ape::read.tree(file="Plantae_Species_Tree_list.tre")
plot(full_plant_tree)
full_plant_tree$tip.label
plantae_full_traits$species

#Sanity check the names are the same between the traits and the tree
match(full_plant_tree$tip.label,plantae_full_traits$species, nomatch = "not a match")
length(match(full_plant_tree$tip.label,plantae_full_traits$species, nomatch = 0))

# Match tree tip labels to the trait lables, should be the same species names (if not edit tree to match gtf)
print(plantae_full_traits$WT_mean_exonPerTranscript)
print(rownames(plantae_full_traits))
q = setNames(plantae_full_traits$WT_mean_exonPerTranscript, rownames(plantae_full_traits))
length(q)#Check the no tips dropped between the trait and tree files

# Likelihood test for rate variation in a continuous trait
brownian_chisqt <- brownie.lite(full_plant_tree, q, maxit=1000, test="chisq", se=NULL)

brownian_sim <- brownie.lite(full_plant_tree, q, maxit=1000, test="simulation", nsim=100, se=NULL)

# Whole-transcriptome trait metric objects
plant_WT_TpG = setNames(plantae_full_traits$WT_mean_transcriptPerGene, rownames(plantae_full_traits))
plant_WT_EpT = setNames(plantae_full_traits$WT_mean_exonPerTranscript, rownames(plantae_full_traits))
plant_WT_EpG = setNames(plantae_full_traits$WT_mean_exonPerGene, rownames(plantae_full_traits))

# Ortholog trait metric objects
plant_ortho_TpG = setNames(plantae_full_traits$ortho_mean_transcriptPerGene, rownames(plantae_full_traits))
plant_ortho_EpT = setNames(plantae_full_traits$ortho_mean_exonPerTranscript, rownames(plantae_full_traits))
plant_ortho_EpG = setNames(plantae_full_traits$ortho_mean_exonPerGene, rownames(plantae_full_traits))

# Create the multitree ofr analysis between ortho vs WT
plant_multiTree = rep(full_plant_tree,2)

# Input for rate test along tree, furthermore post-hoc tests
plant_TpG_inpt = list(plant_WT_TpG,plant_ortho_TpG)
plant_EpT_inpt = list(plant_WT_EpT,plant_ortho_EpT)
plant_EpG_inpt = list(plant_WT_EpG,plant_ortho_EpG)

# Stored object for post-hoc tests
plant_TpG_outpt=ratebytree(plant_multiTree, plant_TpG_inpt)
plant_EpT_outpt=ratebytree(plant_multiTree, plant_EpT_inpt)
plant_EpG_outpt=ratebytree(plant_multiTree, plant_EpG_inpt)

# Post-hoc tests
posthoc(plant_TpG_outpt)
posthoc(plant_EpT_outpt)
posthoc(plant_EpG_outpt)

# Example visualizations of traits on the phylogenies
plant_Trait_Tree_ortho <- contMap(full_plant_tree,plant_ortho_EpT,fsize=c(0.6,1),outline=FALSE)
plant_Trait_Tree_WT <- contMap(full_plant_tree,plant_WT_EpT,fsize=c(0.6,1),outline=FALSE)
contMap(full_plant_tree,plant_ortho_EpT,fsize=c(0.6,1),outline=FALSE)

plotTree.wBars(full_plant_tree, plant_WT_EpT, fsize=0.7,type="fan")

plotTree.wBars(plant_Trait_Tree_WT$tree,plant_WT_EpT,method="plotSimmap",
               tip.labels=TRUE,fsize=0.7,colors=plant_Trait_Tree_WT$cols)
add.color.bar(1.0,plant_Trait_Tree_WT$cols,title="trait value",lims=plant_Trait_Tree_WT$lims,prompt=FALSE,
              x=0.9*par()$usr[1],y=0.9*par()$usr[3])
