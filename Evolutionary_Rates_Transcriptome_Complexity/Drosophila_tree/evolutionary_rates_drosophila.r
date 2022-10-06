library("phytools")
library(dplyr)

setwd("/Volumes/GoogleDrive/My\ Drive/OrthoDB/Ortholog_Final_Files/Drosophila_Complexity/Drosophila_tree")

# Full dataset files for Drosophila, the tree and the traits, to prepare
drosophila_traits_WG = read.table("Drosophila_WG_gtf_species_full_list.csv", sep = ",", row.names=1, header = TRUE)
drosophila_traits_ortho = read.table("Drosophila_Ortholog_gtf_species_full_list.csv", sep = ",", row.names=1, header = TRUE)
colnames(drosophila_traits_WG) <- paste0("WG_", colnames(drosophila_traits_WG))
colnames(drosophila_traits_ortho) <- paste0("ortho_", colnames(drosophila_traits_ortho))

#Prepare the dataframes
drosophila_traits_WG_prepared <- drosophila_traits_WG %>% tibble::rownames_to_column('species') %>% 
  select(-c('WG_commonName', 'WG_group', 'WG_database', 'WG_taxid', 'WG_link', 'WG_filePrefix'))
drosophila_traits_ortho_prepared <- drosophila_traits_ortho %>% tibble::rownames_to_column('species') %>%
  select(-c('ortho_commonName', 'ortho_group', 'ortho_database', 'ortho_taxid', 'ortho_link', 'ortho_filePrefix'))

## Merge the orthologs and WG traits
drosophila_full_traits <- merge(drosophila_traits_WG_prepared,drosophila_traits_ortho_prepared,by = 'species',all = T)
row.names(drosophila_full_traits) <- drosophila_full_traits$species
write.csv(drosophila_full_traits, 'drosophila_full_traits.csv', row.names = FALSE)

#Read trees and check tip labels
full_dros_tree = ape::read.tree(file="Drosophila_Species_Tree_List.tre")
plot(full_dros_tree)
full_dros_tree$tip.label

#Sanity check the names are the same between the traits and the tree
match(full_dros_tree$tip.label,drosophila_full_traits$species, nomatch = "not a match")
length(match(full_dros_tree$tip.label,drosophila_full_traits$species, nomatch = 0))

# Match tree tip labels to the trait labels, should be the same species names (if not edit gtf to match tree)
print(drosophila_full_traits$WG_mean_exonPerTranscript)
print(rownames(drosophila_full_traits))
q = setNames(drosophila_full_traits$WG_mean_exonPerTranscript, rownames(drosophila_full_traits))
length(q)#Check the no tips dropped between the trait and tree files

WG_TpG = setNames(drosophila_full_traits$WG_mean_transcriptPerGene, rownames(drosophila_full_traits))
WG_EpT = setNames(drosophila_full_traits$WG_mean_exonPerTranscript, rownames(drosophila_full_traits))
WG_EpG = setNames(drosophila_full_traits$WG_mean_exonPerGene, rownames(drosophila_full_traits))


# Likelihood test for rate variation in a continuous trait
brownian_chisqt <- brownie.lite(full_dros_tree, q, maxit=1000, test="chisq", se=NULL)

brownian_sim <- brownie.lite(full_dros_tree, q, maxit=1000000, test="simulation", nsim=100, se=NULL)

ortho_TpG = setNames(drosophila_full_traits$ortho_mean_transcriptPerGene, rownames(drosophila_full_traits))
ortho_EpT = setNames(drosophila_full_traits$ortho_mean_exonPerTranscript, rownames(drosophila_full_traits))
ortho_EpG = setNames(drosophila_full_traits$ortho_mean_exonPerGene, rownames(drosophila_full_traits))
multiTree = rep(full_dros_tree,2)

TpG_inpt = list(WG_TpG,ortho_TpG)
EpT_inpt = list(WG_EpT,ortho_EpT)
EpG_inpt = list(WG_EpG,ortho_EpG)
TpG_outpt=ratebytree(multiTree, TpG_inpt)
EpT_outpt=ratebytree(multiTree, EpT_inpt)
EpG_outpt=ratebytree(multiTree, EpG_inpt)

posthoc(TpG_outpt)
posthoc(EpT_outpt)
posthoc(EpG_outpt)

contMap(full_dros_tree,EpT_outpt,outline=FALSE)
contMap(full_dros_tree,ortho_EpT,outline=FALSE)
