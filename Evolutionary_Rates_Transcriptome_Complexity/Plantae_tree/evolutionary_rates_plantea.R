library("phytools")
library(dplyr)

setwd("/Volumes/GoogleDrive/My\ Drive/OrthoDB/Ortholog_Final_Files/Plantae_Complexity/Plantae_tree")

# Full dataset for plantaes, the tree and the traits
plantae_traits_WG = read.csv("plantae_WG_gtf_species_full_list.csv", sep = ",", row.names=1, header = TRUE)
plantae_traits_ortho = read.csv("plantae_Ortholog_gtf_species_full_list.csv", sep = ",", row.names=1, header = TRUE)
colnames(plantae_traits_WG) <- paste0("WG_", colnames(plantae_traits_WG))
colnames(plantae_traits_ortho) <- paste0("ortho_", colnames(plantae_traits_ortho))

#Prepare the dataframes
plantae_traits_WG_prepared <- plantae_traits_WG %>% tibble::rownames_to_column('species') %>% 
  select(-c('WG_commonName', 'WG_group', 'WG_database', 'WG_taxid', 'WG_link', 'WG_filePrefix'))
plantae_traits_ortho_prepared <- plantae_traits_ortho %>% tibble::rownames_to_column('species') %>%
  select(-c('ortho_commonName', 'ortho_group', 'ortho_database', 'ortho_taxid', 'ortho_link', 'ortho_filePrefix'))

## Merge the orthologs and WG traits
plantae_full_traits <- merge(plantae_traits_WG_prepared,plantae_traits_ortho_prepared,by = 'species',all = T)
row.names(plantae_full_traits) <- plantae_full_traits$species
write.csv(plantae_full_traits, 'plantae_full_traits.csv', row.names = FALSE)

#Read trees and check tip labels
full_plant_tree = ape::read.tree(file="Plantae_Species_Tree_list.tre")
plot(full_plant_tree)
full_plant_tree$tip.label
plantae_full_traits$species

#Sanity check the names are the same between the traits and the tree
match(full_plant_tree$tip.label,plantae_full_traits$species, nomatch = "not a match")
length(match(full_plant_tree$tip.label,plantae_full_traits$species, nomatch = 0))

# Match tree tip labels to the trait lables, should be the same species names (if not edit tree to match gtf)
print(plantae_full_traits$WG_mean_exonPerTranscript)
print(rownames(plantae_full_traits))
WG_TpG = setNames(plantae_full_traits$WG_mean_transcriptPerGene, rownames(plantae_full_traits))
WG_EpT = setNames(plantae_full_traits$WG_mean_exonPerTranscript, rownames(plantae_full_traits))
WG_EpG = setNames(plantae_full_traits$WG_mean_exonPerGene, rownames(plantae_full_traits))
length(q)#Check the no tips dropped between the trait and tree files

# Likelihood test for rate variation in a continuous trait
brownian_chisqt <- brownie.lite(full_plant_tree, q, maxit=999999999, test="chisq", se=NULL)

brownian_sim <- brownie.lite(full_plant_tree, q, maxit=1000000, test="simulation", nsim=100, se=NULL)

ortho_TpG = setNames(plantae_full_traits$ortho_mean_transcriptPerGene, rownames(plantae_full_traits))
ortho_EpT = setNames(plantae_full_traits$ortho_mean_exonPerTranscript, rownames(plantae_full_traits))
ortho_EpG = setNames(plantae_full_traits$ortho_mean_exonPerGene, rownames(plantae_full_traits))
multiTree = rep(full_plant_tree,2)
TpG_inpt = list(WG_TpG,ortho_TpG)
EpT_inpt = list(WG_EpT,ortho_EpT)
EpG_inpt = list(WG_EpG,ortho_EpG)
TpG_outpt=ratebytree(multiTree, TpG_inpt)
EpT_outpt=ratebytree(multiTree, EpT_inpt)
EpG_outpt=ratebytree(multiTree, EpG_inpt)
posthoc(TpG_outpt)
posthoc(EpT_outpt)
posthoc(EpG_outpt)


# Only reading it prepared files
#plantae_traits_check = read.csv("plantae_full_traits.csv", sep = ",")
#plantae_traits = read.table("plantae_full_traits.csv", sep = ",", row.names=1, header = TRUE)
#plantae_traits_check$species
