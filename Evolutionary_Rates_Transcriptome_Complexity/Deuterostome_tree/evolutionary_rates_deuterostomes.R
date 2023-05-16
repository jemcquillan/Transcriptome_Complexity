library("phytools")
library(dplyr)

setwd("/Volumes/GoogleDrive/My\ Drive/OrthoDB/Ortholog_Final_Files/Deuterostome_Complexity/Deuterostome_tree")

# Full dataset for Deuterostomes, the tree and the traits
deuterostome_traits_WT = read.csv("Deuterostome_WT_gtf_species_full_list.csv", sep = ",", row.names=1, header = TRUE)
deuterostome_traits_ortho = read.csv("Deuterostome_Ortholog_gtf_species_full_list.csv", sep = ",", row.names=1, header = TRUE)
colnames(deuterostome_traits_WT) <- paste0("WG_", colnames(deuterostome_traits_WG))
colnames(deuterostome_traits_ortho) <- paste0("ortho_", colnames(deuterostome_traits_ortho))

#Prepare the dataframes
deuterostome_traits_WT_prepared <- deuterostome_traits_WG %>% tibble::rownames_to_column('species') %>% 
  select(-c('WG_commonName', 'WG_group', 'WG_database', 'WG_taxid', 'WG_link', 'WG_filePrefix'))
deuterostome_traits_ortho_prepared <- deuterostome_traits_ortho %>% tibble::rownames_to_column('species') %>%
  select(-c('ortho_commonName', 'ortho_group', 'ortho_database', 'ortho_taxid', 'ortho_link', 'ortho_filePrefix'))

## Merge the orthologs and WG traits
deuterostome_full_traits <- merge(deuterostome_traits_WT_prepared,deuterostome_traits_ortho_prepared,by = 'species',all = T)
row.names(deuterostome_full_traits) <- deuterostome_full_traits$species
write.csv(deuterostome_full_traits, 'deuterostome_full_traits.csv', row.names = FALSE)

#Read trees and check tip labels
full_deut_tree = ape::read.tree(file="Deuterostome_Species_TimeTree.tre")
plot(full_deut_tree)
full_deut_tree$tip.label
deuterostome_full_traits$species

#Sanity check the names are the same between the traits and the tree
match(full_deut_tree$tip.label,deuterostome_full_traits$species, nomatch = "not a match")
length(match(full_deut_tree$tip.label,deuterostome_full_traits$species, nomatch = 0))

# Match tree tip labels to the trait lables, should be the same species names (if not edit tree to match gtf)
print(deuterostome_full_traits$WG_mean_exonPerTranscript)
print(rownames(deuterostome_full_traits))

WT_TpG = setNames(deuterostome_full_traits$WG_mean_transcriptPerGene, rownames(deuterostome_full_traits))
WT_EpT = setNames(deuterostome_full_traits$WG_mean_exonPerTranscript, rownames(deuterostome_full_traits))
WT_EpG = setNames(deuterostome_full_traits$WG_mean_exonPerGene, rownames(deuterostome_full_traits))
length(q)#Check the no tips dropped between the trait and tree files

# Likelihood test for rate variation in a continuous trait
brownian_chisqt <- brownie.lite(full_deut_tree, q, maxit=99999999, test="chisq", se=NULL)

brownian_sim <- brownie.lite(full_deut_tree, q, maxit=1000000, test="simulation", nsim=100, se=NULL)

ortho_TpG = setNames(deuterostome_full_traits$ortho_mean_transcriptPerGene, rownames(deuterostome_full_traits))
ortho_EpT = setNames(deuterostome_full_traits$ortho_mean_exonPerTranscript, rownames(deuterostome_full_traits))
ortho_EpG = setNames(deuterostome_full_traits$ortho_mean_exonPerGene, rownames(deuterostome_full_traits))
multiTree = rep(full_deut_tree,2)
TpG_inpt = list(WT_TpG,ortho_TpG)
EpT_inpt = list(WT_EpT,ortho_EpT)
EpG_inpt = list(WT_EpG,ortho_EpG)
TpG_outpt=ratebytree(multiTree, TpG_inpt)
EpT_outpt=ratebytree(multiTree, EpT_inpt)
EpG_outpt=ratebytree(multiTree, EpG_inpt)

posthoc(TpG_outpt)
posthoc(EpT_outpt)
posthoc(EpG_outpt)

Trait_Tree_WT <- contMap(full_deut_tree,ortho_EpT,fsize=c(0.6,1),outline=FALSE)
Trait_Tree_WT <- contMap(full_deut_tree,WT_EpT,fsize=c(0.6,1),outline=FALSE)
contMap(full_deut_tree,ortho_EpT,fsize=c(0.6,1),outline=FALSE)

# TpG Tree
plotTree.wBars(Trait_Tree_WT$tree,WT_TpG,method="plotSimmap",
               tip.labels=TRUE,fsize=0.7,colors=Trait_Tree_WT$cols)
add.color.bar(200.0,Trait_Tree_WT$cols,title="trait value",lims=Trait_Tree_WT$lims,prompt=FALSE,
              x=0.5*par()$usr[1],y=0.5*par()$usr[3])

# EpT Tree
plotTree.wBars(Trait_Tree_WT$tree,WT_EpT,method="plotSimmap",
               tip.labels=TRUE,fsize=0.7,colors=Trait_Tree_WT$cols)
add.color.bar(200.0,Trait_Tree_WT$cols,title="trait value",lims=Trait_Tree_WT$lims,prompt=FALSE,
              x=0.5*par()$usr[1],y=0.5*par()$usr[3])

# EpG Tree
plotTree.wBars(Trait_Tree_WT$tree,WT_EpG,method="plotSimmap",
               tip.labels=TRUE,fsize=0.7,colors=Trait_Tree_WT$cols)
add.color.bar(200.0,Trait_Tree_WT$cols,title="trait value",lims=Trait_Tree_WT$lims,prompt=FALSE,
              x=0.5*par()$usr[1],y=0.5*par()$usr[3])

# Only reading it prepared files
#deuterostome_traits_check = read.csv("deuterostome_full_traits.csv", sep = ",")
#deuterostome_traits = read.table("deuterostome_full_traits.csv", sep = ",", row.names=1, header = TRUE)
#deuterostome_traits_check$species
