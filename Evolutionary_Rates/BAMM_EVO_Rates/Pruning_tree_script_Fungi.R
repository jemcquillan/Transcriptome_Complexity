library(phytools)

setwd("/Volumes/GoogleDrive/My\ Drive/OrthoDB/Ortholog_Final_Files/Fungi_Complexity/Fungi_tree")

fungi_class_tree <- ape::read.tree("fungi_class.tre")
plot(fungi_class_tree)

classes_keep <- c('Saccharomycetes','Chytridiomycetes','Eurotiomycetes','Dothideomycetes','Sordariomycetes','Agaricomycetes','Tremellomycetes','Leotiomycetes','Ustilaginomycetes','Malasseziomycetes','Pucciniomycetes','Neocallimastigomycetes','Pneumocystidomycetes','Microbotryomycetes','Schizosaccharomycetes','Wallemiomycetes')

fungi_pruned.tree <- drop.tip(fungi_class_tree,fungi_class_tree$tip.label[-match(classes_keep, fungi_class_tree$tip.label)])
plot(fungi_pruned.tree)
write.tree(phy = fungi_pruned.tree, file = "fungi_class_pruned.tre")
