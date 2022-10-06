library(ape)

setwd("/Volumes/GoogleDrive/My\ Drive/OrthoDB/Ortholog_Final_Files/Fungi_Complexity/Fungi_tree")

#Read in the tree file from TimeTree (.tre/.nwk)
fungi_order_tree <- read.tree(file = "fungi_order.tre")
#plot(fungi_order_tree)

fungi_class_tree <- read.tree(file = "fungi_class.tre")
plot(fungi_class_tree)

#Read the contents of the file
cat(readLines("fungi_class.tre"), sep = ";")

updated_fungi_tree <- drop.tip(fungi_class_tree, tip = c("Agaricostilbomycetes","Arthoniomycetes","Atractiellomycetes","Basidiobolomycetes","Blastocladiomycetes","Coniocybomycetes","Dacrymycetes","Entomophthoromycetes","Exobasidiomycetes","Geoglossomycetes","Glomeromycetes","Lecanoromycetes","Lichinomycetes","Orbiliomycetes","Pezizomycetes","Taphrinomycetes"))
plot(updated_fungi_tree)
