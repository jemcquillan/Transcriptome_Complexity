library(BAMMtools)
library(viridis)
library(coda)

setwd("/Volumes/GoogleDrive/My Drive/OrthoDB/Ortholog_Final_Files/Plantae_Complexity/Plantae_tree/Plantae_BAMM_Analysis")
#setwd("Your working directory here")

#Sanity check for mcmc convergence
#Convergence WT
plant_mcmcout <- read.csv("plantae_WT_mean_EpT_mcmc_out.txt", header=T)
plot(plant_mcmcout$logLik ~ plant_mcmcout$generation)
plant_burnstart <- floor(0.1 * nrow(plant_mcmcout))
plant_postburn <- plant_mcmcout[plant_burnstart:nrow(plant_mcmcout), ]
# Better visualization of convergence
plot(plant_postburn$logLik ~ plant_postburn$generation, ylim=c(-200,100))

# Check the ESS values and plot the WT prior
effectiveSize(plant_postburn$N_shifts)
effectiveSize(plant_postburn$logLik)
plant_WT_EpT_bfma <- computeBayesFactors(plant_mcmcout, expectedNumberOfShifts=1, burnin=0.1)
plotPrior(plant_mcmcout, expectedNumberOfShifts=1)

#Convergence ortholog
plant_ortho_mcmcout <- read.csv("plantae_ortho_mean_EpT_mcmc_out.txt", header=T)
plot(plant_ortho_mcmcout$logLik ~ plant_ortho_mcmcout$generation)
plant_ortho_burnstart <- floor(0.1 * nrow(plant_ortho_mcmcout))
plant_ortho_postburn <- plant_ortho_mcmcout[plant_burnstart:nrow(plant_ortho_mcmcout), ]
# Better visualization of convergence
plot(plant_ortho_postburn$logLik ~ plant_ortho_postburn$generation, ylim=c(-100,50))

# Check the ESS values and plot the orthologous prior
effectiveSize(plant_ortho_postburn$N_shifts)
effectiveSize(plant_ortho_postburn$logLik)
plant_ortho_EpT_bfma <- computeBayesFactors(plant_ortho_mcmcout, expectedNumberOfShifts=1, burnin=0.1)
plotPrior(plant_ortho_mcmcout, expectedNumberOfShifts=1)

#Analysis of rate shifts in the BAMM framework (WT)
#Marginal shift probabilities
plantae_tre <- read.tree("Plantae_Species_Tree_list_no_polytomy.tre")
plantae_ladder_tree <- ladderize(plantae_tre, right = TRUE)
edata_WT_plant <- getEventData(plantae_ladder_tree, eventdata = "plantae_WT_mean_EpT_event_data.txt", type = "trait", burnin = 0.1)
plant_marg_probs_WT_plant <- marginalShiftProbsTree(edata_WT_plant)
plot.phylo(plant_marg_probs_WT_plant)
plot.bammdata(edata_WT_plant, labels = TRUE, legend = TRUE, lwd = 3, pal = viridis(15))
#plot.bammdata(edata_WT_plant, labels = FALSE, legend = FALSE, lwd = 3, color.interval = c(0,0.21), pal = viridis(15))
plot.bammdata(edata_WT_plant, labels = TRUE, legend = TRUE, lwd = 3, color.interval = c(0,0.21), pal = viridis(15))
summary(edata_WT_plant)

plant_WT_css <- credibleShiftSet(edata_WT_plant, expectedNumberOfShifts = 1, threshold = 100)
plot.credibleshiftset(x = plant_WT_css)
plot.new()
plotRateThroughTime(edata_WT_plant, ratetype = "auto")

write.tree(marg_probs_plant, file = "")

#Analysis of rate shifts in the BAMM framework (Orthologs)
#Marginal shift probabilities
plantae_tre <- read.tree("Plantae_Species_Tree_list_no_polytomy.tre")
plantae_ladder_tree <- ladderize(plantae_tre, right = TRUE)
edata_ortho_plant <- getEventData(plantae_ladder_tree, eventdata = "plantae_ortho_mean_EpT_event_data.txt", type = "trait", burnin = 0.1)
plant_marg_probs_ortho_plant <- marginalShiftProbsTree(edata_ortho_plant)
plot.phylo(plant_marg_probs_ortho_plant)
#plot.bammdata(edata_ortho_plant, labels = TRUE, legend = TRUE, lwd = 3, direction = "leftwards", pal = viridis(15))
#plot.bammdata(edata_ortho_plant, labels = TRUE, legend = TRUE, lwd = 3, direction = "leftwards", color.interval = c(0,0.21), pal = viridis(15))
plot.bammdata(edata_ortho_plant, labels = TRUE, legend = TRUE, lwd = 3, pal = viridis(15))
summary(edata_ortho_plant)

plant_ortho_css <- credibleShiftSet(edata_ortho_plant, expectedNumberOfShifts = 1, threshold = 100)
plot.credibleshiftset(x = plant_ortho_css)
plot.new()
plotRateThroughTime(edata_ortho_plant, ratetype = "auto")

write.tree(marg_probs_ortho_plant, file = "")

# Pulling out the interesting Posterior Probabilities
# From the tree, get the labels
plot(plantae_ladder_tree)
nodelabels()
# Identify the nodes of interest pertaining to the rates shift
# basal node are node = 45
# Arabidopsis thaliana and Camelina sativa is at node = 64
# Quercus suber (cork oak) in WT dataset
# Zea mays high in orthologous dataset

# Whole-transcriptome PP's
Arabidopsis_C.sativa_WT_rates <- getCladeRates(ephy = edata_WT_plant, node = 64)
mean(Arabidopsis_C.sativa_WT_rates$beta)

Basal_viriplantae_WT_rates <- getCladeRates(ephy = edata_WT_plant, node = 45)
mean(Basal_viriplantae_WT_rates$beta)

getTipRates(edata_WT_plant, statistic = "mean")
# Quercus_suber - 0.109740157

# Ortho PP's
getTipRates(edata_ortho_plant, statistic = "mean")
# Zea mayes - 0.277881691

Basal_viriplantae_ortho_rates <- getCladeRates(ephy = edata_ortho_plant, node = 45)
mean(Basal_viriplantae_ortho_rates$beta)
