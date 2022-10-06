library(BAMMtools)
library(viridis)
library(coda)

setwd("/Volumes/GoogleDrive/My Drive/OrthoDB/Ortholog_Final_Files/Fungi_Complexity/Fungi_tree/Fungi_BAMM_Analysis")
#setwd("Your working directory here")

#Sanity check for mcmc convergence
#Convergence WT
fungi_mcmcout <- read.csv("fungi_WT_mean_EpT_mcmc_out.txt", header=T)
plot(fungi_mcmcout$logLik ~ fungi_mcmcout$generation)
fungi_burnstart <- floor(0.1 * nrow(fungi_mcmcout))
fungi_postburn <- fungi_mcmcout[fungi_burnstart:nrow(fungi_mcmcout), ]
# Better visualization of convergence
plot(fungi_postburn$logLik ~ fungi_postburn$generation, ylim=c(-200,75))

effectiveSize(fungi_postburn$N_shifts)
effectiveSize(fungi_postburn$logLik)
fungi_WT_EpT_bfma <- computeBayesFactors(fungi_mcmcout, expectedNumberOfShifts=1, burnin=0.1)
plotPrior(fungi_mcmcout, expectedNumberOfShifts=1)

#Convergence ortholog
fungi_ortho_mcmcout <- read.csv("fungi_ortho_mean_EpT_mcmc_out.txt", header=T)
plot(fungi_ortho_mcmcout$logLik ~ fungi_ortho_mcmcout$generation)
fungi_ortho_burnstart <- floor(0.1 * nrow(fungi_ortho_mcmcout))
fungi_ortho_postburn <- fungi_ortho_mcmcout[fungi_ortho_burnstart:nrow(fungi_ortho_mcmcout), ]
# Better visualization of convergence
plot(fungi_ortho_postburn$logLik ~ fungi_ortho_postburn$generation, ylim=c(-200,75))

effectiveSize(fungi_ortho_postburn$N_shifts)
effectiveSize(fungi_ortho_postburn$logLik)
fungi_ortho_EpT_bfma <- computeBayesFactors(fungi_ortho_mcmcout, expectedNumberOfShifts=1, burnin=0.1)
plotPrior(fungi_ortho_mcmcout, expectedNumberOfShifts=1)

#Analysis of rate shifts in the BAMM framework (WT)
#Marginal shift probabilities
fungi_tre <- read.tree("fungi_class_pruned_no_polytomy.tre")
fungi_ladder_tree <- ladderize(fungi_tre, right = TRUE)
edata_WT_fungi <- getEventData(fungi_ladder_tree, eventdata = "fungi_WT_mean_EpT_event_data.txt", type = "trait", burnin = 0.1)
fungi_marg_probs_WT_fungi <- marginalShiftProbsTree(edata_WT_fungi)
plot.phylo(fungi_marg_probs_WT_fungi)
plot.bammdata(edata_WT_fungi, labels = TRUE, legend = TRUE, lwd = 3, pal = viridis(15))
#plot.bammdata(edata_WT_fungi, labels = FALSE, legend = FALSE, lwd = 3, color.interval = c(0,0.64), pal = viridis(15))
plot.bammdata(edata_WT_fungi, labels = TRUE, legend = TRUE, lwd = 3, color.interval = c(0,0.64), pal = viridis(15))
summary(edata_WT_fungi)

fungi_WT_css <- credibleShiftSet(edata_WT_fungi, expectedNumberOfShifts = 1, threshold = 100)
plot.credibleshiftset(x = fungi_WT_css)
plot.new()
plotRateThroughTime(edata_WT_fungi, ratetype = "auto")

write.tree(marg_probs_fungi, file = "")

#Analysis of rate shifts in the BAMM framework (Orthologs)
#Marginal shift probabilities
fungi_tre <- read.tree("fungi_class_pruned_no_polytomy.tre")
fungi_ladder_tree <- ladderize(fungi_tre, right = TRUE)
edata_ortho_fungi <- getEventData(fungi_ladder_tree, eventdata = "fungi_ortho_mean_EpT_event_data.txt", type = "trait", burnin = 0.1)
fungi_marg_probs_ortho_fungi <- marginalShiftProbsTree(edata_ortho_fungi)
plot.phylo(fungi_marg_probs_ortho_fungi)
plot.bammdata(edata_ortho_fungi, labels = TRUE, legend = TRUE, lwd = 3, pal = viridis(15))
#plot.bammdata(edata_ortho_fungi, labels = TRUE, legend = TRUE, lwd = 3, direction = "leftwards", color.interval = c(0,0.64), pal = viridis(15))
summary(edata_ortho_fungi)

fungi_ortho_css <- credibleShiftSet(edata_ortho_fungi, expectedNumberOfShifts = 1, threshold = 100)
plot.credibleshiftset(x = fungi_ortho_css)
plot.new()
plotRateThroughTime(edata_ortho_fungi, ratetype = "auto")

write.tree(marg_probs_ortho_fungi, file = "")

# Pulling out the interesting Posterior Probabilities
# From the tree, get the labels
plot(fungi_ladder_tree)
nodelabels()
# Identify the nodes of interest pertaining to the rates shift
# basal node are node = 17
# Pneumocystidomycetes WT and ortho dataset

# Whole-transcriptome PP's
basal_fungi_WT_rates <- getCladeRates(ephy = edata_WT_fungi, node = 17)
median(basal_fungi_WT_rates$beta)
#0.03008925
getTipRates(edata_WT_fungi, statistic = "mean")
# Pneumocystidomycetes - 0.3957112664

# Ortho PP's
basal_fungi_ortho_rates <- getCladeRates(ephy = edata_ortho_fungi, node = 17)
mean(basal_fungi_ortho_rates$beta)
# 0.04266893
getTipRates(edata_ortho_fungi, statistic = "mean")
# Pneumocystidomycetes - 0.5699981151


