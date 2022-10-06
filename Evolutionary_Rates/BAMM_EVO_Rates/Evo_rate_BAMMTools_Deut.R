library(BAMMtools)
library(viridis)
library(coda)


setwd("/Volumes/GoogleDrive/My\ Drive/OrthoDB/Ortholog_Final_Files/Deuterostome_Complexity/Deuterostome_tree/Deuterostome_BAMM_Analysis")
#setwd("Your working directory here")

#Sanity check for mcmc convergence
#Convergence WT
deut_WT_mcmcout <- read.csv("deuterostome_WT_mean_EpT_mcmc_out.txt", header=T)
plot(deut_WT_mcmcout$logLik ~ deut_WT_mcmcout$generation)
deut_WT_burnstart <- floor(0.1 * nrow(deut_WT_mcmcout))
deut_WT_postburn <- deut_WT_mcmcout[deut_WT_burnstart:nrow(deut_WT_mcmcout), ]
plot(deut_WT_postburn$logLik ~ deut_WT_postburn$generation, ylim=c(-500,0))

# Check the ESS values and plot the WT prior
effectiveSize(deut_WT_postburn$N_shifts)
effectiveSize(deut_WT_postburn$logLik)
deut_WT_EpT_bfma <- computeBayesFactors(deut_WT_postburn, expectedNumberOfShifts=1, burnin=0.1)
plotPrior(deut_WT_postburn, expectedNumberOfShifts=1)

#Convergence ortholog
deut_ortho_mcmcout <- read.csv("deuterostome_ortho_mean_EpT_mcmc_out.txt", header=T)
plot(deut_ortho_mcmcout$logLik ~ deut_ortho_mcmcout$generation)
deut_ortho_burnstart <- floor(0.1 * nrow(deut_ortho_mcmcout))
deut_ortho_postburn <- deut_ortho_mcmcout[deut_ortho_burnstart:nrow(deut_ortho_mcmcout), ]
# Better visualization of convergence
plot(deut_ortho_postburn$logLik ~ deut_ortho_postburn$generation, ylim=c(-500,0))

# Check the ESS values and plot the orthologous prior
effectiveSize(deut_ortho_postburn$N_shifts)
effectiveSize(deut_ortho_postburn$logLik)
deut_ortho_EpT_bfma <- computeBayesFactors(deut_ortho_mcmcout, expectedNumberOfShifts=1, burnin=0.1)
plotPrior(deut_ortho_mcmcout, expectedNumberOfShifts=1)

#Analysis of rate shifts in the BAMM framework (WT)
#Marginal shift probabilities
deuterostome_tre <- read.tree("Deuterostome_Species_TimeTree_no_polytomy.tre")
deuterostome_ladder_tree <- ladderize(deuterostome_tre, right = TRUE)
edata_WT_deut <- getEventData(deuterostome_ladder_tree, eventdata = "deuterostome_WT_mean_EpT_event_data.txt", type = "trait", burnin = 0.1)
deut_marg_probs_WT_deut <- marginalShiftProbsTree(edata_WT_deut)
plot.phylo(deut_marg_probs_WT_deut)
edata_WT_deut_Best <- getBestShiftConfiguration(edata_WT_deut,expectedNumberOfShifts = 1, threshold = 100)
plot.bammdata(edata_WT_deut_Best, labels = TRUE, legend = TRUE, lwd = 3, pal = viridis(15))

# Plot
plot.bammdata(edata_WT_deut, labels = TRUE, legend = TRUE, lwd = 3, pal = viridis(15), color.interval = c(0,0.21))
summary(edata_WT_deut)

deuterostome_WT_css <- credibleShiftSet(edata_WT_deut, expectedNumberOfShifts = 1, threshold = 100)
plot.new()
par(mfrow=c(3,3))
plot.credibleshiftset(x = deuterostome_WT_css)
plotRateThroughTime(edata_WT_deut, ratetype = "auto")

#Analysis of rate shifts in the BAMM framework (Orthologs)
#Marginal shift probabilities
deuterostome_tre <- read.tree("Deuterostome_Species_TimeTree_no_polytomy.tre")
deuterostome_ladder_tree <- ladderize(deuterostome_tre, right = TRUE)
edata_ortho_deut <- getEventData(deuterostome_ladder_tree, eventdata = "deuterostome_ortho_mean_EpT_event_data.txt", type = "trait", burnin = 0.1)
deut_marg_probs_ortho_deut <- marginalShiftProbsTree(edata_ortho_deut)
plot.phylo(deut_marg_probs_ortho_deut)
edata_ortho_deut_Best <- getBestShiftConfiguration(edata_ortho_deut,expectedNumberOfShifts = 1, threshold = 100)
plot.bammdata(edata_ortho_deut_Best, labels = TRUE, legend = TRUE, lwd = 3, pal = viridis(15))

# Plot
plot.bammdata(edata_ortho_deut, labels = TRUE, legend = TRUE, lwd = 3, pal = viridis(15))
plot.bammdata(edata_ortho_deut, labels = TRUE, legend = TRUE, lwd = 3, pal = viridis(15), color.interval = c(0,0.21), direction = "leftwards")
summary(edata_ortho_deut)

deuterostome_ortho_css <- credibleShiftSet(edata_ortho_deut, expectedNumberOfShifts = 1, threshold = 100)
plot.new()
par(mfrow=c(3,3))
plot.credibleshiftset(x = deuterostome_ortho_css)
plotRateThroughTime(edata_ortho_deut, ratetype = "auto")
# ^^^ If input object is of type 'trait', ratetype can only be 'auto'.

write.tree(marg_probs_ortho_deut, file = "")
write.csv(as.data.frame(edata_ortho_deut), file = "edata_ortho_deut.csv")

# Pulling out the interesting Posterior Probabilities
# From the tree, get the labels
plot(deuterostome_ladder_tree)
nodelabels()
nodelabels(node = 121)
# Identify the nodes of interest pertaining to the rates shift
# simians are node =125


# Whole-transcriptome PP's
simian_WT_rates <- getCladeRates(ephy = edata_WT_deut, node = 125)
mean(simian_WT_rates$beta)
getTipRates(edata_WT_deut, statistic = "mean")# lists out all the species, looking at Branchiostoma_floridae8

# Ortho PP's
crown_primate_ortho_rates <- getCladeRates(ephy = edata_ortho_deut, node = 121)
mean(crown_primate_ortho_rates$beta)
getTipRates(edata_ortho_deut, statistic = "mean")


