library(BAMMtools)
library(viridis)
library(coda)

setwd("/Volumes/GoogleDrive/My Drive/OrthoDB/Ortholog_Final_Files/Drosophila_Complexity/Drosophila_tree/Drosophila_BAMM_Analysis")
#setwd("Your working directory here")

#Sanity check for mcmc convergence
#Convergence WT
dros_mcmcout <- read.csv("drosophila_WT_mean_EpT_mcmc_out.txt", header=T)
plot(dros_mcmcout$logLik ~ dros_mcmcout$generation)
dros_burnstart <- floor(0.1 * nrow(dros_mcmcout))
dros_postburn <- dros_mcmcout[dros_burnstart:nrow(dros_mcmcout), ]
# Better visualization of convergence
plot(dros_postburn$logLik ~ dros_postburn$generation, ylim=c(-200,10))

# Check the ESS values and plot the WT prior
effectiveSize(dros_postburn$N_shifts)
effectiveSize(dros_postburn$logLik)
dros_WT_EpT_bfma <- computeBayesFactors(dros_mcmcout, expectedNumberOfShifts=1, burnin=0.1)
plotPrior(dros_mcmcout, expectedNumberOfShifts=1)

#Convergence ortholog
dros_ortho_mcmcout <- read.csv("drosophila_ortho_mean_EpT_mcmc_out.txt", header=T)
plot(dros_ortho_mcmcout$logLik ~ dros_ortho_mcmcout$generation)
dros_ortho_burnstart <- floor(0.1 * nrow(dros_ortho_mcmcout))
dros_ortho_postburn <- dros_ortho_mcmcout[dros_burnstart:nrow(dros_ortho_mcmcout), ]
# Better visualization of convergence
plot(dros_ortho_postburn$logLik ~ dros_ortho_postburn$generation, ylim=c(-200,10))

# Check the ESS values and plot the orthologous prior
effectiveSize(dros_ortho_postburn$N_shifts)
effectiveSize(dros_ortho_postburn$logLik)
dros_ortho_EpT_bfma <- computeBayesFactors(dros_ortho_mcmcout, expectedNumberOfShifts=1, burnin=0.1)
plotPrior(dros_ortho_mcmcout, expectedNumberOfShifts=1)

#Analysis of rate shifts in the BAMM framework (WT)
#Marginal shift probabilities
drosophila_tre <- read.tree("drosophila_Species_TimeTree_no_polytomy.tre")
drosophila_ladder_tree <- ladderize(drosophila_tre, right = TRUE)
edata_WT_dros <- getEventData(drosophila_ladder_tree, eventdata = "drosophila_WT_mean_EpT_event_data.txt", type = "trait", burnin = 0.1)
dros_marg_probs_WT_dros <- marginalShiftProbsTree(edata_WT_dros)
plot.phylo(dros_marg_probs_WT_dros)
plot.bammdata(edata_WT_dros, labels = TRUE, lwd = 3, pal = viridis(15), legend = TRUE)
plot.bammdata(edata_WT_dros, labels = FALSE, lwd = 3, pal = viridis(15), color.interval = c(0,0.92))
#Temp plot
plot.bammdata(edata_WT_dros, labels = TRUE, legend = TRUE, lwd = 3, color.interval = c(0,0.92), pal = viridis(15))
summary(edata_WT_dros)

drosophila_WT_css <- credibleShiftSet(edata_WT_dros, expectedNumberOfShifts = 1, threshold = 100)
par(mfrow=c(3,3))
plot.credibleshiftset(x = drosophila_WT_css)
plot.new()
plotRateThroughTime(edata_WT_dros, ratetype = "auto")

write.tree(dros_marg_probs_WT_dros, file = "")

#Analysis of rate shifts in the BAMM framework (Orthologs)
#Marginal shift probabilities
drosophila_tre <- read.tree("drosophila_Species_TimeTree_no_polytomy.tre")
drosophila_ladder_tree <- ladderize(drosophila_tre, right = TRUE)
edata_ortho_dros <- getEventData(drosophila_ladder_tree, eventdata = "drosophila_ortho_mean_EpT_event_data.txt", type = "trait", burnin = 0.1)
dros_marg_probs_ortho_dros <- marginalShiftProbsTree(edata_ortho_dros)
plot.phylo(dros_marg_probs_ortho_dros)
plot.bammdata(edata_ortho_dros, labels = TRUE, legend = TRUE, lwd = 3, direction = "leftwards", pal = viridis(15))
plot.bammdata(edata_ortho_dros, labels = TRUE, legend = TRUE, lwd = 3, color.interval = c(0,0.92), direction = "leftwards", pal = viridis(15))
summary(edata_ortho_dros)

drosophila_ortho_css <- credibleShiftSet(edata_ortho_dros, expectedNumberOfShifts = 1, threshold = 100)
par(mfrow=c(3,3))
plot.credibleshiftset(x = drosophila_WT_css)
plot.new()
plotRateThroughTime(edata_WT_dros, ratetype = "auto")

write.tree(dros_marg_probs_ortho_dros, file = "")

# Pulling out the interesting Posterior Probabilities
# From the tree, get the labels
plot(drosophila_ladder_tree)
nodelabels()
# Identify the nodes of interest pertaining to the rates shift
# D. sechellia and D. simulans are node = 20


# Whole-transcriptome PP's
sim_schell_WT_rates <- getCladeRates(ephy = edata_WT_dros, node = 20)
mean(sim_schell_WT_rates$beta)

# Ortho PP's
sim_schell_ortho_rates <- getCladeRates(ephy = edata_ortho_dros, node = 20)
mean(sim_schell_ortho_rates$beta)

