---
title: "Esimating Transcriptome Complexities Supplementary Data"
author: "James Titus-McQuillan"
date: '2022-11-14'
output: html_document
---
# Esimating Transcriptome Complexities Supplementary Data:
# Transcriptome_Complexity
Scripts and code used for transcriptome complexity
---
The contents presented here are graphs from analyses suitable for supplementary material and the raw data used for such analyses.
## Description of the Data and file structure
Data and results presented here are dervived from code and analses found on these two GitHub repositories and from TranD, linked below.
### Mean_Complexity_Metric_Boxplots.pdf
This is a PDF illustrating the results between and among group complexity dynamics and the relative complexty performance.
### Whole-Transcriptome_VS_Ortholog_Density_Plots.pdf
This PDF contains figures of density plots between whole-genome (red) and orthologs (light blue) for
complexity metrics (TpG, EpT, and EpG). Here we show the difference in complexity between whole-transcriptome data and orthologous data and if they are signifcantly different from eachother. The densities on the y-axis and genetic element counts on the x-axis. Plots include every organism used in this study.
### Novel_Gene_Density_Plots.pdf
This PDF contains figures of density plots density plots for all novel gene distributions (light green) for complexity metrics (TpG, EpT, and EpG). The densities on the y-axis and genetic element counts on the xaxis. Plots include every organism used in this study.
### Raw data output from TranD describe comeplexities
Each group (Deuterostoma, *Drosophila*, Planae, and Fungi) contains three csv files for whole-transcriptome complexities, ortholog complexities, and novel genetic element complexities. Each file has meta data for each organisms parsed from RefSeq, GenBank, or FlyBase [species, commonName, group, database, taxid, link, filePrefix]. The database section defines where the data (GTF) was taken and link is the FPT link to download the dataset used in this study. All other metrics are derived from TranD describe complexity function. Note the links in othrologous and novel complexity csv's are directed to the parent GTF files (Whole-transcriptome data found from the corresponding database). All parsed GTF data (orthogous and novel genetic elements) can be found on GitHub Transcriptome Complexity repsository, listed below. 
## Sharing/access Information
### Transcriptome Complexity
https://github.com/jemcquillan/Transcriptome_Complexity
### OrthoDB Parser
https://github.com/jemcquillan/OrthoDB_Parser
### TranD
https://github.com/McIntyre-Lab/TranD/wiki
