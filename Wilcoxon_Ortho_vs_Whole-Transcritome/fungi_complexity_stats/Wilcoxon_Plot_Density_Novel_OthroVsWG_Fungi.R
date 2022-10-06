library(tidyverse)
library(rstatix)
library(ggpubr)
library(stringr)
library(poolr)

setwd("/Volumes/GoogleDrive/My Drive/OrthoDB/Ortholog_Final_Files/Fungi_Complexity/")
#setwd("<your path here>")

# Creates a list of files from the directory that has the appropriate metrics of novel genes
novel_EpT <- list.files('fungi_complexity_stats/fungi_novel_genes/fungi_novel_genes_EpT/')
novel_TpG_EpG <- list.files('fungi_complexity_stats/fungi_novel_genes/fungi_novel_genes_TpG_EpG/')

# Creates a list of files from the directory that has the appropriate metrics of orthologs
ortho_EpT <- list.files('fungi_complexity_stats/fungi_ortholog/fungi_ortho_EpT/')
ortho_TpG_EpG <- list.files('fungi_complexity_stats/fungi_ortholog/fungi_ortho_TpG_EpG/')

# Creates a list of files from the directory that has the appropriate metrics of whole transcriptomes
wt_EpT <- list.files('fungi_complexity_stats/fungi_whole_genome/fungi_whole_genome_EpT/')
wt_TpG_EpG <- list.files('fungi_complexity_stats/fungi_whole_genome/fungi_whole_genome_TpG_EpG/')

# Empty lists that will have all the plots to be printed into a pdf at the end.
plot_EpT <- list()
plot_EpG <- list()
plot_TpG <- list()
plot_EpT_trunc <- list()
plot_EpG_trunc <- list()
plot_TpG_trunc <- list()

plot_novel_EpT <- list()
plot_novel_EpG <- list()
plot_novel_TpG <- list()
plot_novel_EpT_trunc <- list()
plot_novel_EpG_trunc <- list()
plot_novel_TpG_trunc <- list()

#Create your output directory
output_directory <- "<your path here>/fungi_complexity_stats/fungi_density_plots"

pvals_list <- list()

# For loop that iterates over all the files in ortho_EpT.
for(i in 1:length(ortho_EpT)){
  
  # This creates a dataframe with all the data from each EpT metric file
  ortho_EpT_df <- read.csv(file.path('fungi_complexity_stats/fungi_ortholog/fungi_ortho_EpT/',ortho_EpT[i]), sep = ',')
  # This creates and appends a column names group, and labels the orthologs 'ortho'
  ortho_EpT_df$group <- 'ortho'
  # Removes the extra naming from the file giving just the GCF number
  # e.g., "GCF_000001405.39_GRCh38.p13_genomic.transcriptome_counts_transcript_level.csv"
  # becomes "GCF_000001405.39_GRCh38.p13"
  file_name <- stringr::str_remove(ortho_EpT[i],'_genomic.*')
  
  # Same this here as above. Creates a data frame with all the data, and adds column group
  # labeled ortho for all the orthologs for the other 2 metrics
  ortho_TpG_EpG_df <- read.csv(file.path('fungi_complexity_stats/fungi_ortholog/fungi_ortho_TpG_EpG/',ortho_TpG_EpG[i]), sep = ',')
  ortho_TpG_EpG_df$group <- 'ortho'
  
  # These lines are going to do the same as above except with the whole transcriptome. It will also
  # add in the group column but label it 'whole transcriptome'
  wt_EpT_df <- read.csv(file.path('fungi_complexity_stats/fungi_whole_genome/fungi_whole_genome_EpT/',wt_EpT[i]), sep = ',')
  wt_EpT_df$group <- 'whole transcriptome'
  
  wt_TpG_EpG_df <- read.csv(file.path('fungi_complexity_stats/fungi_whole_genome/fungi_whole_genome_TpG_EpG/',wt_TpG_EpG[i]), sep = ',')
  wt_TpG_EpG_df$group <- 'whole transcriptome'
  
  # This creates a dataframe with all the data from each EpT metric file for novel genes
  novel_EpT_df <- read.csv(file.path('fungi_complexity_stats/fungi_novel_genes/fungi_novel_genes_EpT/',novel_EpT[i]), sep = ',')
  # This creates and appends a column names group, and labels the orthologs 'ortho'
  novel_EpT_df$group <- 'novel'
  file_name <- stringr::str_remove(novel_EpT[i],'_genomic.*')
  
  # Same this here as above. Creates a data frame with all the data, and adds column group
  # labeled ortho for all the orthologs for the other 2 metrics
  novel_TpG_EpG_df <- read.csv(file.path('fungi_complexity_stats/fungi_novel_genes/fungi_novel_genes_TpG_EpG/',novel_TpG_EpG[i]), sep = ',')
  novel_TpG_EpG_df$group <- 'novel'
  
  #Make Novel genes density plot
  #EpT plots for novel genes
  novel_ept_plot_trunc <- ggplot(novel_EpT_df, aes_string('num_exon', fill = 'group', color = 'group')) +
    geom_density(bw = "bcv",color=3,fill=3,alpha=0.4) + xlim(0,30) + 
    labs(title = paste0(file_name,'\n','EpT'), subtitle = paste0('Novel Genes'))
  plot_novel_EpT_trunc[[i]] <- novel_ept_plot_trunc
  
  novel_ept_plot <- ggplot(novel_EpT_df, aes_string('num_exon', fill = 'group', color = 'group')) +
    geom_density(bw = "bcv",color=3,fill=3,alpha=0.4) + 
    labs(title = paste0(file_name,'\n','EpT'), subtitle = paste0('Novel Genes'))
  plot_novel_EpT[[i]] <- novel_ept_plot
  
  #Plot EpG and TpG novel genes plot
  novel_tpg_plot_trunc <- ggplot(novel_TpG_EpG_df, aes_string('num_transcript', fill = 'group', color = 'group'), alpha(.25)) +
    geom_density(bw = "bcv",color=3,fill=3,alpha=0.4) + xlim(0,30) +
    labs(title = paste0(file_name,'\n','TpG'))
  plot_novel_TpG_trunc[[i]] <- novel_tpg_plot_trunc
  
  #TpG ploits for novel genes
  novel_tpg_plot <- ggplot(novel_TpG_EpG_df, aes_string('num_transcript', fill = 'group', color = 'group'), alpha(.25)) +
    geom_density(bw = "bcv",color=3,fill=3,alpha=0.4) + 
    labs(title = paste0(file_name,'\n','EpT'), subtitle = paste0('Novel Genes'))
  plot_novel_TpG[[i]] <- novel_tpg_plot
  
  #EpG plots for novel genes
  novel_epg_plot_trunc <- ggplot(novel_TpG_EpG_df, aes_string('num_uniq_exon', fill = 'group', color = 'group')) +
    geom_density(bw = "bcv",color=3,fill=3,alpha=0.4) + xlim(0,30) +
    labs(title = paste0(file_name,'\n','EpT'), subtitle = paste0('Novel Genes'))
  plot_novel_EpG_trunc[[i]] <- novel_epg_plot_trunc
  
  novel_epg_plot <- ggplot(novel_TpG_EpG_df, aes_string('num_uniq_exon', fill = 'group', color = 'group')) +
    geom_density(bw = "bcv",color=3,fill=3,alpha=0.4) + 
    labs(title = paste0(file_name,'\n','EpT'), subtitle = paste0('Novel Genes'))
  plot_novel_EpG[[i]] <- novel_epg_plot
  
  #Combine  the whole transcriptome and orthologs into on df with the metric EpT
  EpT_df <- rbind(wt_EpT_df,ortho_EpT_df)
  EpT_df$group <- factor(EpT_df$group, levels = c("ortho", "whole transcriptome"))
  
  #K-S Between EpT of ortho vs wt
  #pval_ks <-ks.test(EpT_df$num_exon[EpT_df$group=='ortho'], EpT_df$num_exon[EpT_df$group=='whole transcriptome'])
  # Performs a wilcoxon test between values 'ortho' vs 'whole transcriptome'
  pval1 <- wilcox.test(EpT_df$num_exon[EpT_df$group=='ortho'], EpT_df$num_exon[EpT_df$group=='whole transcriptome'], p.adjust.methods = "bonferroni")
  
  ####
  #plot(density(thing1[,1]))# Make a swooshy density plot and make it a proportion (normalize the density).
  ####
  
  ept_plot_trunc <- ggplot(EpT_df, aes_string('num_exon', fill = 'group', color = 'group')) +
    geom_density(alpha=0.4) + xlim(0,30) + 
    labs(title = paste0(file_name,'\n','EpT'),
         subtitle = paste0('Wilcoxon p-value = ',
                           format(pval1$p.value,digit = 5),", W = ",
                           format(pval1$statistic,digit = 4)))
  plot_EpT_trunc[[i]] <- ept_plot_trunc
  
  ept_plot <- ggplot(EpT_df, aes_string('num_exon', fill = 'group', color = 'group')) +
    geom_density(alpha=0.4) + 
    labs(title = paste0(file_name,'\n','EpT'),
         subtitle = paste0('Wilcoxon p-value = ',
                           format(pval1$p.value,digit = 5),", W = ",
                           format(pval1$statistic,digit = 4)))
  plot_EpT[[i]] <- ept_plot
  
  #comine TpG_EpG_df
  TpG_EpG_df <- rbind(ortho_TpG_EpG_df, wt_TpG_EpG_df)
  
  #K-S Between EpG of ortho vs wt
  #pval_ks <-ks.test(TpG_EpG_df$num_transcript[TpG_EpG_df$group=='ortho'], TpG_EpG_df$num_transcript[TpG_EpG_df$group=='whole transcriptome'])
  pval2 <-wilcox.test(TpG_EpG_df$num_transcript[TpG_EpG_df$group=='ortho'], TpG_EpG_df$num_transcript[TpG_EpG_df$group=='whole transcriptome'], p.adjust.methods = "bonferroni")
  
  tpg_plot_trunc <- ggplot(TpG_EpG_df, aes_string('num_transcript', fill = 'group', color = 'group'), alpha(.25)) +
    geom_density(alpha=0.4) + xlim(0,30) +
    labs(title = paste0(file_name,'\n','TpG'),subtitle = paste0('Wilcoxon p-value = ',
                                                                format(pval2$p.value,digit = 5),", W = ",
                                                                format(pval2$statistic,digit = 4)))
  plot_TpG_trunc[[i]] <- tpg_plot_trunc
  
  tpg_plot <- ggplot(TpG_EpG_df, aes_string('num_transcript', fill = 'group', color = 'group'), alpha(.25)) +
    geom_density(alpha=0.4) + 
    labs(title = paste0(file_name,'\n','TpG'),subtitle = paste0('Wilcoxon p-value = ',
                                                                format(pval2$p.value,digit = 5),", W = ",
                                                                format(pval2$statistic,digit = 4)))
  plot_TpG[[i]] <- tpg_plot
  
  #K-S Between TpG of ortho vs wt
  #pval_ks <-ks.test(TpG_EpG_df$num_uniq_exon[TpG_EpG_df$group=='ortho'], TpG_EpG_df$num_uniq_exon[TpG_EpG_df$group=='whole transcriptome'])
  pval3 <-wilcox.test(TpG_EpG_df$num_uniq_exon[TpG_EpG_df$group=='ortho'], TpG_EpG_df$num_uniq_exon[TpG_EpG_df$group=='whole transcriptome'], p.adjust.methods = "bonferroni")
  
  epg_plot_trunc <- ggplot(TpG_EpG_df, aes_string('num_uniq_exon', fill = 'group', color = 'group')) +
    geom_density(alpha=0.4) + xlim(0,30) +
    labs(title = paste0(file_name,'\n','EpG'),subtitle = paste0('Wilcoxon p-value = ',
                                                                format(pval3$p.value,digit = 5),", W = ",
                                                                format(pval3$statistic,digit = 4)))
  plot_EpG_trunc[[i]] <- epg_plot_trunc
  
  epg_plot <- ggplot(TpG_EpG_df, aes_string('num_uniq_exon', fill = 'group', color = 'group')) +
    geom_density(alpha=0.4) + 
    labs(title = paste0(file_name,'\n','EpG'),subtitle = paste0('Wilcoxon p-value = ',
                                                                format(pval3$p.value,digit = 5),", W = ",
                                                                format(pval3$statistic,digit = 4)))
  plot_EpG[[i]] <- epg_plot
  
  pvals_list[[i]] <- data.frame(p_value_Ept = pval1$p.value, p_value_TPG = pval2$p.value, p_value_EpG = pval3$p.value)
  
}

pvals_df <- dplyr::bind_rows(pvals_list)

# Since TpG has NA's when comaring the orhtos vs the WT, they equal, so the stat is NA.
# Given TpG is NAs, we need to drop the Na's from TpG to do the Fisher Method
TpG_drop_NaNs_pavls_df <- pvals_df %>% drop_na(p_value_TPG)

Ept_Bon = p.adjust(pvals_df$p_value_Ept, method = "bonferroni")
EpG_Bon = p.adjust(pvals_df$p_value_EpG, method = "bonferroni")
TpG_Bon = p.adjust(pvals_df$p_value_TPG, method = "bonferroni")

TpG_Bon_Fisher_meth = fisher(p.adjust(TpG_drop_NaNs_pavls_df$p_value_TPG, method = "bonferroni"))
EpT_Bon_Fisher_meth = fisher(p.adjust(pvals_df$p_value_Ept, method = "bonferroni"))
EpG_Bon_Fisher_meth = fisher(p.adjust(pvals_df$p_value_EpG, method = "bonferroni"))

#Print out the plots into PDFs
pdf(paste0(output_directory,'/Fungi_EpT_plot_novel_trunc.pdf'),7,5)
for(i in 1:length(ortho_EpT)){
  print(plot_novel_EpT_trunc[[i]])
}
dev.off()

pdf(paste0(output_directory,'/Fungi_EpT_plot_novel.pdf'),7,5)
for(i in 1:length(ortho_EpT)){
  print(plot_novel_EpT[[i]])
}
dev.off()

pdf(paste0(output_directory,'/Fungi_TpG_plot_novel_trunc.pdf'),7,5)
for(i in 1:length(ortho_EpT)){
  print(plot_novel_TpG_trunc[[i]])
}
dev.off()

pdf(paste0(output_directory,'/Fungi_TpG_plot_novel.pdf'),7,5)
for(i in 1:length(ortho_EpT)){
  print(plot_novel_TpG[[i]])
}
dev.off()

pdf(paste0(output_directory,'/Fungi_EpG_plot_novel_trunc.pdf'),7,5)
for(i in 1:length(ortho_EpT)){
  print(plot_novel_EpG_trunc[[i]])
}
dev.off()

pdf(paste0(output_directory,'/Fungi_EpG_plot_novel.pdf'),7,5)
for(i in 1:length(ortho_EpT)){
  print(plot_novel_EpG[[i]])
}
dev.off()



pdf(paste0(output_directory,'/Fungi_EpT_plot_wilcoxon.pdf'),7,5)
for(i in 1:length(ortho_EpT)){
  print(plot_EpT[[i]])
}
dev.off()

pdf(paste0(output_directory,'/Fungi_EpG_plot_wilcoxon.pdf'),7,5)
for(i in 1:length(ortho_EpT)){
  print(plot_EpG[[i]])
}
dev.off()

pdf(paste0(output_directory,'/Fungi_TpG_plot_wilcoxon.pdf'),7,5)
for(i in 1:length(ortho_EpT)){
  print(plot_TpG[[i]])
}
dev.off()

pdf(paste0(output_directory,'/Fungi_EpT_truncated_plot_wilcoxon.pdf'),7,5)
for(i in 1:length(ortho_EpT)){
  print(plot_EpT_trunc[[i]])
}
dev.off()

pdf(paste0(output_directory,'/Fungi_EpG_truncated_plot_wilcoxon.pdf'),7,5)
for(i in 1:length(ortho_EpT)){
  print(plot_EpG_trunc[[i]])
}
dev.off()

pdf(paste0(output_directory,'/Fungi_TpG_truncated_plot_wilcoxon.pdf'),7,5)
for(i in 1:length(ortho_EpT)){
  print(plot_TpG_trunc[[i]])
}
dev.off()
