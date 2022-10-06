library("ggpubr")
library('ggplot2')

setwd("/Volumes/GoogleDrive/My Drive/OrthoDB/Ortholog_Final_Files/")
#setwd("set your working directory here")

# Input Data Converted to Dataframes

Fungi_data <- read.csv("Fungi_Complexity/Fungi_merged_complexity/Fungi_WT_gtf_species_full_list.csv")
# Fungi_data is not relavant here because of annotation (i.e. TpG is all 1, and is a poor measure currently with the data available.)

# Run ANOVA
Fungi_species <- Fungi_data$species
Fungi_meanTpG <- Fungi_data$mean_transcriptPerGene

Fungi_one_way <- aov(mean_transcriptPerGene ~ Fungi_species, data = Fungi_data)
summary(Fungi_one_way)

# Run Shapiro Wilks test of normality
shapiro.test(Fungi_data$mean_exonPerTranscript)

#Fungi plots

Fungi_fit_TpG_EpT <- summary(lm(mean_exonPerTranscript~mean_transcriptPerGene,data = Fungi_data))
Fungi_plot_TpG_EpT <- ggplot(Fungi_data, aes_string(x = "mean_transcriptPerGene", y = "mean_exonPerTranscript"))+
  geom_point(aes_string(color = "group"))+
  geom_smooth(method = "lm",color = "black")+
  theme_classic()+
  labs(title = 'Fungi TpG vs EpT',x = "Mean Transcript per Gene", y = "Mean Exon per Transcript",
       subtitle = paste0('R2 = ',format(Fungi_fit_TpG_EpT$r.squared,digits = 3),' p = ',format(Fungi_fit_TpG_EpT$coefficients[2,4],digits = 3)))

Fungi_fit_TpG_EpG <- summary(lm(mean_exonPerGene~mean_transcriptPerGene,data = Fungi_data))
Fungi_plot_TpG_EpG <- ggplot(Fungi_data, aes_string(x = "mean_transcriptPerGene", y = "mean_exonPerGene"))+
  geom_point(aes_string(color = "group"))+
  geom_smooth(method = "lm",color = "black")+
  theme_classic()+
  labs(title = 'Fungi TpG vs EpG',x= "Mean Transcript per Gene", y = "Mean Exon per Gene",
       subtitle = paste0('R2 = ',format(Fungi_fit_TpG_EpG$r.squared,digits = 3),' p = ',format(Fungi_fit_TpG_EpG$coefficients[2,4],digits = 3)))

Fungi_fit_EpT_EpG <- summary(lm(mean_exonPerGene~mean_exonPerTranscript,data = Fungi_data))
Fungi_plot_EpT_EpG <- ggplot(Fungi_data, aes_string(x = "mean_exonPerTranscript", y = "mean_exonPerGene"))+
  geom_point(aes_string(color = "group"))+
  geom_smooth(method = "lm",color = "black")+
  theme_classic()+
  labs(title = 'Fungi EpT vs EpG',x = "Mean Exon per Transcript", y = "Mean Exon per Gene",
       subtitle = paste0('R2 = ',format(Fungi_fit_EpT_EpG$r.squared,digits = 3),' p = ',format(Fungi_fit_EpT_EpG$coefficients[2,4],digits = 3)))

plot_list <- list(Fungi_plot_TpG_EpT, Fungi_plot_TpG_EpG, Fungi_plot_EpT_EpG)

pdf('Sumpplimental_Fungi_Scatter_plots.pdf', 15, 15)

print(cowplot::plot_grid(plotlist = plot_list, nrow = 3, ncol = 1))
dev.off()
