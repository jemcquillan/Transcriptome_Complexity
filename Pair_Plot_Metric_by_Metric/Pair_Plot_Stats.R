library("ggpubr")
library('ggplot2')

setwd("/Volumes/GoogleDrive/My Drive/OrthoDB/Ortholog_Final_Files/")
#setwd("set your working directory here")

# Input Data Converted to Dataframes
Deut_data <- read.csv("Deuterostome_Complexity/Deuterostome_merged_complexity/Deuterostome_WT_gtf_species_full_list.csv")
Dros_data <- read.csv("Drosophila_Complexity/Drosophila_merged_complexity/Drosophila_WT_gtf_species_full_list.csv")
Plant_data <- read.csv("Plantae_Complexity/Plantae_merged_complexity/Plantae_WT_gtf_species_full_list.csv")
Fungi_data <- read.csv("Fungi_Complexity/Fungi_merged_complexity/Fungi_WT_gtf_species_full_list.csv")
# Fungi_data is not relevant here because of annotation (i.e. TpG is all 1, and is a poor measure currently with the data available.)

# Run ANOVA Deuterostomes
Deut_species <- Deut_data$species
meanTpG <- Deut_data$mean_transcriptPerGene

Deut_one_way <- aov(mean_transcriptPerGene ~ Deut_species, data = Deut_data)
summary(Deut_one_way)

# Run ANOVA Drosophila
Dros_species <- Dros_data$species
Dros_meanTpG <- Dros_data$mean_transcriptPerGene

Dros_one_way <- aov(mean_transcriptPerGene ~ Dros_species, data = Dros_data)
summary(Dros_one_way)

# Run ANOVA Plantae
Plant_species <- Plant_data$species
Plant_meanTpG <- Plant_data$mean_transcriptPerGene

Plant_one_way <- aov(mean_transcriptPerGene ~ Plant_species, data = Plant_data)
summary(Plant_one_way)

# Run ANOVA Fungi
Fungi_species <- Fungi_data$species
Fungi_meanTpG <- Fungi_data$mean_transcriptPerGene

Fungi_one_way <- aov(mean_transcriptPerGene ~ Fungi_species, data = Fungi_data)
summary(Fungi_one_way)

# Tests of normality of the data

shapiro.test(Deut_data$mean_exonPerTranscript)
shapiro.test(Dros_data$mean_exonPerTranscript)
shapiro.test(Plant_data$mean_exonPerTranscript)
shapiro.test(Fungi_data$mean_exonPerTranscript)

###Create correlation plots with stats
all_data <- list(Deut_data,Dros_data,Plant_data)
plot_list <- list()
index <- 1

for (data in all_data){
  
  fit <- summary(lm(mean_exonPerTranscript~mean_transcriptPerGene,data = data))
  plot1<-ggplot(data, aes_string(x = "mean_transcriptPerGene", y = "mean_exonPerTranscript"))+
    geom_point(aes_string(color = "group"))+
    geom_smooth(method = "lm",color = "black")+
    theme_classic(base_size = 20)+
    labs(title = 'TpG vs EpT',x = "Mean Transcript per Gene", y = "Mean Exon per Transcript",
         subtitle = paste0('R2 = ',format(fit$r.squared,digits = 3),' p = ',format(fit$coefficients[2,4],digits = 3)))
  
  plot_list[[index]] <- plot1
  index <- index + 1
  
  fit <- summary(lm(mean_exonPerGene~mean_transcriptPerGene,data = data))
  plot2<-ggplot(data, aes_string(x = "mean_transcriptPerGene", y = "mean_exonPerGene"))+
    geom_point(aes_string(color = "group"))+
    geom_smooth(method = "lm",color = "black")+
    theme_classic(base_size = 20)+
    labs(title = 'TpG vs EpG',x= "Mean Transcript per Gene", y = "Mean Exon per Gene",
         subtitle = paste0('R2 = ',format(fit$r.squared,digits = 3),' p = ',format(fit$coefficients[2,4],digits = 3)))
  
  plot_list[[index]] <- plot2
  index <- index + 1
  
  fit <- summary(lm(mean_exonPerGene~mean_exonPerTranscript,data = data))
  plot3<-ggplot(data, aes_string(x = "mean_exonPerTranscript", y = "mean_exonPerGene"))+
    geom_point(aes_string(color = "group"))+
    geom_smooth(method = "lm",color = "black")+
    theme_classic2(base_size = 20)+
    labs(title = 'EpT vs EpG',x = "Mean Exon per Transcript", y = "Mean Exon per Gene",
         subtitle = paste0('R2 = ',format(fit$r.squared,digits = 3),' p = ',format(fit$coefficients[2,4],digits = 3)))
  
  plot_list[[index]] <- plot3
  index <- index + 1
  
}


pdf('Fig2.Correlation_Pair_Plot.pdf',20,20)

print(cowplot::plot_grid(plotlist = plot_list, nrow = 3, ncol = 3))
dev.off()
