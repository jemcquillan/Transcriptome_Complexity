library(tidyverse)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(stringr)

setwd('set your working directory here - /Ortholog_Final_Files/Metric_Analysis/')
output_directory <- "Set your output directory here"


complexity_metrics <- read.csv("Complexity_Data/All_Complexity_metrics.csv")

complexity_metrics_plot <- complexity_metrics %>% 
  tidyr::gather(key = 'group', value = 'variation',TpG.Normalized.Variation,EpT.Normalized.Variation,EpG.Normalized.Variation) 

complexity_metrics_plot$group <- stringr::str_remove(complexity_metrics_plot$group,'\\.Normalized\\.Variation')
complexity_metrics_plot$group <- factor(complexity_metrics_plot$group,levels = c('TpG','EpT','EpG'))
complexity_metrics_plot$Higher.Order <- factor(complexity_metrics_plot$Higher.Order,levels = c('Fungi','Plantae','Flies','Deuterostomia'))

ggplot(complexity_metrics_plot, aes(x = group, y = variation)) +
  geom_boxplot(aes(color = Higher.Order, fill = Higher.Order), alpha = 0.5) +
  ggtitle("Normalized Variance Complexity Metrics") +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), aes(color = Higher.Order), alpha = 0.6, size = 0.7)+
  scale_color_manual(values = c("#D95F02","#1B9E77","#7570B3","#E7298A"))+
  scale_fill_manual(values = c("#D95F02","#1B9E77","#7570B3","#E7298A"))+
  labs(y = 'Normalized variation',x ='', color = "Higher Order", fill = "Higher Order") +
  theme_classic()


# Filter the groups
Deuterostomes <- dplyr::filter(complexity_metrics, Higher.Order %in% "Deuterostomia")
Fungus <- dplyr::filter(complexity_metrics, Higher.Order %in% "Fungi")
Drosophila <- dplyr::filter(complexity_metrics, Higher.Order %in% "Flies")
Plants <- dplyr::filter(complexity_metrics, Higher.Order %in% "Plantae")


# Within group normalized variance EpT vs EpG, Wilcoxon test
wilcox.test(Deuterostomes$EpT.Normalized.Variation,Deuterostomes$EpG.Normalized.Variation)
wilcox.test(Fungus$EpT.Normalized.Variation,Fungus$EpG.Normalized.Variation)
wilcox.test(Drosophila$EpT.Normalized.Variation,Drosophila$EpG.Normalized.Variation)
wilcox.test(Plants$EpT.Normalized.Variation,Plants$EpG.Normalized.Variation)

# Within group mean EpT vs EpG, Wilcoxon test
wilcox.test(Deuterostomes$mean_exonPerGene,Deuterostomes$mean_exonPerTranscript)
wilcox.test(Fungus$mean_exonPerGene,Fungus$mean_exonPerTranscript)
wilcox.test(Drosophila$mean_exonPerGene,Drosophila$mean_exonPerTranscript)
wilcox.test(Plants$mean_exonPerGene,Plants$mean_exonPerTranscript)


pdf(paste0(output_directory,'/Mean_Complexity_Metrics_Boxplots.pdf'),8,6)
#Kruskal-Wallis Test
#All group TpG
test <- kruskal.test(mean_transcriptPerGene ~ Higher.Order, data = complexity_metrics)
boxplot(complexity_metrics$mean_transcriptPerGene ~ factor(complexity_metrics$Higher.Order), main = paste0("All Groups Mean TpG\nKruskal-Wallis rank sum test\nChi-Sq = ",format(test$statistic,digits = 4)," p= ",format(test$p.value,digits = 4)), cex.axis = 0.7, las = 2, xlab="", ylab="Mean")
#All group EpT
test <- kruskal.test(complexity_metrics$mean_exonPerTranscript ~ Higher.Order, data = complexity_metrics)
boxplot(complexity_metrics$mean_exonPerTranscript ~ factor(complexity_metrics$Higher.Order), main = paste0("All Groups Mean EpT\nKruskal-Wallis rank sum test\nChi-Sq = ",format(test$statistic,digits = 4)," p= ",format(test$p.value,digits = 4)), cex.axis = 0.7, las = 2, xlab="", ylab="Mean")
#All group EpG
test <- kruskal.test(complexity_metrics$mean_exonPerGene ~ Higher.Order, data = complexity_metrics)
boxplot(complexity_metrics$mean_exonPerGene ~ factor(complexity_metrics$Higher.Order), main = paste0("All Groups Mean EpG\nKruskal-Wallis rank sum test\nChi-Sq = ",format(test$statistic,digits = 4)," p= ",format(test$p.value,digits = 4)), cex.axis = 0.7, las = 2, xlab="", ylab="Mean")

#Deuterostome Class TpG
test <- kruskal.test(Deuterostomes$mean_transcriptPerGene ~ group, data = Deuterostomes)
boxplot(Deuterostomes$mean_transcriptPerGene ~ factor(Deuterostomes$group), main = paste0("Deuterostome Class Mean TpG\nKruskal-Wallis rank sum test\nChi-Sq = ",format(test$statistic,digits = 4)," p= ",format(test$p.value,digits = 4)), cex.axis = 0.65, las = 2, xlab="", ylab="Mean")
#Deuterostome Class EpT
test <- kruskal.test(Deuterostomes$mean_exonPerTranscript ~ group, data = Deuterostomes)
boxplot(Deuterostomes$mean_exonPerTranscript ~ factor(Deuterostomes$group), main = paste0("Deuterostome Class Mean EpT\nKruskal-Wallis rank sum test\nChi-Sq = ",format(test$statistic,digits = 4)," p= ",format(test$p.value,digits = 4)), cex.axis = 0.65, las = 2, xlab="", ylab="Mean")
#Deuterostome Class EpG
test <- kruskal.test(Deuterostomes$mean_exonPerGene ~ group, data = Deuterostomes)
boxplot(Deuterostomes$mean_exonPerGene ~ factor(Deuterostomes$group), main = paste0("Deuterostome Class Mean EpG\nKruskal-Wallis rank sum test\nChi-Sq = ",format(test$statistic,digits = 4)," p= ",format(test$p.value,digits = 4)), cex.axis = 0.65, las = 2, xlab="", ylab="Mean")

#Drosophila Class TpG
shapiro.test(Drosophila$mean_transcriptPerGene)
test <- t.test(Drosophila$mean_transcriptPerGene, data = Drosophila)
boxplot(Drosophila$mean_transcriptPerGene ~ factor(Drosophila$species), main = paste0("Drosophila Class Mean TpG\nKruskal-Wallis rank sum test\nChi-Sq = ",format(test$statistic,digits = 4)," p= ",format(test$p.value,digits = 4)), cex.axis = 0.65, las = 2, xlab="", ylab="Mean")
#Drosophila Class EpT
shapiro.test(Drosophila$mean_exonPerTranscript)
test <- wilcox.test(Drosophila$mean_exonPerTranscript, data = Drosophila)
boxplot(Drosophila$mean_exonPerTranscript ~ factor(Drosophila$species), main = paste0("Drosophila Class Mean EpT\nKruskal-Wallis rank sum test\nChi-Sq = ",format(test$statistic,digits = 4)," p= ",format(test$p.value,digits = 4)), cex.axis = 0.65, las = 2, xlab="", ylab="Mean")
#Drosophila Class EpG
shapiro.test(Drosophila$mean_exonPerGene)
test <- wilcox.test(Drosophila$mean_exonPerGene, data = Drosophila)
boxplot(Drosophila$mean_exonPerGene ~ factor(Drosophila$species), main = paste0("Drosophila Class Mean EpG\nKruskal-Wallis rank sum test\nChi-Sq = ",format(test$statistic,digits = 4)," p= ",format(test$p.value,digits = 4)), cex.axis = 0.65, las = 2, xlab="", ylab="Mean")


#Fungus Class TpG
test <- kruskal.test(Fungus$mean_transcriptPerGene ~ group, data = Fungus)
boxplot(Fungus$mean_transcriptPerGene ~ factor(Fungus$group), main = paste0("Fungus Class Mean TpG\nKruskal-Wallis rank sum test\nChi-Sq = ",format(test$statistic,digits = 4)," p= ",format(test$p.value,digits = 4)), cex.axis = 0.45, las = 2, xlab="", ylab="Mean")
#Fungus Class EpT
test <- kruskal.test(Fungus$mean_exonPerTranscript ~ group, data = Fungus)
boxplot(Fungus$mean_exonPerTranscript ~ factor(Fungus$group), main = paste0("Fungus Class Mean EpT\nKruskal-Wallis rank sum test\nChi-Sq = ",format(test$statistic,digits = 4)," p= ",format(test$p.value,digits = 4)), cex.axis = 0.45, las = 2, xlab="", ylab="Mean")
#Fungus Class EpG
test <- kruskal.test(Fungus$mean_exonPerGene ~ group, data = Fungus)
boxplot(Fungus$mean_exonPerGene ~ factor(Fungus$group), main = paste0("Fungus Class Mean EpG\nKruskal-Wallis rank sum test\nChi-Sq = ",format(test$statistic,digits = 4)," p= ",format(test$p.value,digits = 4)), cex.axis = 0.45, las = 2, xlab="", ylab="Mean")

#Fungus Phylum TpG
test <- kruskal.test(Fungus$mean_transcriptPerGene ~ Phylum, data = Fungus)
boxplot(Fungus$mean_transcriptPerGene ~ factor(Fungus$Phylum), main = paste0("Fungus Phylum Mean TpG\nKruskal-Wallis rank sum test\nChi-Sq = ",format(test$statistic,digits = 4)," p= ",format(test$p.value,digits = 4)), cex.axis = 0.68, las = 2, xlab="", ylab="Mean")
#Fungus Phylum EpT
test <- kruskal.test(Fungus$mean_exonPerTranscript ~ Phylum, data = Fungus)
boxplot(Fungus$mean_exonPerTranscript ~ factor(Fungus$Phylum), main = paste0("Fungus Phylum Mean EpT\nKruskal-Wallis rank sum test\nChi-Sq = ",format(test$statistic,digits = 4)," p= ",format(test$p.value,digits = 4)), cex.axis = 0.68, las = 2, xlab="", ylab="Mean")
#Fungus Phylum EpG
test <- kruskal.test(Fungus$mean_exonPerGene ~ Phylum, data = Fungus)
boxplot(Fungus$mean_exonPerGene ~ factor(Fungus$Phylum), main = paste0("Fungus Phylum Mean EpG\nKruskal-Wallis rank sum test\nChi-Sq = ",format(test$statistic,digits = 4)," p= ",format(test$p.value,digits = 4)), cex.axis = 0.68, las = 2, xlab="", ylab="Mean")

#Fungus Phylum EpT sans Ascomycota
fungi_sans_ascomycota <- Fungus[Fungus$Phylum != "Ascomycota", ]
test <- kruskal.test(fungi_sans_ascomycota$mean_exonPerTranscript ~ Phylum, data = fungi_sans_ascomycota)
boxplot(fungi_sans_ascomycota$mean_exonPerGene ~ factor(fungi_sans_ascomycota$Phylum), main = paste0("Fungus Phylum sans Ascomycota Mean EpT\nKruskal-Wallis rank sum test\nChi-Sq = ",format(test$statistic,digits = 4)," p= ",format(test$p.value,digits = 4)), cex.axis = 0.68, las = 2, xlab="", ylab="Mean")

#Fungus Phylum EpG sans Ascomycota
fungi_sans_ascomycota <- Fungus[Fungus$Phylum != "Ascomycota", ]
test <- kruskal.test(fungi_sans_ascomycota$mean_exonPerGene ~ Phylum, data = fungi_sans_ascomycota)
boxplot(fungi_sans_ascomycota$mean_exonPerGene ~ factor(fungi_sans_ascomycota$Phylum), main = paste0("Fungus Phylum sans Ascomycota Mean EpG\nKruskal-Wallis rank sum test\nChi-Sq = ",format(test$statistic,digits = 4)," p= ",format(test$p.value,digits = 4)), cex.axis = 0.68, las = 2, xlab="", ylab="Mean")


#Plants Order TpG
test <- kruskal.test(Plants$mean_transcriptPerGene ~ group, data = Plants)
boxplot(Plants$mean_transcriptPerGene ~ factor(Plants$group), main = paste0("Plants Order Mean TpG\nKruskal-Wallis rank sum test\nChi-Sq = ",format(test$statistic,digits = 4)," p= ",format(test$p.value,digits = 4)), cex.axis = 0.72, las = 2, xlab="", ylab="Mean")
#Plants Order EpT
test <- kruskal.test(Plants$mean_exonPerTranscript ~ group, data = Plants)
boxplot(Plants$mean_exonPerTranscript ~ factor(Plants$group), main = paste0("Plants Order Mean EpT\nKruskal-Wallis rank sum test\nChi-Sq = ",format(test$statistic,digits = 4)," p= ",format(test$p.value,digits = 4)), cex.axis = 0.72, las = 2, xlab="", ylab="Mean")
#Plants Order EpG
test <- kruskal.test(Plants$mean_exonPerGene ~ group, data = Plants)
boxplot(Plants$mean_exonPerGene ~ factor(Plants$group), main = paste0("Plants Order Mean EpG\nKruskal-Wallis rank sum test\nChi-Sq = ",format(test$statistic,digits = 4)," p= ",format(test$p.value,digits = 4)), cex.axis = 0.72, las = 2, xlab="", ylab="Mean")

dev.off()


