#Xiaoyu Zhao 
library(ggplot2)
library(data.table) ## for read.table("xxx.xls",sep = "\t",header = T) 
library(pheatmap) 
library(Cairo) #https://stackoverflow.com/questions/27690729/greek-letters-symbols-and-line-breaks-inside-a-ggplot-legend-label
library(VennDiagram)
library(ggpubr)
library(tidyr)
library(UpSetR)
library(dplyr)
library(scales)
setwd("~/Downloads/5_Combinatorial_CROP_seq/")
rnaGI_S2 =read.csv("CROP_Python/sequencing/outs_S2/results/Relationship_of_GI_score_and_Transcriptomic_differences_S2.csv", row.names = 1)
rnaGI_S3 =read.csv("CROP_Python/sequencing/outs_S3/results/Relationship_of_GI_score_and_Transcriptomic_differences_S3.csv", row.names = 1)
rnaGI_S4 =read.csv("CROP_Python/sequencing/outs_S4/results/Relationship_of_GI_score_and_Transcriptomic_differences_S4.csv", row.names = 1)

nm = unique(c(rownames(rnaGI_S2), rownames(rnaGI_S3), rownames(rnaGI_S4)))
tIs = as.data.frame(matrix(NA, length(nm), 6))
rownames(tIs) = nm
colnames(tIs) = c("GI_S2","tI_S2", "GI_S3","tI_S3","GI_S4","tI_S4")

tIs[,1:2]= c(rnaGI_S2[match(nm, rownames(rnaGI_S2)),c("GI_score","transcriptomic_diff_var_lm")])
tIs[,3:4]= c(rnaGI_S3[match(nm, rownames(rnaGI_S3)),c("GI_score","transcriptomic_diff_var_lm")])
tIs[,5:6]= c(rnaGI_S4[match(nm, rownames(rnaGI_S4)),c("GI_score","transcriptomic_diff_var_lm")])

#--------------------- I. Relationship among transcriptomic interactions between S2 and S4----------------------------
ggscatter(tIs, x = "GI_S4" , y = "tI_S2",  alpha =0.8, size = 4,
          palette = c("grey58", "blue"),
          #label = "perturbed",
          #label.select = c("NF2-PTEN", "NF2-TP53", "PTEN-TP53"),
          font.label = c(20, "plain", "blue"),
          label.rectangle = FALSE,
          repel = TRUE,
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "#007AFF", fill = "lightgray"), # Customize reg. line
          cor.coef = TRUE,rug = TRUE,
          cor.coeff.args = list(method = "pearson", label.x = -0.1, label.y =0.028,label.sep = "\n", size = 7))+
  scale_size(range = c(0, 4.0))+
  theme(#panel.border = element_blank(), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    #axis.ticks.y = element_line(colour = "black"),
    legend.position = 'none',
    legend.text = element_text(size=21),
    axis.text.y  =element_text(size=28,color = 'black'),
    axis.text.x =element_text(size=28,color = 'black', margin = margin(5,0,0,0)),
    axis.title=element_text(size=28,margin = margin(10,0,0,0)),#face="bold")
    plot.title = element_text(size=28, hjust= 0.5, margin = margin(0,0,10,0))
  )+
  scale_x_continuous(limits = c(-0.1,0.18), breaks = seq(-0.1,0.2,0.1))+
  scale_y_continuous(limits = c(-0.002,0.03), breaks = seq(-0.01,0.03,0.01))+
  labs(title="",x = "GI score (Minimal)", y = "TI score (Full)")
ggsave("CROP_R/4_Transcriptional_interactions/correlation_S4_GI_and_S2_TI.pdf", width = 7, height =5.5)


#--------------------- II. Relationship among transcriptomic interactions between S2 and S4----------------------------
ggscatter(tIs, x = "GI_S3" , y = "tI_S2",  alpha =0.8, size = 4,
          palette = c("grey58", "blue"),
          #label = "perturbed",
          #label.select = c("NF2-PTEN", "NF2-TP53", "PTEN-TP53"),
          font.label = c(20, "plain", "blue"),
          label.rectangle = FALSE,
          repel = TRUE,
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "#007AFF", fill = "lightgray"), # Customize reg. line
          cor.coef = TRUE,rug = TRUE,
          cor.coeff.args = list(method = "pearson", label.x = -0.25, label.y =0.028,label.sep = "\n", size = 7))+
  scale_size(range = c(0, 4.0))+
  theme(#panel.border = element_blank(), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    #axis.ticks.y = element_line(colour = "black"),
    legend.position = 'none',
    legend.text = element_text(size=21),
    axis.text.y  =element_text(size=28,color = 'black'),
    axis.text.x =element_text(size=28,color = 'black', margin = margin(5,0,0,0)),
    axis.title=element_text(size=28,margin = margin(10,0,0,0)),#face="bold")
    plot.title = element_text(size=28, hjust= 0.5, margin = margin(0,0,10,0))
  )+
  scale_x_continuous(limits = c(-0.28,0.11), breaks = seq(-0.2,0.1,0.1))+
  scale_y_continuous(limits = c(-0.002,0.03), breaks = seq(-0.01,0.03,0.01))+
  labs(title="",x = "GI score (TGF-β1)", y = "TI score (Full)")
ggsave("CROP_R/4_Transcriptional_interactions/correlation_S3_GI_and_S2_TI.pdf", width = 7, height =5.5, device = cairo_pdf)
























#--------------------- III. Relationship among transcriptomic interactions between S3 and S4----------------------------
ggscatter(tIs, x = "tI_S3" , y = "tI_S4",  alpha =0.8, size = 4,
          palette = c("grey58", "blue"),
          #label = "perturbed",
          #label.select = c("NF2-PTEN", "NF2-TP53", "PTEN-TP53"),
          font.label = c(20, "plain", "blue"),
          label.rectangle = FALSE,
          repel = TRUE,
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "#007AFF", fill = "lightgray"), # Customize reg. line
          cor.coef = TRUE,rug = TRUE,
          cor.coeff.args = list(method = "pearson", label.x = 0.018, label.y =0.018,label.sep = "\n", size = 7))+
  scale_size(range = c(0, 4.0))+
  theme(#panel.border = element_blank(), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    #axis.ticks.y = element_line(colour = "black"),
    legend.position = 'none',
    legend.text = element_text(size=21),
    axis.text.y  =element_text(size=28,color = 'black'),
    axis.text.x =element_text(size=28,color = 'black', margin = margin(5,0,0,0)),
    axis.title=element_text(size=28,margin = margin(10,0,0,0)),#face="bold")
    plot.title = element_text(size=28, hjust= 0.5, margin = margin(0,0,10,0))
  )+
  scale_x_continuous(limits = c(-0.002,0.035), breaks = seq(-0.01,0.03,0.01))+
  scale_y_continuous(limits = c(-0.002,0.02), breaks = seq(-0.01,0.02,0.01))+
  labs(title="",x = "TI scores (TGF-β1)", y = "TI scores (Minimal)")
ggsave("CROP_R/4_Transcriptional_interactions/correlation_of_transcription_interactions_in_S3_and_S4.pdf", width = 7, height =5.5, device = cairo_pdf)

