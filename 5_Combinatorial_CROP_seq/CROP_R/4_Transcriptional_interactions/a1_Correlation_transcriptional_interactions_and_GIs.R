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
#--------------------- I. Scatter plot of GI scores and transcriptomic interactions in S2, S3 and S4----------------------------
rnaGI_S2 =read.csv("CROP_Python/sequencing/outs_S2/results/Relationship_of_GI_score_and_Transcriptomic_differences_S2.csv", row.names = 1)
rnaGI_S3 =read.csv("CROP_Python/sequencing/outs_S3/results/Relationship_of_GI_score_and_Transcriptomic_differences_S3.csv", row.names = 1)
rnaGI_S4 =read.csv("CROP_Python/sequencing/outs_S4/results/Relationship_of_GI_score_and_Transcriptomic_differences_S4.csv", row.names = 1)

##ggscater for S2
ggscatter(rnaGI_S2, x = "GI_score" , y = "transcriptomic_diff_var_lm", alpha = 0.8, size = 4,
          #palette = c("grey58", "black"),
          #label = "perturbed",
          #label.select = label.select,
          #font.label = c(16, "plain", "blue"),
          #label.rectangle = TRUE,
          repel = TRUE,
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "#007AFF", fill = "lightgray"), # Customize reg. line
          cor.coef = TRUE,rug = TRUE,
          cor.coeff.args = list(method = "pearson", label.x = -0.09, label.y =0.027, label.sep = "\n", size = 8))+
  theme(#panel.border = element_blank(), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    #axis.ticks.y = element_line(colour = "black"),
    legend.position = 'none',
    legend.text = element_text(size=21, color = 'black'),
    axis.text.y  =element_text(size=28,color = 'black'),
    axis.text.x =element_text(size=28,color = 'black'),
    axis.title=element_text(size=28,color = 'black'),#face="bold")
    plot.title = element_text(size=28, hjust= 0.5,color = 'black')
  )+
  scale_x_continuous(limits = c(-0.1,0.06), breaks = seq(-0.1,0.1,0.05))+
  scale_y_continuous(limits = c(-0.01,0.03), breaks = seq(-0.01,0.03,0.01))+
  labs(title="Full",x = "GI scores", y = "Transcriptomic interactions")
ggsave("CROP_R/4_Transcriptional_interactions/correlation_GI_and_transcription_interactions_S2.pdf", width = 7, height =5.5)

##ggscater for S3
ggscatter(rnaGI_S3, x = "GI_score" , y = "transcriptomic_diff_var_lm", alpha = 0.8, size = 4,
          #palette = c("grey58", "black"),
          #label = "perturbed",
          #label.select = label.select,
          #font.label = c(16, "plain", "blue"),
          #label.rectangle = TRUE,
          repel =TRUE,
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "#007AFF", fill = "lightgray"), # Customize reg. line
          cor.coef = TRUE,rug = TRUE,
          cor.coeff.args = list(method = "pearson", label.x =  -0.25, label.y =0.032, label.sep = "\n", size = 8))+
  theme(#panel.border = element_blank(), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    #axis.ticks.y = element_line(colour = "black"),
    legend.position = 'none',
    legend.text = element_text(size=21,color = 'black'),
    axis.text.y  =element_text(size=28,color = 'black'),
    axis.text.x =element_text(size=28,color = 'black'),
    axis.title=element_text(size=28, color = 'black'),#face="bold")
    plot.title = element_text(size=28,hjust= 0.5,color = 'black')
  )+
  scale_x_continuous(limits = c(-0.25,0.18), breaks = seq(-0.3,0.2,0.15))+
  scale_y_continuous(limits = c(-0.01,0.035), breaks = seq(-0.01,0.03,0.01))+
  labs(title="TGF-β1",x = "GI scores", y = "Transcriptomic interactions")
ggsave("CROP_R/4_Transcriptional_interactions/correlation_GI_and_transcription_interactions_S3.pdf", width = 7, height =5.5, device = cairo_pdf)

##ggscater for S4
rnaGI_S4$hL = "No"
rnaGI_S4["NF2-PTEN", 'hL'] = "Yes"
rnaGI_S4["NF2-TP53", 'hL'] = "Yes"
rnaGI_S4["PTEN-TP53", 'hL'] = "Yes"
ggscatter(rnaGI_S4, x = "GI_score" , y = "transcriptomic_diff_var_lm",  color = 'hL', alpha =0.8, size = 4,
          palette = c("grey58", "blue"),
          label = "perturbed",
          label.select = c("NF2-PTEN", "NF2-TP53", "PTEN-TP53"),
          font.label = c(20, "plain", "blue"),
          label.rectangle = FALSE,
          repel = TRUE,
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "#007AFF", fill = "lightgray"), # Customize reg. line
          cor.coef = TRUE,rug = TRUE,
          cor.coeff.args = list(method = "pearson", label.x = -0.1, label.y =0.018,label.sep = "\n", size = 7))+
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
  scale_x_continuous(limits = c(-0.1,0.2), breaks = seq(-0.1,0.2,0.1))+
  scale_y_continuous(limits = c(-0.01,0.02), breaks = seq(-0.01,0.02,0.01))+
  labs(title="Minimal",x = "GI scores", y = "Transcriptomic interactions")
ggsave("CROP_R/4_Transcriptional_interactions/correlation_GI_and_transcription_interactions_S4.pdf", width = 7, height =5.5)

#------------------------------------II. Distribution of transcriptomic differences in DKOs----------------------------------------
tI_1 =read.csv("CROP_Python/sequencing/outs_S4/results/prediction_and_transcriptomic_differences/diff_NF2-PTEN.csv", row.names = 1)
tI_2 =read.csv("CROP_Python/sequencing/outs_S4/results/prediction_and_transcriptomic_differences/diff_NF2-TP53.csv", row.names = 1)
tI_3 =read.csv("CROP_Python/sequencing/outs_S4/results/prediction_and_transcriptomic_differences/diff_PTEN-TP53.csv", row.names = 1)

tI_4 =read.csv("CROP_Python/sequencing/outs_S4/results/prediction_and_transcriptomic_differences/diff_CDH1-TP53.csv", row.names = 1)
tI_5 =read.csv("CROP_Python/sequencing/outs_S4/results/prediction_and_transcriptomic_differences/diff_CDH1-RB1.csv", row.names = 1)
tI_6 =read.csv("CROP_Python/sequencing/outs_S4/results/prediction_and_transcriptomic_differences/diff_CBFB-CDH1.csv", row.names = 1)
#tI_7 =read.csv("CROP_Python/sequencing/outs_S4/results/prediction_and_transcriptomic_differences/diff_NF1-TP53.csv", row.names = 1)

#Group 1:NF2, PTEN, TP53 - GI with high GI score and transcriptomic interactions
#data = rbind.data.frame(tI_1, tI_2, tI_3)
#Dual_sgRNAs = factor(data$name, levels = c("NF2-TP53", "NF2-PTEN","PTEN-TP53"))

#Group 2:NF2, PTEN, SMAD4, CBFB-TP53 - GI with low GI score and transcriptomic interactions
data = rbind.data.frame(tI_4, tI_5, tI_6)
Dual_sgRNAs = factor(data$name, levels = c("CDH1-TP53", "CDH1-RB1","CBFB-CDH1"))

ggplot(data, aes(x=diff_lm, color = Dual_sgRNAs)) +
  geom_density(alpha=.6, position="identity",na.rm = TRUE, size =2)+
  scale_fill_manual(values=c("blue", "green", "red"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line = element_line(colour = "black"),
    legend.position = c(0.77,0.75),
    legend.title = element_text(size = 18, colour = "black"),
    legend.text = element_text(size=18, color = "black" ),
    axis.text.y  =element_text(size=23, color = "black"),
    axis.text.x =element_text(size=24, color = "black"),
    axis.title=element_text(size=24),
    plot.title = element_text(size=36,hjust=0.5)
  )+
  labs(title="",x ="LFCs(DKO) - LFCs(fit)", y = "Density")+
  scale_x_continuous(limits= c(-0.2,0.2))+
  scale_y_continuous(limits= c(0,20), breaks= seq(0, 20, 5))+
  #scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
  #scale_fill_manual(values=c("red", "blue", "green3"))+
  geom_vline(xintercept = 0, linetype="dotted", color = "black", size=1)

#ggsave("CROP_R/4_Transcriptional_interactions/Transcriptomic_differences_NF2_PTEN_TP53.pdf", width = 5.6, height =5)
ggsave("CROP_R/4_Transcriptional_interactions/Transcriptomic_differences_CDH1_RB1_CBFB_TP53_scientific_y.pdf", width = 5.8, height =5)




