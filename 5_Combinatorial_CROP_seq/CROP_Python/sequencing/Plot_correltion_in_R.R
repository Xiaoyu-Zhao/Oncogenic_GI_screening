#Xiaoyu Zhao 20200908
library(ggplot2)
library(data.table) ## for read.table("xxx.xls",sep = "\t",header = T) 
library(pheatmap) 
library(Cairo) #https://stackoverflow.com/questions/27690729/greek-letters-symbols-and-line-breaks-inside-a-ggplot-legend-label
library(VennDiagram)
library(ggpubr)
library(tidyr)
library(UpSetR)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#--------------------- I. Correlation of fitness and number of DEGs in S2----------------------------
setwd("~/Desktop/Paper_writing_XYZ_20200908/5_CROPseq_20200908/Python_code/sequencing/")
df2 = read.csv('outs_S2/results/correlation_fitness_and_transcription/correlation_fitness_and_transcription_S2.csv')
dim(df2)
df2['log2_de'] = log2(df2['DESeq2_de']+1)
df2_SKO = df2[df2$knock_outs == "SKO",]

##ggscater for df2 or df2_SKO
ggscatter(df2_SKO, x = "log2_de" , y = "fitness", size = 3,
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "#007AFF", fill = "lightgray"), # Customize reg. line
          cor.coef = TRUE,rug = FALSE,
          cor.coeff.args = list(method = "pearson", label.x = 0.5,label.y =0.3, label.sep = "\n", size = 8))+
  theme(#panel.border = element_blank(), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    #axis.ticks.y = element_line(colour = "black"),
    legend.position = 'none',
    #legend.position = c(0.75,0.9),
    #legend.title = element_text(size = 20, face= "bold"),
    legend.text = element_text(size=21),
    axis.text.y  =element_text(size=28,color = 'black'),
    axis.text.x =element_text(size=28,color = 'black', margin = margin(5,0,0,0)),
    axis.title=element_text(size=28,margin = margin(10,0,0,0)),#face="bold")
    plot.title = element_text(size=23, hjust= 1, margin = margin(0,0,10,0))
    #strip.text  = element_text(size = 16, face="bold")
    ##strip.text is for facet
  )+
  scale_x_continuous(limits = c(0,9), breaks = seq(0,8,2))+
  scale_y_continuous(limits = c(-0.2,0.4), breaks = seq(-0.2,0.4,0.2))+
  labs(title="",x = expression(log[2]~(no.~of~DEGs)), y = "Fitness (SKOs) in S2")
ggsave("outs_S2/results/correlation_fitness_and_transcription/correlation_fitness_and_transcription_SKOs_S2.pdf", width = 4.8, height =4.2)
#ggsave("outs_S2/results/correlation_fitness_and_transcription/correlation_fitness_and_transcription_SKOs_S2.pdf", width = 4.8, height =4.4)

#--------------------- II. Correlation of fitness and number of DEGs in S3----------------------------
df3 = read.csv('outs_S3/results/correlation_fitness_and_transcription/correlation_fitness_and_transcription_S3.csv')
dim(df3)
df3['log2_de'] = log2(df3['DESeq2_de']+1)
df3_SKO = df3[df3$knock_outs == "SKO",]

##ggscater for df3 or df3_SKO
ggscatter(df3_SKO, x = "log2_de" , y = "fitness", size = 3,
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "#007AFF", fill = "lightgray"), # Customize reg. line
          cor.coef = TRUE,rug = FALSE,
          cor.coeff.args = list(method = "pearson", label.x = 0.5,label.y =0.3, label.sep = "\n", size = 8))+
  theme(#panel.border = element_blank(), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    #axis.ticks.y = element_line(colour = "black"),
    legend.position = 'none',
    #legend.position = c(0.75,0.9),
    #legend.title = element_text(size = 20, face= "bold"),
    legend.text = element_text(size=21),
    axis.text.y  =element_text(size=28,color = 'black'),
    axis.text.x =element_text(size=28,color = 'black', margin = margin(5,0,0,0)),
    axis.title=element_text(size=28,margin = margin(10,0,0,0)),#face="bold")
    plot.title = element_text(size=23, hjust= 1, margin = margin(0,0,10,0))
    #strip.text  = element_text(size = 16, face="bold")
    ##strip.text is for facet
  )+
  scale_x_continuous(limits = c(0,9), breaks = seq(0,8,2))+
  scale_y_continuous(limits = c(-0.2,0.4), breaks = seq(-0.2,0.4,0.2))+
  labs(title="",x = expression(log[2]~(no.~of~DEGs)), y = "Fitness (SKOs) in S3")
ggsave("outs_S3/results/correlation_fitness_and_transcription/correlation_fitness_and_transcription_SKOs_S3.pdf", width = 4.8, height = 4.2)
#ggsave("outs_S3/results/correlation_fitness_and_transcription/correlation_fitness_and_transcription_SKOs_S3.pdf", width = 4.8, height = 4.4)

#--------------------- III. Correlation of fitness and number of DEGs in S4----------------------------
df4 = read.csv('outs_S4/results/correlation_fitness_and_transcription/correlation_fitness_and_transcription_S4.csv')
df4['log2_de'] = log2(df4['DESeq2_de']+1)
df4_SKO = df4[df4$knock_outs == "SKO",]

##ggscater for df4 or df4_SKO
ggscatter(df4_SKO, x = "log2_de" , y = "fitness", size = 3,
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "#007AFF", fill = "lightgray"), # Customize reg. line
          cor.coef = TRUE,rug = FALSE,
          cor.coeff.args = list(method = "pearson", label.x = 4.5,label.y =-0.1, label.sep = "\n", size = 8))+
  theme(#panel.border = element_blank(), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    #axis.ticks.y = element_line(colour = "black"),
    legend.position = 'none',
    #legend.position = c(0.75,0.9),
    #legend.title = element_text(size = 20, face= "bold"),
    legend.text = element_text(size=21),
    axis.text.y  =element_text(size=28,color = 'black'),
    axis.text.x =element_text(size=28,color = 'black', margin = margin(5,0,0,0)),
    axis.title=element_text(size=28,margin = margin(10,0,0,0)),#face="bold")
    plot.title = element_text(size=23, hjust= 1, margin = margin(0,0,10,0))
    #strip.text  = element_text(size = 16, face="bold")
    ##strip.text is for facet
  )+
  scale_x_continuous(limits = c(0,8), breaks = seq(0,8,2))+
  scale_y_continuous(limits = c(-0.2,0.4), breaks = seq(-0.2,0.4,0.2))+
  labs(title="",x = expression(log[2]~(no.~of~DEGs)), y = "Fitness (SKOs) in S4")
ggsave("outs_S4/results/correlation_fitness_and_transcription/correlation_fitness_and_transcription_SKOs_S4.pdf", width = 4.8, height =4.2)





#--------------------- IV. Number of DEGs induced by SKOs in S2, S3 and S4----------------------------
num_DEGs = read.csv('Results_S1234/num_DEGs_20200908.csv')
num_DEGs[num_DEGs == 0] = NA
SGS = factor(num_DEGs$SKOs, levels= c("NF2-CTRL1",  "PTEN-CTRL1","SMAD4-CTRL1", "CBFB-CTRL1", "NF1-CTRL1",  "TP53-CTRL1", "RB1-CTRL1", "CDH1-CTRL1", "CASP8-CTRL1", "TBX3-CTRL1", "USP9X-CTRL1"))
SGS = factor(num_DEGs$SKOs, levels = rev(levels(SGS)))
Conditions = factor(num_DEGs$conditions, levels = c("S2_full","S3_tgf","S4_min" ))
ggplot(num_DEGs, aes(Conditions, SGS, size= DESeq2_de))+ 
  geom_point(color = "blue")+ scale_size_area(max_size = 15)+
  #scale_fill_gradient(color = "blue")+
  theme(#panel.border = element_blank(), 
    #panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.grid= element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none",
    legend.title = element_text(size = 20, face= "bold"),
    legend.text = element_text(size=20 ),
    axis.text.y  = element_text(size=18,  color = 'black'),
    axis.text.x =element_text(size=20, color = 'black',vjust=0),
    axis.title.y=element_text(size=22, color = 'black', vjust=-1.5, margin = margin(t =0, r =0, b =1, l=0, unit = "cm")),#face="bold")
    plot.title = element_text(size=26,hjust= 0.5, vjust = 2)
    #strip.text.y  = element_text(size = 26, face="bold")
    ##strip.text is for facet
  )+
  scale_y_discrete(labels = c("sgUSP9X","sgTBX3","sgCASP8","sgCDH1", "sgRB1","sgTP53", "sgNF1", "sgCBFB", "sgSMAD4", "sgPTEN","sgNF2"))+
  labs(title="No. of DEGs",x ="", y = "SKOs")
ggsave("Results_S1234/Num_DEGs_SKOs_20200908.pdf", width = 6., height = 6.5)








