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
tIs = read.csv("CROP_R/4_Transcriptional_interactions/Groups_of_TI_categorization_5069_genes_based_on_TDs_3DKOs.csv", row.names = 1)
#-------------------------------I. Common transcriptional differences-------------------------------------
#Generate matrices that contain (DKO_LFC and diff_lm) of common transcriptional differences
#Add_up
comAdd_up = tIs[which(tIs$Group_Additive_Up == "Additive_Up: Common"),]
write.csv(comAdd_up, "CROP_R/4_Transcriptional_interactions/74_common_additive_upregulated_genes.csv")
#Syn_up
comSyn_up = tIs[which(tIs$Group_Synergistic_Up == "Synergistic_Up: Common"),]
write.csv(comSyn_up, "CROP_R/4_Transcriptional_interactions/148_common_synergistic_upregulated_genes.csv")
#Buff_up
comBuff_up = tIs[which(tIs$Group_Buffering_Up == "Buffering_Up: Common"),]
write.csv(comBuff_up, "CROP_R/4_Transcriptional_interactions/1_common_buffering_upregulated_genes.csv")
#-------------------------------------------------------------------------------------------------------------
#Add_down
comAdd_dn = tIs[which(tIs$Group_Additive_Dn== "Additive_Dn: Common"),]
write.csv(comAdd_dn, "CROP_R/4_Transcriptional_interactions/28_common_additive_downregulated_genes.csv")
#Syn_down
comSyn_dn = tIs[which(tIs$Group_Synergistic_Dn == "Synergistic_Dn: Common"),]
write.csv(comSyn_dn, "CROP_R/4_Transcriptional_interactions/9_common_synergistic_downregulated_genes.csv")
#Buff_up
comBuff_dn = tIs[which(tIs$Group_Buffering_Dn == "Buffering_Dn: Common"),]
write.csv(comBuff_dn, "CROP_R/4_Transcriptional_interactions/2_common_buffering_downregulated_genes.csv")

#-------------------------------II. Pather Enrichment analysis of common synergistic upregulated tIs-------------------------------------
d1 = cbind.data.frame(Pathways = "Cytoskeletal regulation by Rho GTPase",
                      FDR = -log10(5.38e-4))

d2 = cbind.data.frame(Pathways = c("RNA splicing factor", "Acyltransferase", "Translational protein", "Cytoskeletal protein"),
                      FDR = c(-log10(6.31e-8), -log10(1.48e-2), -log10(1.63e-4), -log10(7.69e-3) ))
Pathways0 = factor(d2$Pathways, levels = c("Acyltransferase", "Cytoskeletal protein", "Translational protein",  "RNA splicing factor"))

ggplot(d2, aes(x = Pathways0, y = FDR, color = Pathways0))+
  geom_bar(stat="identity", width = 0.7, fill = "#5AC8FA", size=0.6, alpha = 0.6)+
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        axis.line = element_line(colour = "black"),
        # legend.position = c(2,0.9),
        legend.position = "None",
        legend.title = element_text(size = 26, color = "black"),
        legend.text = element_text(size=26 , color = "black"),
        axis.text.y  =element_text(size=16, color = "black"),
        axis.text.x =element_text(size=16,vjust =1, hjust = 1, color = "black", margin = margin(t = 0, r = 0, b = 0, l = 10)),
        axis.title=element_text(size=18),
        plot.title = element_text(size=18, hjust = 0.5)
        #strip.text.y  = element_text(size = 26, face="bold")
        ##strip.text is for facet
  )+
  coord_flip()+
  scale_y_continuous(position = "left")+
  scale_color_manual(values=c('limegreen','blue','orange','red'))+
  labs(title="",x ="Protein class", y = expression(-log[10]~FDR))

ggsave("CROP_R/4_Transcriptional_interactions/Enriched_pather_protein_class_Common_Syn_UP_tIs.pdf", width = 5.6, height = 2.6)

#---------------------------------------III. Figures on enriched_pather_protein_class------------------------------------------------------
##1. Enriched RNA splicing factors
splicing = c("SRSF2", "SRSF11", "SRSF3", "SRSF7", "SNRPE", "SNRPF", "SNRPG", "SNRPB", "LSM4", "LSM2", "LSM5","NAA38")
#splicing = c("SRSF2", "SRSF11", "SRSF3", "SRSF7", "SNRPE", "SNRPF", "SNRPG", "SNRPB")
splicing_M = comSyn_up[splicing, 2:7]
colnames(splicing_M) = c("NF2","PTEN", "TP53", "NF2-PTEN", "NF2-TP53", "PTEN-TP53")
d1 = melt(splicing_M)
d1$Group = c(rep("SKOs", nrow(d1)/2), rep("DKOs", nrow(d1)/2) )

vars = factor(d1$variable, levels = c("NF2", "PTEN", "TP53", "NF2-PTEN", "PTEN-TP53", "NF2-TP53"))
ggplot(d1, aes(x = vars, y = value, fill = Group))+
  geom_boxplot(width = 0.7)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.2), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line = element_line(colour = "black"),
    # legend.position = c(0.8,0.7),
    #legend.position = "None",
    legend.title = element_text(size = 16, color = "black"),
    legend.text = element_text(size=16, color = "black" ),
    axis.text.y  =element_text(size=20, color = "black"),
    axis.text.x =element_text(size=16, angle = 45, vjust =1,hjust =1, color = "black", margin = margin(t = 0, r = 10, b = 10, l = 0)),
    axis.title=element_text(size=20, color = "black"),
    plot.title = element_text(size=19,hjust=0.5, color = "black")
  )+
  labs(title="RNA splicing facor",x ="", y = "LFCs")
ggsave("CROP_R/4_Transcriptional_interactions/RNA_splicing_factor_enriched_in_3_DKOs.pdf", width = 4.2, height =4.6)
t.test(splicing_M$`NF2-TP53`, splicing_M$NF2,  alternative = "greater")

##2.Enriched Translation protein
translation = c("MRPL14", "EIF4EBP1", "MRPL20", "EIF2S2", "MRPL37", "EIF1AX", "MRPL18", "MRPL51", "ABCF1", "MRPS12", "MRPS18C")
translation_M = comSyn_up[translation, 2:7]
colnames(translation_M) = c("NF2","PTEN", "TP53", "NF2-PTEN", "NF2-TP53", "PTEN-TP53")
d2 = melt(translation_M)
d2$Group = c(rep("SKOs", nrow(d2)/2), rep("DKOs", nrow(d2)/2) )

vars = factor(d2$variable, levels = c("NF2", "PTEN", "TP53", "NF2-PTEN", "PTEN-TP53", "NF2-TP53"))
ggplot(d2, aes(x = vars, y = value, fill = Group))+
  geom_boxplot(width = 0.7)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.2), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line = element_line(colour = "black"),
    # legend.position = c(0.8,0.7),
    #legend.position = "None",
    legend.title = element_text(size = 16, color = "black"),
    legend.text = element_text(size=16, color = "black" ),
    axis.text.y  =element_text(size=20, color = "black"),
    axis.text.x =element_text(size=16, angle = 45, vjust =1,hjust =1, color = "black", margin = margin(t = 0, r = 10, b = 10, l = 0)),
    axis.title=element_text(size=20, color = "black"),
    plot.title = element_text(size=19,hjust=0.5, color = "black")
  )+
  labs(title="Translational protein",x ="", y = "LFCs")+
  scale_y_continuous(limits = c(-0.1, 0.5), breaks = seq(-0.2,0.6,0.2))
ggsave("CROP_R/4_Transcriptional_interactions/Translational_protein_enriched_in_3_DKOs.pdf", width = 4.2, height =4.6)
t.test(translation_M$`NF2-TP53`, translation_M$NF2, alternative = "greater")

##3. Enriched Cytoskeleton protein
cytoskeleton = c("MYL6", "ACTB", "TPM2", "ATG1", "PLEC", "CFL1", "TUBB", "TPM3", "NUDCD2", "ARPC5", "ARPC4", "WDR1")
#cytoskeleton = c("MYL6", "ACTB", "TPM2", "ATG1",  "CFL1",  "TPM3",  "ARPC5", "ARPC4", "WDR1")
cytoskeleton_M = comSyn_up[cytoskeleton, 2:7]
colnames(cytoskeleton_M) = c("NF2","PTEN", "TP53", "NF2-PTEN", "NF2-TP53", "PTEN-TP53")
d3 = melt(cytoskeleton_M)
d3$Group = c(rep("SKOs", nrow(d3)/2), rep("DKOs", nrow(d3)/2) )

vars = factor(d3$variable, levels = c("NF2", "PTEN", "TP53", "NF2-PTEN", "PTEN-TP53", "NF2-TP53"))
ggplot(d3, aes(x = vars, y = value, fill = Group))+
  geom_boxplot(width = 0.7)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.2), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line = element_line(colour = "black"),
    # legend.position = c(0.8,0.7),
    #legend.position = "None",
    legend.title = element_text(size = 16, color = "black"),
    legend.text = element_text(size=16, color = "black" ),
    axis.text.y  =element_text(size=20, color = "black"),
    axis.text.x =element_text(size=16, angle = 45, vjust =1,hjust =1, color = "black", margin = margin(t = 0, r = 10, b = 10, l = 0)),
    axis.title=element_text(size=20, color = "black"),
    plot.title = element_text(size=19,hjust=0.5, color = "black")
  )+
  labs(title="Cytoskeleton protein",x ="", y = "LFCs")+
  scale_y_continuous(limits = c(-0.1, 0.7), breaks = seq(-0.2,0.6,0.2))
ggsave("CROP_R/4_Transcriptional_interactions/Cytoskeleton_protein_enriched_in_3_DKOs.pdf", width = 4.2, height =4.6)
t.test(cytoskeleton_M $`NF2-TP53`, cytoskeleton_M$NF2,alternative = "greater" )
t.test(cytoskeleton_M $`NF2-PTEN`, cytoskeleton_M$NF2,alternative = "greater" )

##4. Enriched Acyltransferase
acy = c("SOAT2","ACAT2", "RABGGTB", "FDPS")
acy_M = comSyn_up[acy, 2:7]
colnames(acy_M) = c("NF2","PTEN", "TP53", "NF2-PTEN", "NF2-TP53", "PTEN-TP53")
d4 = melt(acy_M)
d4$Group = c(rep("SKOs", nrow(d4)/2), rep("DKOs", nrow(d4)/2) )

vars = factor(d4$variable, levels = c("NF2", "PTEN", "TP53", "NF2-PTEN", "PTEN-TP53", "NF2-TP53"))
ggplot(d4, aes(x = vars, y = value, fill = Group))+
  geom_boxplot(width = 0.7)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.2),  
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line = element_line(colour = "black"),
    # legend.position = c(0.8,0.7),
    #legend.position = "None",
    legend.title = element_text(size = 16, color = "black"),
    legend.text = element_text(size=16, color = "black" ),
    axis.text.y  =element_text(size=20, color = "black"),
    axis.text.x =element_text(size=16, angle = 45, vjust =1,hjust =1, color = "black", margin = margin(t = 0, r = 10, b = 10, l = 0)),
    axis.title=element_text(size=20, color = "black"),
    plot.title = element_text(size=19,hjust=0.5, color = "black")
  )+
  labs(title="Acyltransferase",x ="", y = "LFCs")+
  scale_y_continuous(limits = c(0, 0.61), breaks = seq(0.0, 0.6, 0.2))
ggsave("CROP_R/4_Transcriptional_interactions/Acyltransferase_enriched_in_3_DKOs.pdf", width = 4.2, height =4.6)

#---------------------------------------IV. Figures on common synergistic downregulated DEGs------------------------------------------------------
syn_dn_com_M = comSyn_dn[, 2:7]
colnames(syn_dn_com_M) = c("NF2","PTEN", "TP53", "NF2-PTEN", "NF2-TP53", "PTEN-TP53")
d5 = melt(syn_dn_com_M)
d5$Group = c(rep("SKOs", nrow(d5)/2), rep("DKOs", nrow(d5)/2) )

vars = factor(d5$variable, levels = c("NF2", "PTEN", "TP53", "NF2-PTEN", "PTEN-TP53", "NF2-TP53"))
ggplot(d5, aes(x = vars, y = -value, fill = Group))+
  geom_boxplot(width = 0.7)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.2),  
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line = element_line(colour = "black"),
    # legend.position = c(0.8,0.7),
    #legend.position = "None",
    legend.title = element_text(size = 16, color = "black"),
    legend.text = element_text(size=16, color = "black" ),
    axis.text.y  =element_text(size=20, color = "black"),
    axis.text.x =element_text(size=16, angle = 45, vjust =1,hjust =1, color = "black", margin = margin(t = 0, r = 10, b = 10, l = 0)),
    axis.title=element_text(size=20, color = "black"),
    plot.title = element_text(size=20,hjust=0.5, color = "black")
  )+
  labs(title="9 common synergies",x ="", y = "- LFCs")
ggsave("CROP_R/4_Transcriptional_interactions/9_common_synergistic_downregulated_DEGs_in_3_DKOs.pdf", width = 4.0, height =5)
t.test(splicing_M$`NF2-TP53`, splicing_M$NF2,  alternative = "greater")





















