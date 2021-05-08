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
#---------------------------------------------------------------------------------------------------------------------------------
tI = read.csv("CROP_R/4_Transcriptional_interactions/Groups_of_TI_categorization_5069_genes_based_on_TDs_3DKOs.csv", row.names = 1)
sname = rownames(tI)
#--------------------------------------I. Determine which TF targets do the 5069 genes belong to?-------------------------------------------------------------------------------------------
threshold = 0.01
#To determine NF2/YAP1 targets
chip_yap1 = read.delim("CROP_R/4_Transcriptional_interactions/YAP1_binding_ChipAtlas.txt", row.names = 1) #Downloaded from CHIP-Atlas
col_sub = c("SRX4213898.MCF_10A", "SRX4213899.MCF_10A", "SRX5287688.MCF.7", "SRX5287705.MCF.7", "SRX883576.MDA.MB.231", "SRX883577.MDA.MB.231")
chip_yap1_sub = chip_yap1[,col_sub]

chip_yap1_subM = chip_yap1_sub[match(sname, rownames(chip_yap1_sub)),]
rownames(chip_yap1_subM) = sname
chip_yap1_subM[is.na(chip_yap1_subM)] = 0

tgts_yap1 = cbind.data.frame(ensembl_id = tI$Ensembl_ID, lfc_NF2 = tI$lfc_NF2, chip_yap1_subM[sname, ])
tgts_yap1$score.MCF10A = rowSums(tgts_yap1[,c("SRX4213898.MCF_10A", "SRX4213899.MCF_10A")])
tgts_yap1$score.MCF7_MDA231 = rowSums(tgts_yap1[,c("SRX5287688.MCF.7", "SRX5287705.MCF.7", "SRX883576.MDA.MB.231", "SRX883577.MDA.MB.231")])

tgts_yap1$YAP1_targets = "no"
tgts_yap1$YAP1_targets[which(tgts_yap1$lfc_NF2 > threshold & tgts_yap1$score.MCF10A > 0 & tgts_yap1$score.MCF7_MDA231 > 0)] = "YAP1_activator"
sum(tgts_yap1$lfc_NF2 > threshold & tgts_yap1$score.MCF10A > 0 & tgts_yap1$score.MCF7_MDA231 > 0)
tgts_yap1$YAP1_targets[which(tgts_yap1$lfc_NF2 < -threshold & tgts_yap1$score.MCF10A > 0 & tgts_yap1$score.MCF7_MDA231 > 0)] = "YAP1_repressor"
sum(tgts_yap1$lfc_NF2 < -threshold & tgts_yap1$score.MCF10A > 0 & tgts_yap1$score.MCF7_MDA231 > 0)
which(rownames(tgts_yap1) == "CKS1B")
#---------------------------------------------------------------------------------------------------------------------------------
#To determine PTEN/FOXO1 targets 
chip_foxo1 = read.delim("CROP_R/4_Transcriptional_interactions/FOXO1_binding_ChipAtlas.txt", row.names = 1)  #Downloaded from CHIP-Atlas
chip_foxo1_sub = data.frame(SRX3230384.Hep_G2 = chip_foxo1$SRX3230384.Hep_G2, row.names = rownames(chip_foxo1))

chip_foxo1_subM = chip_foxo1_sub[match(sname, rownames(chip_foxo1_sub)),]
names(chip_foxo1_subM) = sname
chip_foxo1_subM[is.na(chip_foxo1_subM)] = 0

tgts_foxo1 = cbind.data.frame(ensembl_id = tI$Ensembl_ID, lfc_PTEN = tI$lfc_PTEN, SRX3230384.Hep_G2 = chip_foxo1_subM[sname])
tgts_foxo1$FOXO1_targets = "no"
tgts_foxo1$FOXO1_targets[which(tgts_foxo1$lfc_PTEN > threshold & tgts_foxo1$SRX3230384.Hep_G2 > 0)] = "FOXO1_activator"
sum(tgts_foxo1$lfc_PTEN > threshold & tgts_foxo1$SRX3230384.Hep_G2 > 0)
tgts_foxo1$FOXO1_targets[which(tgts_foxo1$lfc_PTEN < -threshold & tgts_foxo1$SRX3230384.Hep_G2  > 0)] = "FOXO1_repressor"
sum(tgts_foxo1$lfc_PTEN < -threshold & tgts_foxo1$SRX3230384.Hep_G2  > 0)
#---------------------------------------------------------------------------------------------------------------------------------
#To determine TP53 targets 
library(readxl)
activator_tp53 = data.frame(read_excel("CROP_R/4_Transcriptional_interactions/TP53_targets_Fischer_2017_PMC5511239.xlsx")) #Downloaded from PMC5511239
repressor_tp53 = data.frame(read_excel("CROP_R/4_Transcriptional_interactions/DREAM_targets_Fischer_2016_PMC4994865.xlsx")) #Downloaded from PMC5511239

tgts_tp53 = cbind.data.frame(ensembl_id = tI$Ensembl_ID, lfc_TP53 = tI$lfc_TP53, row.names = sname)
tgts_tp53$TP53_targets = "no"

tgts_tp53$TP53_targets[tgts_tp53$ensembl_id  %in% activator_tp53$ensembl.ID] = "TP53_activator"
#tgts_tp53$TP53_targets[sname %in% activator_tp53$Gene.Symbol] = "TP53_activator"
sum(tgts_tp53$ensembl_id  %in% activator_tp53$ensembl.ID)
#sum(sname %in% activator_tp53$Gene.Symbol)
tgts_tp53$TP53_targets[tgts_tp53$ensembl_id  %in% repressor_tp53$ensembl.ID] = "TP53_repressor"
#tgts_tp53$TP53_targets[sname %in% repressor_tp53$Gene.Symbol] = "TP53_repressor"
sum(tgts_tp53$ensembl_id  %in% repressor_tp53$ensembl.ID)
#sum(sname %in% repressor_tp53$Gene.Symbol)
#---------------------------------------------------------------------------------------------------------------------------------
tI_M = cbind.data.frame(tI, YAP1_targets = tgts_yap1[sname, "YAP1_targets"], FOXO1_targets = tgts_foxo1[sname, "FOXO1_targets"], TP53_targets = tgts_tp53[sname, "TP53_targets"])
which(rownames(tI_M) == "CKS1B")
write.csv(tI_M, "CROP_R/4_Transcriptional_interactions/TF_targets_and_Groups_of_TI_categorization_5069_genes_based_on_TDs_3DKOs.csv")


#--------------------------------------II. Determine which TF targets does gene-specific DEGs belong to?-------------------------------------------------------------------------------------------
#---1. 442 gene-specific synergistic upregulated DEGs---
tI_syn_up_spec = tI_M[which(tI_M$Group_Synergistic_Up %in% c("Synergistic_Up: NF2_specific", "Synergistic_Up: PTEN_specific", "Synergistic_Up: TP53_specific")),]
write.csv(tI_syn_up_spec, "CROP_R/4_Transcriptional_interactions/TF_targets_442_gene_specific_synergistic_upregulated_DEGs.csv")

#---------------------------------------------------------------------------------------------------------------------------------
#Total number of TP53_repressor_targets (424), FOXO1_activator_targets (156), YAP1_activator_targets (146) among 5069 genes
tgts_num = c(sum(tgts_tp53$ensembl_id %in% repressor_tp53$ensembl.ID), sum(tgts_foxo1$lfc_PTEN > threshold & tgts_foxo1$SRX3230384.Hep_G2 > 0),
                   sum(tgts_yap1$lfc_NF2 > threshold & tgts_yap1$score.MCF10A > 0 & tgts_yap1$score.MCF7_MDA231 > 0))
#---------------------------------------------------------------------------------------------------------------------------------
#Among TP53_specific synergistic genes, number of TP53_repressors, FOXO1_activators, YAP1_activators
num_TP53_specific = sum(tI_syn_up_spec$Group_Synergistic_Up == "Synergistic_Up: TP53_specific")
tgts_num_TP53_specific = c(sum(tI_syn_up_spec$Group_Synergistic_Up == "Synergistic_Up: TP53_specific" & tI_syn_up_spec$TP53_targets == "TP53_repressor"), 
                           sum(tI_syn_up_spec$Group_Synergistic_Up == "Synergistic_Up: TP53_specific" & tI_syn_up_spec$FOXO1_targets == "FOXO1_activator"), 
                           sum(tI_syn_up_spec$Group_Synergistic_Up == "Synergistic_Up: TP53_specific" & tI_syn_up_spec$YAP1_targets == "YAP1_activator")) 
#Among PTEN_specific synergistic genes, number of TP53_repressors, FOXO1_activators, YAP1_activators
num_PTEN_specific = sum(tI_syn_up_spec$Group_Synergistic_Up == "Synergistic_Up: PTEN_specific")
tgts_num_PTEN_specific = c(sum(tI_syn_up_spec$Group_Synergistic_Up == "Synergistic_Up: PTEN_specific" & tI_syn_up_spec$TP53_targets == "TP53_repressor"), 
                           sum(tI_syn_up_spec$Group_Synergistic_Up == "Synergistic_Up: PTEN_specific" & tI_syn_up_spec$FOXO1_targets == "FOXO1_activator"), 
                           sum(tI_syn_up_spec$Group_Synergistic_Up == "Synergistic_Up: PTEN_specific" & tI_syn_up_spec$YAP1_targets == "YAP1_activator")) 
#Among NF2_specific synergistic genes, number of TP53_repressors, FOXO1_activators, YAP1_activators
num_NF2_specific = sum(tI_syn_up_spec$Group_Synergistic_Up == "Synergistic_Up: NF2_specific")
tgts_num_NF2_specific = c(sum(tI_syn_up_spec$Group_Synergistic_Up == "Synergistic_Up: NF2_specific" & tI_syn_up_spec$TP53_targets == "TP53_repressor"), 
                          sum(tI_syn_up_spec$Group_Synergistic_Up == "Synergistic_Up: NF2_specific" & tI_syn_up_spec$FOXO1_targets == "FOXO1_activator"), 
                          sum(tI_syn_up_spec$Group_Synergistic_Up == "Synergistic_Up: NF2_specific" & tI_syn_up_spec$YAP1_targets == "YAP1_activator")) 

#A dataframe: Among each gene_specific synergistic genes, number of TP53_repressors, FOXO1_activators, YAP1_activators
df_up = cbind.data.frame(tgts_num_TP53_specific = tgts_num_TP53_specific, num_TP53_specific = num_TP53_specific, 
                         tgts_num_PTEN_specific = tgts_num_PTEN_specific, num_PTEN_specific = num_PTEN_specific, 
                         tgts_num_NF2_specific = tgts_num_NF2_specific, num_NF2_specific = num_NF2_specific, 
                         tgts_num = c(sum(tgts_tp53$ensembl_id %in% repressor_tp53$ensembl.ID), sum(tgts_foxo1$lfc_PTEN > threshold & tgts_foxo1$SRX3230384.Hep_G2 > 0),
                              sum(tgts_yap1$lfc_NF2 > threshold & tgts_yap1$score.MCF10A > 0 & tgts_yap1$score.MCF7_MDA231 > 0)))
write.csv(df_up, "CROP_R/4_Transcriptional_interactions/Number_TF_targets_in_gene_specific_synergistic_upregulated_DEGs.csv")

#Enrichment scores of TFs: Normalized to percentages of targets in all genes
tf_Tgts_up = cbind.data.frame(df_up$tgts_num_TP53_specific/df_up$num_TP53_specific/(df_up$tgts_num/5069),
                              df_up$tgts_num_PTEN_specific/df_up$num_PTEN_specific/(df_up$tgts_num/5069),
                              df_up$tgts_num_NF2_specific/df_up$num_NF2_specific/(df_up$tgts_num/5069),
                              df_up$tgts_num/5069/(df_up$tgts_num/5069))
rownames(tf_Tgts_up) = c("TP53_repressor", "FOXO1_activator","YAP1_activator")
colnames(tf_Tgts_up) = c("TP53 specific", "PTEN specific", "NF2 specific", "All genes")

tf_Tgts_up$TF_targets = rownames(tf_Tgts_up)
tf_enrich_score_up = melt(tf_Tgts_up)
tf_enrich_score_up$tf_targets = factor(tf_enrich_score_up$TF_targets, levels = c("TP53_repressor", "YAP1_activator", "FOXO1_activator"))
tf_enrich_score_up$group = factor(tf_enrich_score_up$variable, levels = c("All genes", "TP53 specific", "PTEN specific", "NF2 specific"))

ggplot(tf_enrich_score_up , aes(x = group, y = value, fill = tf_targets))+
  geom_bar(stat="identity", width = 0.7, position=position_dodge(),colour="black")+
theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        axis.line = element_line(colour = "black"),
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size=15 , color = "black"),
        axis.text.y  =element_text(size=15, color = "black"),
        axis.text.x =element_text(size=16, angle = 45, vjust = 1, hjust= 1, color = "black"),
        axis.title.x = element_text(size=16),
        axis.title.y  =element_text(size=17),
        plot.title = element_text(size=16, hjust = 0.5)
  )+
  scale_fill_manual(values = c("mediumorchid1","orange","#5AC8FA"),
                    name="TF targets",
                    labels=c( "TP53 repressor", "YAP1 activator", "FOXO1 activator"))+
  scale_y_continuous(limits = c(0,2.8), breaks = seq(0,3,0.5))+
  labs(title="Upregulated synergy",x ="", y = "Enrichment scores of TF targets")
ggsave("CROP_R/4_Transcriptional_interactions/Enrichment_scores_TF_targets_gene_specific_upregulated_synergy.pdf", width=5.8, height = 4)

#---------------------------------------------------------------------------------------------------------------------------------
#---2. 51 gene_specific synergistic upregulated DEGs---
tI_syn_dn_spec = tI_M[which(tI_M$Group_Synergistic_Dn %in% c("Synergistic_Dn: NF2_specific", "Synergistic_Dn: PTEN_specific", "Synergistic_Dn: TP53_specific")),]
write.csv(tI_syn_dn_spec, "CROP_R/4_Transcriptional_interactions/TF_targets_51_gene_specific_synergistic_downregulated_DEGs.csv")

#---------------------------------------------------------------------------------------------------------------------------------
#Total number of TP53_activators (151), FOXO1_repressors (73), YAP1_repressors (70) among 5069 genes
tgts_num = c(sum(tgts_tp53$ensembl_id %in% activator_tp53$ensembl.ID), sum(tgts_foxo1$lfc_PTEN < -threshold & tgts_foxo1$SRX3230384.Hep_G2 > 0),
             sum(tgts_yap1$lfc_NF2 < -threshold & tgts_yap1$score.MCF10A > 0 & tgts_yap1$score.MCF7_MDA231 > 0))
tgts_num

#---------------------------------------------------------------------------------------------------------------------------------
#Among TP53_specific synergistic genes, number of TP53_activators, FOXO1_repressors, YAP1_repressors
num_TP53_specific = sum(tI_syn_dn_spec$Group_Synergistic_Dn == "Synergistic_Dn: TP53_specific")
tgts_num_TP53_specific = c(sum(tI_syn_dn_spec$Group_Synergistic_Dn == "Synergistic_Dn: TP53_specific" & tI_syn_dn_spec$TP53_targets == "TP53_activator"), 
                           sum(tI_syn_dn_spec$Group_Synergistic_Dn == "Synergistic_Dn: TP53_specific" & tI_syn_dn_spec$FOXO1_targets == "FOXO1_repressor"), 
                           sum(tI_syn_dn_spec$Group_Synergistic_Dn == "Synergistic_Dn: TP53_specific" & tI_syn_dn_spec$YAP1_targets == "YAP1_repressor")) 
#Among PTEN_specific synergistic genes, number of TP53_activators, FOXO1_repressors, YAP1_repressors
num_PTEN_specific = sum(tI_syn_dn_spec$Group_Synergistic_Dn == "Synergistic_Dn: PTEN_specific")
tgts_num_PTEN_specific = c(sum(tI_syn_dn_spec$Group_Synergistic_Dn == "Synergistic_Dn: PTEN_specific" & tI_syn_dn_spec$TP53_targets == "TP53_activator"), 
                           sum(tI_syn_dn_spec$Group_Synergistic_Dn == "Synergistic_Dn: PTEN_specific" & tI_syn_dn_spec$FOXO1_targets == "FOXO1_repressor"), 
                           sum(tI_syn_dn_spec$Group_Synergistic_Dn == "Synergistic_Dn: PTEN_specific" & tI_syn_dn_spec$YAP1_targets == "YAP1_repressor")) 
#Among NF2_specific synergistic genes, number of TP53_activators, FOXO1_repressors, YAP1_repressors
num_NF2_specific = sum(tI_syn_dn_spec$Group_Synergistic_Dn == "Synergistic_Dn: NF2_specific")
tgts_num_NF2_specific = c(sum(tI_syn_dn_spec$Group_Synergistic_Dn == "Synergistic_Dn: NF2_specific" & tI_syn_dn_spec$TP53_targets == "TP53_activator"), 
                          sum(tI_syn_dn_spec$Group_Synergistic_Dn == "Synergistic_Dn: NF2_specific" & tI_syn_dn_spec$FOXO1_targets == "FOXO1_repressor"), 
                          sum(tI_syn_dn_spec$Group_Synergistic_Dn == "Synergistic_Dn: NF2_specific" & tI_syn_dn_spec$YAP1_targets == "YAP1_repressor")) 

#A dataframe: Among each gene_specific synergistic genes, number of TP53_activators, FOXO1_repressors, YAP1_repressors
df_dn = cbind.data.frame(tgts_num_TP53_specific = tgts_num_TP53_specific, num_TP53_specific = num_TP53_specific, 
                         tgts_num_PTEN_specific = tgts_num_PTEN_specific, num_PTEN_specific = num_PTEN_specific, 
                         tgts_num_NF2_specific = tgts_num_NF2_specific, num_NF2_specific = num_NF2_specific, 
                         tgts_num = c(sum(tgts_tp53$ensembl_id %in% activator_tp53$ensembl.ID), sum(tgts_foxo1$lfc_PTEN < -threshold & tgts_foxo1$SRX3230384.Hep_G2 > 0),
                                      sum(tgts_yap1$lfc_NF2 < -threshold & tgts_yap1$score.MCF10A > 0 & tgts_yap1$score.MCF7_MDA231 > 0)))
write.csv(df_dn, "CROP_R/4_Transcriptional_interactions/Number_TF_targets_in_gene_specific_synergistic_downregulated_DEGs.csv")

tf_Tgts_dn = cbind.data.frame(df_dn$tgts_num_TP53_specific/df_dn$num_TP53_specific/(df_dn$tgts_num/5069),
                        df_dn$tgts_num_PTEN_specific/df_dn$num_PTEN_specific/(df_dn$tgts_num/5069),
                        df_dn$tgts_num_NF2_specific/df_dn$num_NF2_specific/(df_dn$tgts_num/5069),
                        df_dn$tgts_num/5069/(df_dn$tgts_num/5069))
rownames(tf_Tgts_dn) = c("TP53_activator", "FOXO1_repressor","YAP1_repressor")
colnames(tf_Tgts_dn) = c("TP53 specific", "PTEN specific", "NF2 specific", "All genes")

tf_Tgts_dn$TF_targets = rownames(tf_Tgts_dn)
tf_enrich_score_dn = melt(tf_Tgts_dn)
tf_enrich_score_dn$tf_targets = factor(tf_enrich_score_dn$TF_targets, levels = c("TP53_activator", "YAP1_repressor", "FOXO1_repressor"))
tf_enrich_score_dn$group = factor(tf_enrich_score_dn$variable, levels = c("All genes", "TP53 specific", "PTEN specific", "NF2 specific"))

ggplot(tf_enrich_score_dn , aes(x = group, y = value, fill = tf_targets))+
  geom_bar(stat="identity", width = 0.7, position=position_dodge(),colour="black")+
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        axis.line = element_line(colour = "black"),
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size=15 , color = "black"),
        axis.text.y  =element_text(size=15, color = "black"),
        axis.text.x =element_text(size=16, angle = 45, vjust = 1, hjust= 1, color = "black"),
        axis.title.x = element_text(size=16),
        axis.title.y  =element_text(size=17),
        plot.title = element_text(size=16, hjust = 0.5)
  )+
  scale_fill_manual(values = c("mediumorchid1","orange","#5AC8FA"),
                    name="TF targets",
                    labels=c( "TP53 activator", "YAP1 repressor", "FOXO1 repressor"))+
  scale_y_continuous(limits = c(0,13.5), breaks = seq(0,14,2))+
  labs(title="Downregulated synergy",x ="", y = "Enrichment scores of\n TF targets")
ggsave("CROP_R/4_Transcriptional_interactions/Enrichment_scores_TF_targets_gene_specific_downregulated_synergy.pdf", width=5.8, height = 4)


#--------------------------------------III. Expression of TF targets in DKOs-------------------------------------------------------------------------------------------
#1. Expression of 22 enriched targets of TP53 repressor in TP53-specific upregulated synergistic genes
tp53_repressor_up = tI_syn_up_spec[tI_syn_up_spec$Group_Synergistic_Up == "Synergistic_Up: TP53_specific" & tI_syn_up_spec$TP53_targets == "TP53_repressor" ,2:7]
colnames(tp53_repressor_up) = c("NF2","PTEN", "TP53", "NF2-PTEN", "NF2-TP53", "PTEN-TP53")
lfc_tp53_repressor_up = melt(tp53_repressor_up)

lfc_tp53_repressor_up$group = factor(lfc_tp53_repressor_up$variable, levels = c("NF2", "PTEN", "TP53", "NF2-PTEN","PTEN-TP53", "NF2-TP53"))
ggplot(lfc_tp53_repressor_up, aes(x = group, y = value, color = group))+
  geom_boxplot(width = 0.5, size = 1)+
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        legend.position = "None",
        axis.line = element_line(colour = "black"),
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size=15 , color = "black"),
        axis.text.y  =element_text(size=15, color = "black"),
        axis.text.x =element_text(size=15, angle = 45, vjust = 1, hjust= 1, color = "black"),
        axis.title.x = element_text(size=16),
        axis.title.y  =element_text(size=17),
        plot.title = element_text(size=16, hjust = 0.5)
  )+
  labs(title="",x ="", y = "LFCs of TP53 repressors")
ggsave("CROP_R/4_Transcriptional_interactions/Expression_of_targets_TP53_repressors_in_TP53_specific_upregulated_synergy.pdf", width = 2.7, height =4)

#---------------------------------------------------------------------------------------------------------------------------------
#2. Expression of 13 enriched targets of FOXO1 activator in PTEN-specific upregulated synergistic genes
foxo1_activator_up = tI_syn_up_spec[tI_syn_up_spec$Group_Synergistic_Up == "Synergistic_Up: PTEN_specific" & tI_syn_up_spec$FOXO1_targets == "FOXO1_activator" ,2:7]
colnames(foxo1_activator_up) = c("NF2","PTEN", "TP53", "NF2-PTEN", "NF2-TP53", "PTEN-TP53")
lfc_foxo1_activator_up = melt(foxo1_activator_up)

lfc_foxo1_activator_up$group = factor(lfc_foxo1_activator_up$variable, levels = c("NF2", "PTEN", "TP53", "NF2-PTEN","PTEN-TP53", "NF2-TP53"))
ggplot(lfc_foxo1_activator_up, aes(x = group, y = value, color = group))+
  geom_boxplot(width = 0.5, size = 1)+
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        legend.position = "None",
        axis.line = element_line(colour = "black"),
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size=15 , color = "black"),
        axis.text.y  =element_text(size=15, color = "black"),
        axis.text.x =element_text(size=15, angle = 45, vjust = 1, hjust= 1, color = "black"),
        axis.title.x = element_text(size=16),
        axis.title.y  =element_text(size=17),
        plot.title = element_text(size=16, hjust = 0.5)
  )+
  labs(title="",x ="", y = "LFCs of FOXO1 activators")
ggsave("CROP_R/4_Transcriptional_interactions/Expression_of_targets_FOXO1_activators_in_PTEN_specific_upregulated_synergy.pdf", width = 2.7, height =4)

#---------------------------------------------------------------------------------------------------------------------------------
#3. Expression of 11 enriched targets of YAP1 activator in NF2-specific upregulated synergistic genes
yap1_activator_up = tI_syn_up_spec[tI_syn_up_spec$Group_Synergistic_Up == "Synergistic_Up: NF2_specific" & tI_syn_up_spec$YAP1_targets == "YAP1_activator" ,2:7]

colnames(yap1_activator_up) = c("NF2","PTEN", "TP53", "NF2-PTEN", "NF2-TP53", "PTEN-TP53")
lfc_yap1_activator_up = melt(yap1_activator_up)

lfc_yap1_activator_up$group = factor(lfc_yap1_activator_up$variable, levels = c("NF2", "PTEN", "TP53", "NF2-PTEN","PTEN-TP53", "NF2-TP53"))
ggplot(lfc_yap1_activator_up, aes(x = group, y = value, color = group))+
  geom_boxplot(width = 0.5, size = 1)+
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        legend.position = "None",
        axis.line = element_line(colour = "black"),
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size=15 , color = "black"),
        axis.text.y  =element_text(size=15, color = "black"),
        axis.text.x =element_text(size=15, angle = 45, vjust = 1, hjust= 1, color = "black"),
        axis.title.x = element_text(size=16),
        axis.title.y  =element_text(size=17),
        plot.title = element_text(size=16, hjust = 0.5)
  )+
  labs(title="",x ="", y = "LFCs of YAP1 activators")
ggsave("CROP_R/4_Transcriptional_interactions/Expression_of_targets_YAP1_activators_in_NF2_specific_upregulated_synergy.pdf", width = 2.7, height =4)

#---------------------------------------------------------------------------------------------------------------------------------
#4. Expression CDK4 in TP53-specific upregulated synergistic genes
tp53_CDK4 = tI_syn_up_spec["CDK4",2:7]
lfc_tp53_CDK4$group = factor(lfc_tp53_CDK4$variable, levels = c("NF2", "PTEN", "TP53","NF2-TP53", "PTEN-TP53","NF2-PTEN"))
ggplot(lfc_tp53_CDK4, aes(x = group, y = value, fill = group))+
  geom_bar(stat="identity", width = 0.7, position=position_dodge(),colour="black")+
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        legend.position = "None",
        axis.line = element_line(colour = "black"),
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size=15 , color = "black"),
        axis.text.y  =element_text(size=15, color = "black"),
        axis.text.x =element_text(size=15, angle = 45, vjust = 1, hjust= 1, color = "black"),
        axis.title.x = element_text(size=16),
        axis.title.y  =element_text(size=17),
        plot.title = element_text(size=16, hjust = 0.5)
  )+
  labs(title="",x ="", y = "LFCs of CDK4")+
  scale_y_continuous(limits = c(0, 0.48), breaks = seq(0, 0.5,0.1))
ggsave("CROP_R/4_Transcriptional_interactions/CDK4_synergy_TP53_specific_upregulated_synergy.pdf", width = 2.7, height =4)

#---------------------------------------------------------------------------------------------------------------------------------
#5. Expression of 8 enriched targets of TP53 activators in TP53-specific downregulated synergistic genes
tp53_activator_dn = tI_syn_dn_spec[tI_syn_dn_spec$Group_Synergistic_Dn == "Synergistic_Dn: TP53_specific" & tI_syn_dn_spec$TP53_targets == "TP53_activator" ,2:7]
colnames(tp53_activator_dn) = c("NF2","PTEN", "TP53", "NF2-PTEN", "NF2-TP53", "PTEN-TP53")
lfc_tp53_activator_dn = melt(tp53_activator_dn)

lfc_tp53_activator_dn$group = factor(lfc_tp53_activator_dn$variable, levels = c("NF2", "PTEN", "TP53", "NF2-PTEN","PTEN-TP53", "NF2-TP53"))
ggplot(lfc_tp53_activator_dn, aes(x = group, y = -value, color = group))+
  geom_boxplot(width = 0.5, size = 1)+
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        legend.position = "None",
        axis.line = element_line(colour = "black"),
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size=15 , color = "black"),
        axis.text.y  =element_text(size=15, color = "black"),
        axis.text.x =element_text(size=15, angle = 45, vjust = 1, hjust= 1, color = "black"),
        axis.title.x = element_text(size=16),
        axis.title.y  =element_text(size=17),
        plot.title = element_text(size=16, hjust = 0.5)
  )+
  scale_y_continuous(limits = c(-0.3, 0.85), breaks = seq(-0.2, 0.8, 0.2))+
  labs(title="",x ="", y = "- LFCs of TP53 activators")
  
ggsave("CROP_R/4_Transcriptional_interactions/Expression_of_targets_TP53_activators_in_TP53_specific_downregulated_synergy.pdf", width = 2.7, height =4)

#---------------------------------------------------------------------------------------------------------------------------------
#6. Expression of 3 enriched targets of FOXO1 repressors in PTEN-specific downregulated synergistic genes
foxo1_repressor_dn =  tI_syn_dn_spec[ tI_syn_dn_spec$Group_Synergistic_Dn == "Synergistic_Dn: PTEN_specific" & tI_syn_dn_spec$FOXO1_targets == "FOXO1_repressor" ,2:7]
colnames(foxo1_repressor_dn) = c("NF2","PTEN", "TP53", "NF2-PTEN", "NF2-TP53", "PTEN-TP53")
lfc_foxo1_repressor_dn = melt(foxo1_repressor_dn)

lfc_foxo1_repressor_dn$group = factor(lfc_foxo1_repressor_dn$variable, levels = c("NF2", "PTEN", "TP53", "NF2-PTEN","PTEN-TP53", "NF2-TP53"))
ggplot(lfc_foxo1_repressor_dn, aes(x = group, y = -value, color = group))+
  geom_boxplot(width = 0.5, size = 1)+
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        legend.position = "None",
        axis.line = element_line(colour = "black"),
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size=15 , color = "black"),
        axis.text.y  =element_text(size=14, color = "black"),
        axis.text.x =element_text(size=15, angle = 45, vjust = 1, hjust= 1, color = "black"),
        axis.title.x = element_text(size=16),
        axis.title.y  =element_text(size=17),
        plot.title = element_text(size=16, hjust = 0.5)
  )+
  scale_y_continuous(limits = c(-0.7, 0.4), breaks = seq(-0.8, 0.4, 0.2))+
  labs(title="",x ="", y = "- LFCs of FOXO1 repressors")

ggsave("CROP_R/4_Transcriptional_interactions/Expression_of_targets_FOXO1_repressor_in_PTEN_specific_downregulated_synergy.pdf", width = 2.7, height =4)

#---------------------------------------------------------------------------------------------------------------------------------
#7. Expression of 3 enriched targets of FOXO1 repressors in PTEN-specific downregulated synergistic genes
yap1_repressor_dn =  tI_syn_dn_spec[ tI_syn_dn_spec$Group_Synergistic_Dn == "Synergistic_Dn: NF2_specific" & tI_syn_dn_spec$YAP1_targets == "YAP1_repressor" ,2:7]
colnames(yap1_repressor_dn) = c("NF2","PTEN", "TP53", "NF2-PTEN", "NF2-TP53", "PTEN-TP53")
lfc_yap1_repressor_dn = melt(yap1_repressor_dn)

lfc_yap1_repressor_dn$group = factor(lfc_yap1_repressor_dn$variable, levels = c("NF2", "PTEN", "TP53", "NF2-PTEN","PTEN-TP53", "NF2-TP53"))
ggplot(lfc_yap1_repressor_dn, aes(x = group, y = -value, fill = group))+
  geom_bar(stat="identity", width = 0.7, position=position_dodge(), color = "black")+
  #geom_boxplot(width = 0.5, size = 1)+
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        legend.position = "None",
        axis.line = element_line(colour = "black"),
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size=15 , color = "black"),
        axis.text.y  =element_text(size=15, color = "black"),
        axis.text.x =element_text(size=15, angle = 45, vjust = 1, hjust= 1, color = "black"),
        axis.title.x = element_text(size=16),
        axis.title.y  =element_text(size=17),
        plot.title = element_text(size=16, hjust = 0.5)
  )+
  #scale_y_continuous(limits = c(-0.85, 0.3), breaks = seq(-0.8, 0.2, 0.2))+
  labs(title="",x ="", y = "- LFCs of YAP1 repressors")

ggsave("CROP_R/4_Transcriptional_interactions/Expression_of_targets_YAP1_repressor_in_NF2_specific_downregulated_synergy.pdf", width = 2.7, height =4)


