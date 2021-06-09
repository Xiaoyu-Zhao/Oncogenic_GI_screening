#Xiaoyu Zhao 
library(ggplot2)
library(data.table) ## for read.table("xxx.xls",sep = "\t",header = T) 
library(pheatmap) 
library(Cairo) #https://stackoverflow.com/questions/27690729/greek-letters-symbols-and-line-breaks-inside-a-ggplot-legend-label
library(VennDiagram)
library(ggpubr)
library(tidyr)
library(UpSetR)

setwd("./5_Combinatorial_CROP_seq/")
#-------------------------------- I. Common DEGs induced by sgNF2/sgPTEN/sgSMAD4 in S2, S3 and S4--------------------------------------
#------Number of DEGs of sgNF2 in S2 and S4--------
lfc_S2_nf2 =read.csv("Python_code/sequencing/outs_S2/results/res_Deseq2/res_NF2-CTRL1.csv", row.names = 1)
lfc_S4_nf2 =read.csv("Python_code/sequencing/outs_S4/results/res_Deseq2/res_NF2-CTRL1.csv", row.names = 1)
Nm_S2_nf2 = lfc_S2_nf2[which(as.numeric(lfc_S2_nf2$padj) < 0.05),]
Nm_S4_nf2 = lfc_S4_nf2[which(as.numeric(lfc_S4_nf2$padj) < 0.05),]

Nms_nf2 = unique(unlist(c(rownames(Nm_S2_nf2), rownames(Nm_S4_nf2))))
totalM_nf2 = data.frame(matrix(0, nrow = length(Nms_nf2), ncol = 2))
rownames(totalM_nf2) = Nms_nf2
colnames(totalM_nf2) = c("S2",  "S4")
totalM_nf2[match(rownames(Nm_S2_nf2), rownames(totalM_nf2)), "S2"] = 1
#totalM[match(rownames(Nm_S3), rownames(totalM)), "S3"] = 1
totalM_nf2[match(rownames(Nm_S4_nf2), rownames(totalM_nf2)), "S4"] = 1

#------Number of DEGs of sgPTEN in S2 and S4--------
lfc_S2_pten =read.csv("Python_code/sequencing/outs_S2/results/res_Deseq2/res_PTEN-CTRL1.csv", row.names = 1)
lfc_S4_pten =read.csv("Python_code/sequencing/outs_S4/results/res_Deseq2/res_PTEN-CTRL1.csv", row.names = 1)
Nm_S2_pten = lfc_S2_pten[which(as.numeric(lfc_S2_pten$padj) < 0.05),]
Nm_S4_pten = lfc_S4_pten[which(as.numeric(lfc_S4_pten$padj) < 0.05),]

Nms_pten = unique(unlist(c(rownames(Nm_S2_pten), rownames(Nm_S4_pten))))
totalM_pten = data.frame(matrix(0, nrow = length(Nms_pten), ncol = 2))
rownames(totalM_pten) = Nms_pten
colnames(totalM_pten) = c("S2",  "S4")
totalM_pten[match(rownames(Nm_S2_pten), rownames(totalM_pten)), "S2"] = 1
#totalM[match(rownames(Nm_S3), rownames(totalM)), "S3"] = 1
totalM_pten[match(rownames(Nm_S4_pten), rownames(totalM_pten)), "S4"] = 1

#------Number of DEGs of sgSMAAD4 in S2 and S3--------
lfc_S2_smad4 =read.csv("Python_code/sequencing/outs_S2/results/res_Deseq2/res_SMAD4-CTRL1.csv", row.names = 1)
lfc_S3_smad4 =read.csv("Python_code/sequencing/outs_S3/results/res_Deseq2/res_SMAD4-CTRL1.csv", row.names = 1)
Nm_S2_smad4 = lfc_S2_smad4[which(as.numeric(lfc_S2_smad4$padj) < 0.05),]
Nm_S3_smad4 = lfc_S3_smad4[which(as.numeric(lfc_S3_smad4$padj) < 0.05),]

Nms_smad4 = unique(unlist(c(rownames(Nm_S2_smad4), rownames(Nm_S3_smad4))))
totalM_smad4 = data.frame(matrix(0, nrow = length(Nms_smad4), ncol = 2))
rownames(totalM_smad4) = Nms_smad4
colnames(totalM_smad4) = c("S2",  "S3")
totalM_smad4[match(rownames(Nm_S2_smad4), rownames(totalM_smad4)), "S2"] = 1
#totalM[match(rownames(Nm_S3), rownames(totalM)), "S3"] = 1
totalM_smad4[match(rownames(Nm_S3_smad4), rownames(totalM_smad4)), "S3"] = 1


##Function to find overlapping DEGs
overLap<- function(cond) {
  ppl <- totalM
  
  for (i in 1:length(cond)) {
    ppl <- subset(ppl, ppl[cond[i]] == T)
  }
  nrow(ppl)
}

a = c("S2", "S4")
b = c("S2", "S3")

#Venn diagram of sgNF2 in S2 and S4
pdf("R_code/3_Diffexpr_by_CRISPR/Comparison_SKO_induced_diffexpr_in_S2_S3_S4/Comparison_DEGs_sgNF2_S24.pdf", height = 3.5, width = 3.5)
totalM = totalM_nf2
grid.newpage()
overrideTriple=T  
##https://stackoverflow.com/questions/11727068/scaling-triple-venn-diagram-in-r-with-venndiagram-package
draw.pairwise.venn(overLap(a[1]), overLap(a[2]), overLap(a),
                 #category = c("PTEN-/-", "PIK3CA", "HAMyc"), 
                 lty = "blank", alpha = 0.5,
                 fill = c("#5AC8FA","mediumorchid"),
                 scaled=TRUE,fontfamily ="Helvetica", cat.fontfamily="Helvetica",
                 cex = 2, cat.cex = 2)
#To change the main font and category font
#https://stackoverflow.com/questions/56583545/how-to-change-venn-diagram-font
dev.off()

#Venn diagram of sgPTEN in S2 and S4
pdf("R_code/3_Diffexpr_by_CRISPR/Comparison_SKO_induced_diffexpr_in_S2_S3_S4/Comparison_DEGs_sgPTEN_S24.pdf", height = 3.5, width = 3.5)
totalM = totalM_pten
grid.newpage()
overrideTriple=T  
##https://stackoverflow.com/questions/11727068/scaling-triple-venn-diagram-in-r-with-venndiagram-package
draw.pairwise.venn(overLap(a[1]), overLap(a[2]), overLap(a),
                   #category = c("PTEN-/-", "PIK3CA", "HAMyc"), 
                   lty = "blank", alpha = 0.5,
                   fill = c("#5AC8FA","mediumorchid"),
                   scaled=TRUE,fontfamily ="Helvetica", cat.fontfamily="Helvetica",
                   cex = 2, cat.cex = 2)
#To change the main font and category font
#https://stackoverflow.com/questions/56583545/how-to-change-venn-diagram-font
dev.off()

#Venn diagram of sgSMAD4 in S2 and S3
pdf("R_code/3_Diffexpr_by_CRISPR/Comparison_SKO_induced_diffexpr_in_S2_S3_S4/Comparison_DEGs_sgSMAD4_S23.pdf", height = 3.5, width = 3.5)
totalM = totalM_smad4
grid.newpage()
overrideTriple=T  
##https://stackoverflow.com/questions/11727068/scaling-triple-venn-diagram-in-r-with-venndiagram-package
draw.pairwise.venn(overLap(b[1]), overLap(b[2]), overLap(b),
                   #category = c("PTEN-/-", "PIK3CA", "HAMyc"), 
                   lty = "blank", alpha = 0.5,
                   fill = c("#5AC8FA","mediumorchid"),
                   scaled=TRUE,fontfamily ="Helvetica", cat.fontfamily="Helvetica",
                   cex = 2, cat.cex = 2)
#To change the main font and category font
#https://stackoverflow.com/questions/56583545/how-to-change-venn-diagram-font
dev.off()

#--------------------- II. Enrichment of E2F_targets and G2M checkpoints of NF2 in S2, S3 and S4----------------------------
source("R_code/2_DE_conditions_20200908/gseaPlot_no_legend_20200908.R")
gs0_S2<- data.frame(read.table(file = "R_code/3_Diffexpr_by_CRISPR/S2_GSEA/NF2-CTRL1_vs_Control/NF2_hallmarks.Gsea.1603489868012/ranked_gene_list_NF2-CTRL1_versus_control_1603489868012.tsv",sep = "\t",header = T))
gs0_S3<- data.frame(read.table(file = "R_code/3_Diffexpr_by_CRISPR/S3_GSEA/NF2-CTRL1_vs_Control/NF2_hallmarks.Gsea.1603490524994/ranked_gene_list_NF2-CTRL1_versus_control_1603490524994.tsv",sep = "\t",header = T))
gs0_S4<- data.frame(read.table(file = "R_code/3_Diffexpr_by_CRISPR/S4_GSEA/NF2-CTRL1_vs_Control/NF2_hallmarks.Gsea.1603490996038/ranked_gene_list_NF2-CTRL1_versus_control_1603490996038.tsv",sep = "\t",header = T))

#E2F targets
gs1_S2 <- data.frame(read.table(file = "R_code/3_Diffexpr_by_CRISPR/S2_GSEA/NF2-CTRL1_vs_Control/NF2_hallmarks.Gsea.1603489868012/HALLMARK_E2F_TARGETS.tsv",sep = "\t",header = T))
gs1_S3 <- data.frame(read.table(file = "R_code/3_Diffexpr_by_CRISPR/S3_GSEA/NF2-CTRL1_vs_Control/NF2_hallmarks.Gsea.1603490524994/HALLMARK_E2F_TARGETS.tsv",sep = "\t",header = T))
gs1_S4 <- data.frame(read.table(file = "R_code/3_Diffexpr_by_CRISPR/S4_GSEA/NF2-CTRL1_vs_Control/NF2_hallmarks.Gsea.1603490996038/HALLMARK_E2F_TARGETS.tsv",sep = "\t",header = T))

#G2M checkpoint
gs2_S2 <- data.frame(read.table(file = "R_code/3_Diffexpr_by_CRISPR/S2_GSEA/NF2-CTRL1_vs_Control/NF2_hallmarks.Gsea.1603489868012/HALLMARK_G2M_CHECKPOINT.tsv",sep = "\t",header = T))
gs2_S3 <- data.frame(read.table(file = "R_code/3_Diffexpr_by_CRISPR/S3_GSEA/NF2-CTRL1_vs_Control/NF2_hallmarks.Gsea.1603490524994/HALLMARK_G2M_CHECKPOINT.tsv",sep = "\t",header = T))
gs2_S4 <- data.frame(read.table(file = "R_code/3_Diffexpr_by_CRISPR/S4_GSEA/NF2-CTRL1_vs_Control/NF2_hallmarks.Gsea.1603490996038/HALLMARK_G2M_CHECKPOINT.tsv",sep = "\t",header = T))

#----------------------S2--------------------------
#E2F targets and G2M checkpoints of NF2-CTRL1 in S2
gs0_S2_1 = gs0_S2
gs0_S2_2 = gs0_S2

gs0_S2_1$position = 0
gs0_S2_2$position = 0

gs0_S2_1[match(gs1_S2$SYMBOL, gs0_S2_1$NAME), "position"] = 1
gs0_S2_2[match(gs2_S2$SYMBOL, gs0_S2_2$NAME), "position"] = 1

gs0_S2_1$Description = NULL
gs0_S2_2$Description = NULL

gs0_S2_1[match(gs1_S2$SYMBOL, gs0_S2_1$NAME), "Description"] = "E2F targets"
gs0_S2_2[match(gs2_S2$SYMBOL, gs0_S2_2$NAME), "Description"] = "G2M checkpoint"

gs0_S2_1$ES = 0
gs0_S2_2$ES = 0

gs0_S2_1[match(gs1_S2$SYMBOL, gs0_S2_1$NAME), "ES"] = gs1_S2$RUNNING.ES
gs0_S2_2[match(gs2_S2$SYMBOL, gs0_S2_2$NAME), "ES"] = gs2_S2$RUNNING.ES

gsdataX = rbind.data.frame(gs0_S2_1, gs0_S2_2)

gsdataY_nf2_S2 = data.frame(Name = gsdataX$NAME,
                     x = c(as.numeric(rownames(gs0_S2_1)),as.numeric(rownames(gs0_S2_2))),
                     y= as.numeric(gsdataX$SCORE),
                     Description = gsdataX$Description,
                     runningScore = gsdataX$ES,
                     position = gsdataX$position,
                     ymin = rep(0, nrow(gsdataX)),
                     ymax = rep(0, nrow(gsdataX)))

pos <- gsdataY_nf2_S2$position == 1
h <- diff(range(gsdataY_nf2_S2$runningScore))/20
gsdataY_nf2_S2$ymin[pos] <- -h
gsdataY_nf2_S2$ymax[pos] <- h

#source("R_code/2_DE_conditions_20200908/gseaPlot_20200908.R")
gseaPlot(gsdataY_nf2_S2, geneSetID = 1:4, title = "", color=c("red", "blue", "green", "purple"), dotsize = 4, base_size =40,
         rel_heights=c(1.5, 0.5, 0.5), subplots = 1:2, pvalue_table = FALSE, ES_geom="dot",
         y_min = -0.6, y_max =0.1, legend.position = c(.68, .94))
ggsave("R_code/3_Diffexpr_by_CRISPR/Comparison_SKO_induced_diffexpr_in_S2_S3_S4/Enrichment_E2F_G2M_sgNF2_S2.pdf", width = 8.2, height = 6.8)

#----------------------S3--------------------------
#E2F targets and G2M checkpoints of NF2-CTRL1 in S3
gs0_S3_1 = gs0_S3
gs0_S3_2 = gs0_S3

gs0_S3_1$position = 0
gs0_S3_2$position = 0

gs0_S3_1[match(gs1_S3$SYMBOL, gs0_S3_1$NAME), "position"] = 1
gs0_S3_2[match(gs2_S3$SYMBOL, gs0_S3_2$NAME), "position"] = 1

gs0_S3_1$Description = NULL
gs0_S3_2$Description = NULL

gs0_S3_1[match(gs1_S3$SYMBOL, gs0_S3_1$NAME), "Description"] = "E2F targets"
gs0_S3_2[match(gs2_S3$SYMBOL, gs0_S3_2$NAME), "Description"] = "G2M checkpoint"

gs0_S3_1$ES = 0
gs0_S3_2$ES = 0

gs0_S3_1[match(gs1_S3$SYMBOL, gs0_S3_1$NAME), "ES"] = gs1_S3$RUNNING.ES
gs0_S3_2[match(gs2_S3$SYMBOL, gs0_S3_2$NAME), "ES"] = gs2_S3$RUNNING.ES

gsdataX = rbind.data.frame(gs0_S3_1, gs0_S3_2)

gsdataY_nf2_S3 = data.frame(Name = gsdataX$NAME,
                     x = c(as.numeric(rownames(gs0_S3_1)),as.numeric(rownames(gs0_S3_2))),
                     y= as.numeric(gsdataX$SCORE),
                     Description = gsdataX$Description,
                     runningScore = gsdataX$ES,
                     position = gsdataX$position,
                     ymin = rep(0, nrow(gsdataX)),
                     ymax = rep(0, nrow(gsdataX)))

pos <- gsdataY_nf2_S3$position == 1
h <- diff(range(gsdataY_nf2_S3$runningScore))/20
gsdataY_nf2_S3$ymin[pos] <- -h
gsdataY_nf2_S3$ymax[pos] <- h

gseaPlot(gsdataY_nf2_S3, geneSetID = 1:4, title = "", color=c("red", "blue", "green", "purple"), dotsize = 4, base_size =40,
         rel_heights=c(1.5, 0.5, 0.5), subplots = 1:2, pvalue_table = FALSE, ES_geom="dot",
         y_min = -0.1, y_max = 0.6, legend.position = c(.71, .94))
ggsave("R_code/3_Diffexpr_by_CRISPR/Comparison_SKO_induced_diffexpr_in_S2_S3_S4/Enrichment_E2F_G2M_sgNF2_S3.pdf", width = 8.2, height = 6.8)

#----------------------S4--------------------------
#E2F targets and G2M checkpoints of NF2-CTRL1 in S4
gs0_S4_1 = gs0_S4
gs0_S4_2 = gs0_S4

gs0_S4_1$position = 0
gs0_S4_2$position = 0

gs0_S4_1[match(gs1_S4$SYMBOL, gs0_S4_1$NAME), "position"] = 1
gs0_S4_2[match(gs2_S4$SYMBOL, gs0_S4_2$NAME), "position"] = 1

gs0_S4_1$Description = NULL
gs0_S4_2$Description = NULL

gs0_S4_1[match(gs1_S4$SYMBOL, gs0_S4_1$NAME), "Description"] = "E2F targets"
gs0_S4_2[match(gs2_S4$SYMBOL, gs0_S4_2$NAME), "Description"] = "G2M checkpoint"

gs0_S4_1$ES = 0
gs0_S4_2$ES = 0

gs0_S4_1[match(gs1_S4$SYMBOL, gs0_S4_1$NAME), "ES"] = gs1_S4$RUNNING.ES
gs0_S4_2[match(gs2_S4$SYMBOL, gs0_S4_2$NAME), "ES"] = gs2_S4$RUNNING.ES

gsdataX = rbind.data.frame(gs0_S4_1, gs0_S4_2)

gsdataY_nf2_S4 = data.frame(Name = gsdataX$NAME,
                     x = c(as.numeric(rownames(gs0_S4_1)),as.numeric(rownames(gs0_S4_2))),
                     y= as.numeric(gsdataX$SCORE),
                     Description = gsdataX$Description,
                     runningScore = gsdataX$ES,
                     position = gsdataX$position,
                     ymin = rep(0, nrow(gsdataX)),
                     ymax = rep(0, nrow(gsdataX)))

pos <- gsdataY_nf2_S4$position == 1
h <- diff(range(gsdataY_nf2_S4$runningScore))/20
gsdataY_nf2_S4$ymin[pos] <- -h
gsdataY_nf2_S4$ymax[pos] <- h

gseaPlot(gsdataY_nf2_S4, geneSetID = 1:4, title = "", color=c("red", "blue", "green", "purple"), dotsize = 4, base_size =40,
         rel_heights=c(1.5, 0.5, 0.5), subplots = 1:2, pvalue_table = FALSE, ES_geom="dot",
         y_min = -0.1, y_max = 0.5, legend.position = c(.71, .94))
ggsave("R_code/3_Diffexpr_by_CRISPR/Comparison_SKO_induced_diffexpr_in_S2_S3_S4/Enrichment_E2F_G2M_sgNF2_S4.pdf", width = 8.2, height = 6.8)


#--------------------- II. Enrichment of E2F_targets and G2M checkpoints of PTEN in S2, S3 and S4----------------------------
gs0_S2<- data.frame(read.table(file = "R_code/3_Diffexpr_by_CRISPR/S2_GSEA/PTEN-CTRL1_vs_Control/PTEN_hallmarks.Gsea.1603739247127/ranked_gene_list_PTEN-CTRL1_versus_control_1603739247127.tsv",sep = "\t",header = T))
gs0_S4<- data.frame(read.table(file = "R_code/3_Diffexpr_by_CRISPR/S4_GSEA/PTEN-CTRL1_vs_Control/PTEN_hallmarks.Gsea.1603491024205/ranked_gene_list_PTEN-CTRL1_versus_control_1603491024205.tsv",sep = "\t",header = T))

#E2F targets
gs1_S2 <- data.frame(read.table(file = "R_code/3_Diffexpr_by_CRISPR/S2_GSEA/PTEN-CTRL1_vs_Control/PTEN_hallmarks.Gsea.1603739247127/HALLMARK_E2F_TARGETS.tsv",sep = "\t",header = T))
gs1_S4 <- data.frame(read.table(file = "R_code/3_Diffexpr_by_CRISPR/S4_GSEA/PTEN-CTRL1_vs_Control/PTEN_hallmarks.Gsea.1603491024205/HALLMARK_E2F_TARGETS.tsv",sep = "\t",header = T))

#G2M checkpoint
gs2_S2 <- data.frame(read.table(file = "R_code/3_Diffexpr_by_CRISPR/S2_GSEA/PTEN-CTRL1_vs_Control/PTEN_hallmarks.Gsea.1603739247127/HALLMARK_G2M_CHECKPOINT.tsv",sep = "\t",header = T))
gs2_S4 <- data.frame(read.table(file = "R_code/3_Diffexpr_by_CRISPR/S4_GSEA/PTEN-CTRL1_vs_Control/PTEN_hallmarks.Gsea.1603491024205/HALLMARK_G2M_CHECKPOINT.tsv",sep = "\t",header = T))

#----------------------S2--------------------------
#E2F targets and G2M checkpoints of PTEN-CTRL1 in S2
gs0_S2_1 = gs0_S2
gs0_S2_2 = gs0_S2

gs0_S2_1$position = 0
gs0_S2_2$position = 0

gs0_S2_1[match(gs1_S2$SYMBOL, gs0_S2_1$NAME), "position"] = 1
gs0_S2_2[match(gs2_S2$SYMBOL, gs0_S2_2$NAME), "position"] = 1

gs0_S2_1$Description = NULL
gs0_S2_2$Description = NULL

gs0_S2_1[match(gs1_S2$SYMBOL, gs0_S2_1$NAME), "Description"] = "E2F targets"
gs0_S2_2[match(gs2_S2$SYMBOL, gs0_S2_2$NAME), "Description"] = "G2M checkpoint"

gs0_S2_1$ES = 0
gs0_S2_2$ES = 0

gs0_S2_1[match(gs1_S2$SYMBOL, gs0_S2_1$NAME), "ES"] = gs1_S2$RUNNING.ES
gs0_S2_2[match(gs2_S2$SYMBOL, gs0_S2_2$NAME), "ES"] = gs2_S2$RUNNING.ES

gsdataX = rbind.data.frame(gs0_S2_1, gs0_S2_2)

gsdataY_pten_S2 = data.frame(Name = gsdataX$NAME,
                     x = c(as.numeric(rownames(gs0_S2_1)),as.numeric(rownames(gs0_S2_2))),
                     y= as.numeric(gsdataX$SCORE),
                     Description = gsdataX$Description,
                     runningScore = gsdataX$ES,
                     position = gsdataX$position,
                     ymin = rep(0, nrow(gsdataX)),
                     ymax = rep(0, nrow(gsdataX)))

pos <- gsdataY_pten_S2$position == 1
h <- diff(range(gsdataY_pten_S2$runningScore))/20
gsdataY_pten_S2$ymin[pos] <- -h
gsdataY_pten_S2$ymax[pos] <- h

gseaPlot(gsdataY_pten_S2, geneSetID = 1:4, title = "", color=c("red", "blue", "green", "purple"), dotsize = 4, base_size =40,
         rel_heights=c(1.5, 0.5, 0.5), subplots = 1:2, pvalue_table = FALSE, ES_geom="dot",
         y_min = -0.2, y_max =0.6, legend.position = c(.68, .94))
ggsave("R_code/3_Diffexpr_by_CRISPR/Comparison_SKO_induced_diffexpr_in_S2_S3_S4/Enrichment_E2F_G2M_sgPTEN_S2.pdf", width = 8.2, height = 6.8)


#----------------------S4--------------------------
#E2F targets and G2M checkpoints of PTEN-CTRL1 in S4
gs0_S4_1 = gs0_S4
gs0_S4_2 = gs0_S4

gs0_S4_1$position = 0
gs0_S4_2$position = 0

gs0_S4_1[match(gs1_S4$SYMBOL, gs0_S4_1$NAME), "position"] = 1
gs0_S4_2[match(gs2_S4$SYMBOL, gs0_S4_2$NAME), "position"] = 1

gs0_S4_1$Description = NULL
gs0_S4_2$Description = NULL

gs0_S4_1[match(gs1_S4$SYMBOL, gs0_S4_1$NAME), "Description"] = "E2F targets"
gs0_S4_2[match(gs2_S4$SYMBOL, gs0_S4_2$NAME), "Description"] = "G2M checkpoint"

gs0_S4_1$ES = 0
gs0_S4_2$ES = 0

gs0_S4_1[match(gs1_S4$SYMBOL, gs0_S4_1$NAME), "ES"] = gs1_S4$RUNNING.ES
gs0_S4_2[match(gs2_S4$SYMBOL, gs0_S4_2$NAME), "ES"] = gs2_S4$RUNNING.ES

gsdataX = rbind.data.frame(gs0_S4_1, gs0_S4_2)

gsdataY_nf2_S4 = data.frame(Name = gsdataX$NAME,
                     x = c(as.numeric(rownames(gs0_S4_1)),as.numeric(rownames(gs0_S4_2))),
                     y= as.numeric(gsdataX$SCORE),
                     Description = gsdataX$Description,
                     runningScore = gsdataX$ES,
                     position = gsdataX$position,
                     ymin = rep(0, nrow(gsdataX)),
                     ymax = rep(0, nrow(gsdataX)))

pos <- gsdataY_nf2_S4$position == 1
h <- diff(range(gsdataY_nf2_S4$runningScore))/20
gsdataY_nf2_S4$ymin[pos] <- -h
gsdataY_nf2_S4$ymax[pos] <- h

gseaPlot(gsdataY_nf2_S4, geneSetID = 1:4, title = "", color=c("red", "blue", "green", "purple"), dotsize = 4, base_size =40,
         rel_heights=c(1.5, 0.5, 0.5), subplots = 1:2, pvalue_table = FALSE, ES_geom="dot",
         y_min = -0.2, y_max = 0.6, legend.position = c(.71, .94))
ggsave("R_code/3_Diffexpr_by_CRISPR/Comparison_SKO_induced_diffexpr_in_S2_S3_S4/Enrichment_E2F_G2M_sgPTEN_S4.pdf",  width = 8.2, height = 6.8)

#--------------------- III. Enrichment of E2F_targets and G2M checkpoints of SMAD4 in S2, S3----------------------------
gs0_S2<- data.frame(read.table(file = "R_code/3_Diffexpr_by_CRISPR/S2_GSEA/SMAD4-CTRL1_vs_Control/SMAD4_hallmarks.Gsea.1603729922934/ranked_gene_list_SMAD4-CTRL1_versus_control_1603729922934.tsv",sep = "\t",header = T))
gs0_S3<- data.frame(read.table(file = "R_code/3_Diffexpr_by_CRISPRS3_GSEA/SMAD4-CTRL1_vs_Control/SMAD4_hallmarks.Gsea.1603490632019/ranked_gene_list_SMAD4-CTRL1_versus_control_1603490632019.tsv",sep = "\t",header = T))

#E2F targets
gs1_S2 <- data.frame(read.table(file = "R_code/3_Diffexpr_by_CRISPR/S2_GSEA/SMAD4-CTRL1_vs_Control/SMAD4_hallmarks.Gsea.1603729922934/HALLMARK_E2F_TARGETS.tsv",sep = "\t",header = T))
gs1_S3 <- data.frame(read.table(file = "R_code/3_Diffexpr_by_CRISPR/S3_GSEA/SMAD4-CTRL1_vs_Control/SMAD4_hallmarks.Gsea.1603490632019/HALLMARK_E2F_TARGETS.tsv",sep = "\t",header = T))

#G2M checkpoint
gs2_S2 <- data.frame(read.table(file = "R_code/3_Diffexpr_by_CRISPR/S2_GSEA/SMAD4-CTRL1_vs_Control/SMAD4_hallmarks.Gsea.1603729922934/HALLMARK_G2M_CHECKPOINT.tsv",sep = "\t",header = T))
gs2_S3 <- data.frame(read.table(file = "R_code/3_Diffexpr_by_CRISPR/S3_GSEA/SMAD4-CTRL1_vs_Control/SMAD4_hallmarks.Gsea.1603490632019/HALLMARK_G2M_CHECKPOINT.tsv",sep = "\t",header = T))

#----------------------S2--------------------------
#E2F targets and G2M checkpoints of SMAD4-CTRL1 in S2
gs0_S2_1 = gs0_S2
gs0_S2_2 = gs0_S2

gs0_S2_1$position = 0
gs0_S2_2$position = 0

gs0_S2_1[match(gs1_S2$SYMBOL, gs0_S2_1$NAME), "position"] = 1
gs0_S2_2[match(gs2_S2$SYMBOL, gs0_S2_2$NAME), "position"] = 1

gs0_S2_1$Description = NULL
gs0_S2_2$Description = NULL

gs0_S2_1[match(gs1_S2$SYMBOL, gs0_S2_1$NAME), "Description"] = "E2F targets"
gs0_S2_2[match(gs2_S2$SYMBOL, gs0_S2_2$NAME), "Description"] = "G2M checkpoint"

gs0_S2_1$ES = 0
gs0_S2_2$ES = 0

gs0_S2_1[match(gs1_S2$SYMBOL, gs0_S2_1$NAME), "ES"] = gs1_S2$RUNNING.ES
gs0_S2_2[match(gs2_S2$SYMBOL, gs0_S2_2$NAME), "ES"] = gs2_S2$RUNNING.ES

gsdataX = rbind.data.frame(gs0_S2_1, gs0_S2_2)

gsdataY_smad4_S2 = data.frame(Name = gsdataX$NAME,
                     x = c(as.numeric(rownames(gs0_S2_1)),as.numeric(rownames(gs0_S2_2))),
                     y= as.numeric(gsdataX$SCORE),
                     Description = gsdataX$Description,
                     runningScore = gsdataX$ES,
                     position = gsdataX$position,
                     ymin = rep(0, nrow(gsdataX)),
                     ymax = rep(0, nrow(gsdataX)))

pos <- gsdataY_smad4_S2$position == 1
h <- diff(range(gsdataY_smad4_S2$runningScore))/20
gsdataY_smad4_S2$ymin[pos] <- -h
gsdataY_smad4_S2$ymax[pos] <- h

gseaPlot(gsdataY_smad4_S2, geneSetID = 1:4, title = "", color=c("red", "blue", "green", "purple"), dotsize = 4, base_size =40,
         rel_heights=c(1.5, 0.5, 0.5), subplots = 1:2, pvalue_table = FALSE, ES_geom="dot",
         y_min = -0.1, y_max =0.7, legend.position = c(.68, .94))
ggsave("R_code/3_Diffexpr_by_CRISPR/Comparison_SKO_induced_diffexpr_in_S2_S3_S4/Enrichment_E2F_G2M_sgSMAD4_S2.pdf",width = 8.2, height = 6.8)


#----------------------S3--------------------------
#E2F targets and G2M checkpoints of SAMD4-CTRL1 in S3
gs0_S3_1 = gs0_S3
gs0_S3_2 = gs0_S3

gs0_S3_1$position = 0
gs0_S3_2$position = 0

gs0_S3_1[match(gs1_S3$SYMBOL, gs0_S3_1$NAME), "position"] = 1
gs0_S3_2[match(gs2_S3$SYMBOL, gs0_S3_2$NAME), "position"] = 1

gs0_S3_1$Description = NULL
gs0_S3_2$Description = NULL

gs0_S3_1[match(gs1_S3$SYMBOL, gs0_S3_1$NAME), "Description"] = "E2F targets"
gs0_S3_2[match(gs2_S3$SYMBOL, gs0_S3_2$NAME), "Description"] = "G2M checkpoint"

gs0_S3_1$ES = 0
gs0_S3_2$ES = 0

gs0_S3_1[match(gs1_S3$SYMBOL, gs0_S3_1$NAME), "ES"] = gs1_S3$RUNNING.ES
gs0_S3_2[match(gs2_S3$SYMBOL, gs0_S3_2$NAME), "ES"] = gs2_S3$RUNNING.ES

gsdataX = rbind.data.frame(gs0_S3_1, gs0_S3_2)

gsdataY_smad4_S3 = data.frame(Name = gsdataX$NAME,
                     x = c(as.numeric(rownames(gs0_S3_1)),as.numeric(rownames(gs0_S3_2))),
                     y= as.numeric(gsdataX$SCORE),
                     Description = gsdataX$Description,
                     runningScore = gsdataX$ES,
                     position = gsdataX$position,
                     ymin = rep(0, nrow(gsdataX)),
                     ymax = rep(0, nrow(gsdataX)))

pos <- gsdataY_smad4_S3$position == 1
h <- diff(range(gsdataY_smad4_S3$runningScore))/20
gsdataY_smad4_S3$ymin[pos] <- -h
gsdataY_smad4_S3$ymax[pos] <- h

gseaPlot(gsdataY_smad4_S3, geneSetID = 1:4, title = "", color=c("red", "blue", "green", "purple"), dotsize = 4, base_size =40,
         rel_heights=c(1.5, 0.5, 0.5), subplots = 1:2, pvalue_table = FALSE, ES_geom="dot",
         y_min = -0.1, y_max = 0.7, legend.position = c(.71, .94))
ggsave("R_code/3_Diffexpr_by_CRISPR/Comparison_SKO_induced_diffexpr_in_S2_S3_S4/Enrichment_E2F_G2M_sgSMAD4_S3.pdf", width = 8.2, height = 6.8)


#To get the legend
source("R_code/2_Diffexpr_by_conditions/gseaPlot.R")
gseaPlot(gsdataY_smad4_S3, geneSetID = 1:4, title = "", color=c("red", "blue", "green", "purple"), dotsize = 4, base_size =40,
         rel_heights=c(1.5, 0.5, 0.5), subplots = 1:2, pvalue_table = FALSE, ES_geom="dot",
         y_min = -0.1, y_max = 0.7, legend.position = c(.71, .94))
ggsave("R_code/3_Diffexpr_by_CRISPR/Comparison_SKO_induced_diffexpr_in_S2_S3_S4/Enrichment_E2F_G2M_sgSMAD4_S3_with legend.pdf", width = 8.2, height = 6.8)






