#Xiaoyu Zhao 20200908
library(ggplot2)
library(data.table) ## for read.table("xxx.xls",sep = "\t",header = T) 
library(pheatmap) 
library(Cairo) #https://stackoverflow.com/questions/27690729/greek-letters-symbols-and-line-breaks-inside-a-ggplot-legend-label
library(VennDiagram)
library(ggpubr)
library(tidyr)
library(UpSetR)
setwd("~/Downloads/5_Combinatorial_CROP_seq/")
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#--------------------- I. Hallmarks induced by GF removal: Kras_DN and TNFα_signaling_via NFKB--------------------------
gs0<- data.frame(read.table(file = "CROP_R/2_Diffexpr_by_conditions/S4_vs_S2/S4_vs_S2_phenotype_permutation.Gsea.1602035227339/ranked_gene_list_S4_versus_S2_1602035227339.tsv",sep = "\t",header = T))
gs1<- data.frame(read.table(file = "CROP_R/2_Diffexpr_by_conditions/S4_vs_S2/S4_vs_S2_phenotype_permutation.Gsea.1602035227339/HALLMARK_KRAS_SIGNALING_DN.tsv",sep = "\t",header = T))
gs2<- data.frame(read.table(file = "CROP_R/2_Diffexpr_by_conditions/S4_vs_S2/S4_vs_S2_phenotype_permutation.Gsea.1602035227339/HALLMARK_TNFA_SIGNALING_VIA_NFKB.tsv",sep = "\t",header = T))

gs0_1 = gs0
gs0_2 = gs0

gs0_1$position = 0
gs0_2$position = 0

gs0_1[match(gs1$SYMBOL, gs0_1$NAME), "position"] = 1
gs0_2[match(gs2$SYMBOL, gs0_2$NAME), "position"] = 1

gs0_1$Description = NULL
gs0_2$Description = NULL

gs0_1[match(gs1$SYMBOL, gs0_1$NAME), "Description"] = "Kras_signaling_DN"
gs0_2[match(gs2$SYMBOL, gs0_2$NAME), "Description"] = "TNFα_signaling_via_NFκB"

gs0_1$ES = 0
gs0_2$ES = 0

gs0_1[match(gs1$SYMBOL, gs0_1$NAME), "ES"] = gs1$RUNNING.ES
gs0_2[match(gs2$SYMBOL, gs0_2$NAME), "ES"] = gs2$RUNNING.ES

gsdataX = rbind.data.frame(gs0_1, gs0_2)

gsdataY = data.frame(Name = gsdataX$NAME,
                     x = c(as.numeric(rownames(gs0_1)),as.numeric(rownames(gs0_2))),
                     y= as.numeric(gsdataX$SCORE),
                     Description = gsdataX$Description,
                     runningScore = gsdataX$ES,
                     position = gsdataX$position,
                     ymin = rep(0, nrow(gsdataX)),
                     ymax = rep(0, nrow(gsdataX)))

pos <- gsdataY$position == 1
h <- diff(range(gsdataY$runningScore))/20
gsdataY$ymin[pos] <- -h
gsdataY$ymax[pos] <- h

source("CROP_R/2_Diffexpr_by_conditions/gseaPlot.R")
gseaPlot(gsdataY, geneSetID = 1:4, title = "", color=c("red", "blue", "green", "purple"), dotsize = 4, base_size =40,
         rel_heights=c(1.5, 0.5, 0.5), subplots = 1:2, pvalue_table = FALSE, ES_geom="dot",
         y_min = -0.2, y_max =0.66, legend.position = c(.5, .99))
ggsave("CROP_R/2_Diffexpr_by_conditions/S4_vs_S2/gseaPlot_S4_vs_S2_enriched.pdf", width = 8.2, height = 6.8, device = cairo_pdf)


#-------------------II. Hallmarks downregulated by GF removal: Myc_targets, MTORC1_signaling, Glycolysis, Oxidative_phosphorylation, Adipogenesis and Fatty_acid_metabolish----------------------------
gs0<- data.frame(read.table(file = "CROP_R/2_Diffexpr_by_conditions/S4_vs_S2/S4_vs_S2_phenotype_permutation.Gsea.1602035227339/ranked_gene_list_S4_versus_S2_1602035227339.tsv",sep = "\t",header = T))
gs1<- data.frame(read.table(file = "CROP_R/2_Diffexpr_by_conditions/S4_vs_S2/S4_vs_S2_phenotype_permutation.Gsea.1602035227339/HALLMARK_MYC_TARGETS_V1.tsv",sep = "\t",header = T))
gs2<- data.frame(read.table(file = "CROP_R/2_Diffexpr_by_conditions/S4_vs_S2/S4_vs_S2_phenotype_permutation.Gsea.1602035227339/HALLMARK_MTORC1_SIGNALING.tsv",sep = "\t",header = T))
gs3<- data.frame(read.table(file = "CROP_R/2_Diffexpr_by_conditions/S4_vs_S2/S4_vs_S2_phenotype_permutation.Gsea.1602035227339/HALLMARK_GLYCOLYSIS.tsv",sep = "\t",header = T))
gs4<- data.frame(read.table(file = "CROP_R/2_Diffexpr_by_conditions/S4_vs_S2/S4_vs_S2_phenotype_permutation.Gsea.1602035227339/HALLMARK_OXIDATIVE_PHOSPHORYLATION.tsv",sep = "\t",header = T))
gs5<- data.frame(read.table(file = "CROP_R/2_Diffexpr_by_conditions/S4_vs_S2/S4_vs_S2_phenotype_permutation.Gsea.1602035227339/HALLMARK_ADIPOGENESIS.tsv",sep = "\t",header = T))
gs6<- data.frame(read.table(file = "CROP_R/2_Diffexpr_by_conditions/S4_vs_S2/S4_vs_S2_phenotype_permutation.Gsea.1602035227339/HALLMARK_FATTY_ACID_METABOLISM.tsv",sep = "\t",header = T))

gs0_1 = gs0
gs0_2 = gs0
gs0_3 = gs0
gs0_4 = gs0
gs0_5 = gs0
gs0_6 = gs0

gs0_1$position = 0
gs0_2$position = 0
gs0_3$position = 0
gs0_4$position = 0
gs0_5$position = 0
gs0_6$position = 0

gs0_1[match(gs1$SYMBOL, gs0_1$NAME), "position"] = 1
gs0_2[match(gs2$SYMBOL, gs0_2$NAME), "position"] = 1
gs0_3[match(gs3$SYMBOL, gs0_3$NAME), "position"] = 1
gs0_4[match(gs4$SYMBOL, gs0_4$NAME), "position"] = 1
gs0_5[match(gs5$SYMBOL, gs0_5$NAME), "position"] = 1
gs0_6[match(gs6$SYMBOL, gs0_6$NAME), "position"] = 1

gs0_1$Description = NULL
gs0_2$Description = NULL
gs0_3$Description = NULL
gs0_4$Description = NULL
gs0_5$Description = NULL
gs0_6$Description = NULL

gs0_1[match(gs1$SYMBOL, gs0_1$NAME), "Description"] = "Myc targets"
gs0_2[match(gs2$SYMBOL, gs0_2$NAME), "Description"] = "mTORC1_signaling"
gs0_3[match(gs3$SYMBOL, gs0_3$NAME), "Description"] = "Glycolysis"
gs0_4[match(gs4$SYMBOL, gs0_4$NAME), "Description"] = "OXPHOS"
gs0_5[match(gs5$SYMBOL, gs0_5$NAME), "Description"] = "Adipogenesis"
gs0_6[match(gs6$SYMBOL, gs0_6$NAME), "Description"] = "Fatty_acid_metabolism"

gs0_1$ES = 0
gs0_2$ES = 0
gs0_3$ES = 0
gs0_4$ES = 0
gs0_5$ES = 0
gs0_6$ES = 0

gs0_1[match(gs1$SYMBOL, gs0_1$NAME), "ES"] = gs1$RUNNING.ES
gs0_2[match(gs2$SYMBOL, gs0_2$NAME), "ES"] = gs2$RUNNING.ES
gs0_3[match(gs3$SYMBOL, gs0_3$NAME), "ES"] = gs3$RUNNING.ES
gs0_4[match(gs4$SYMBOL, gs0_4$NAME), "ES"] = gs4$RUNNING.ES
gs0_5[match(gs5$SYMBOL, gs0_5$NAME), "ES"] = gs5$RUNNING.ES
gs0_6[match(gs6$SYMBOL, gs0_6$NAME), "ES"] = gs6$RUNNING.ES

gsdataX = rbind.data.frame(gs0_1, gs0_2, gs0_3, gs0_4, gs0_5, gs0_6)

gsdataY = data.frame(Name = gsdataX$NAME,
                     x = c(as.numeric(rownames(gs0_1)),as.numeric(rownames(gs0_2)), as.numeric(rownames(gs0_3)), as.numeric(rownames(gs0_5)), as.numeric(rownames(gs0_5)), as.numeric(rownames(gs0_6))),
                     y= as.numeric(gsdataX$SCORE),
                     Description = gsdataX$Description,
                     runningScore = gsdataX$ES,
                     position = gsdataX$position,
                     ymin = rep(0, nrow(gsdataX)),
                     ymax = rep(0, nrow(gsdataX)))

pos <- gsdataY$position == 1
h <- diff(range(gsdataY$runningScore))/20
gsdataY$ymin[pos] <- -h
gsdataY$ymax[pos] <- h

source("CROP_R/2_Diffexpr_by_conditions/gseaPlot.R")
gseaPlot(gsdataY, geneSetID = 1:4, title = "", color=c("red", "blue", "green", "purple","pink","yellow"), dotsize = 4, base_size =40,
         rel_heights=c(1.5, 0.5, 0.5), subplots = 1:2, pvalue_table = FALSE, ES_geom="dot",
         y_min = -0.8, y_max = 0.4, legend.position = c(.56, .9))
ggsave("CROP_R/2_Diffexpr_by_conditions/S4_vs_S2/gseaPlot_S4_vs_S2_underrepresented.pdf", width = 8.2, height = 6.8, device = cairo_pdf)

#------------------------- III. Scatter plot of LFCs induced by GF removal and rescued by sgRNAs-----------------------------
setwd("~/Downloads/5_Combinatorial_CROP_seq/")
ls = c("NF2", "SMAD4", "PTEN","NF1", "RB1", "TP53", "CBFB", "CDH1", "CASP8", "TBX3", "USP9X" )
ref = read.csv("CROP_Python/sequencing/outs_S4/results/res_Deseq2/res_S4_vs_S2.csv", row.names = 1)
lfc =read.csv("CROP_Python/sequencing/outs_S4/results/lfc_Deseq2_S4.csv", row.names = 2)
comM1 = intersect(rownames(ref), rownames(lfc))
read = data.frame(S4_vs_S2 = ref[comM1, "LFC"], lfc[comM1,])
write.csv(read, "CROP_R/2_Diffexpr_by_conditions/S4_vs_S2/All_DEGs_induced_by_GF_removal_and_recovered_by_perturbations.csv")
#Scatter plot of all 11 SKOs
threshold = 0.1
for (sg in ls){
  pb = paste0(sg, ".CTRL1")
  data = data.frame(x = read[[pb]], y = read$S4_vs_S2)
  df = data[abs(data$y) > threshold,]
  df$DEGs = ifelse(df$x*df$y < 0, "Yes", "No")
  #scatter plots
  q<-ggscatter(df, x = "x" , y = "y", size = 2, color = "DEGs",
            add = "reg.line", conf.int = FALSE, 
            add.params = list(color = "black", linetype = "dashed"), # Customize reg. line
            cor.coef = TRUE,rug = FALSE,
            cor.coeff.args = list(method = "pearson", label.x = 0.1,  label.y =3.8, label.sep = "\n", size = 8.5, color = "blue"))+
          theme(panel.border = element_rect(colour = "black",size = 1,fill=NA), 
          #panel.grid = element_blank(),
          panel.background = element_blank(),
          #panel.grid= element_line(linetype = "dashed"),
          axis.line = element_line(colour = "black"),
          legend.position = "none",
          axis.text.y  =element_text(size=30),
          axis.text.x =element_text(size=30,margin = margin(t = 1, r = 1, b = 1, l = 1)),
          axis.title=element_text(size=30),#face="bold")
          plot.title = element_text(face="bold", size=20,hjust=0.5) #
          #strip.text.y  = element_text(size = 12, face="bold")
          ##strip.text is for facet
    )+labs(title="",x =paste0("LFCs recovered\n by sg", sg), y = paste0("LFCs induced\n by GF removal"))+
    scale_x_continuous(limits = c(-2,2), breaks = seq(-2, 2, by = 1)) +
    scale_y_continuous(limits = c(-5,5),breaks = seq(-4, 4, by = 2))+
    geom_hline(yintercept=0, linetype="dashed", color = "grey")+
    geom_vline(xintercept=0, linetype="dashed", color = "grey")+
    scale_color_manual(values = c('grey', 'blue'))
  q
  path0 = "CROP_R/2_Diffexpr_by_conditions/S4_vs_S2/"
  ggsave(q, file = paste(path0, "DEs_", sg, "_S4",".pdf", sep =""), width = 5.4, height = 5.4)
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#----------------------------------IV. Heatmap of SKOs on selective removal_of_GFs induced DEGs----------------------------------------
read = data.frame(S4_vs_S2 = ref[comM1, "LFC"], lfc[comM1,])
ls = c("NF2", "PTEN","CBFB","SMAD4","NF1","TP53","CDH1","RB1","CASP8", "TBX3", "USP9X" )
SKOs = paste0(ls, ".CTRL1")
readM = read[,c("S4_vs_S2",SKOs, c("NF2.PTEN", "NF2.TP53", "PTEN.TP53"))]
colnames(readM) = c("S4_vs_S2","NF2", "PTEN", "CBFB", "SMAD4", "NF1", "TP53", "CASP8", "CDH1", "RB1", "TBX3", "USP9X", "NF2-PTEN", "NF2-TP53", "PTEN-TP53")

threshold = 0.1
#pheatmap(readM)
readM_L = readM[which(abs(readM[,1])  > threshold ), ]
dim(readM_L) # 2871 genes

#Identify which gene expression were reversed by each sgRNA
a = rownames(readM_L[which(abs(readM_L[,2]) > 0.1 & readM_L[,2]*readM_L[,1] < 0),])
b = rownames(readM_L[which(abs(readM_L[,3]) > 0.1 & readM_L[,3]*readM_L[,1] < 0),])
c = rownames(readM_L[which(abs(readM_L[,4]) > 0.1 & readM_L[,4]*readM_L[,1] < 0),])
d = rownames(readM_L[which(abs(readM_L[,5]) > 0.1 & readM_L[,5]*readM_L[,1] < 0),])
e = rownames(readM_L[which(abs(readM_L[,6]) > 0.1 & readM_L[,6]*readM_L[,1] < 0),])
f = rownames(readM_L[which(abs(readM_L[,7]) > 0.1 & readM_L[,7]*readM_L[,1] < 0),])
g = rownames(readM_L[which(abs(readM_L[,8]) > 0.1 & readM_L[,8]*readM_L[,1] < 0),])
h = rownames(readM_L[which(abs(readM_L[,9]) > 0.1 & readM_L[,9]*readM_L[,1] < 0),])
i = rownames(readM_L[which(abs(readM_L[,10]) > 0.1 & readM_L[,10]*readM_L[,1] < 0),])
j = rownames(readM_L[which(abs(readM_L[,11]) > 0.1 & readM_L[,11]*readM_L[,1] < 0),])
k = rownames(readM_L[which(abs(readM_L[,12]) > 0.1 & readM_L[,12]*readM_L[,1] < 0),])

comM2 = unique(unlist(c(a, b, c, d, e, f, g, h, i, j, k)))
length(comM2)  #401 growth-factor related genes were reversed by sgRNAs

#Enriched hallmarks in CTRL cells of S4
gs1<- data.frame(read.table(file = "CROP_R/2_Diffexpr_by_conditions/S4_vs_S2/S4_vs_S2_phenotype_permutation.Gsea.1602035227339/HALLMARK_TNFA_SIGNALING_VIA_NFKB.tsv",sep = "\t",header = T, row.names = 2))
gs2<- data.frame(read.table(file = "CROP_R/2_Diffexpr_by_conditions/S4_vs_S2/S4_vs_S2_phenotype_permutation.Gsea.1602035227339/HALLMARK_KRAS_SIGNALING_DN.tsv",sep = "\t",header = T, row.names = 2))

#Enriched hallmarks in CTRL cells of S2
gs3<- data.frame(read.table(file = "CROP_R/2_Diffexpr_by_conditions/S4_vs_S2/S4_vs_S2_phenotype_permutation.Gsea.1602035227339/HALLMARK_MYC_TARGETS_V1.tsv",sep = "\t",header = T, row.names = 2))
gs4<- data.frame(read.table(file = "CROP_R/2_Diffexpr_by_conditions/S4_vs_S2/S4_vs_S2_phenotype_permutation.Gsea.1602035227339/HALLMARK_MTORC1_SIGNALING.tsv",sep = "\t",header = T, row.names = 2))
gs5<- data.frame(read.table(file = "CROP_R/2_Diffexpr_by_conditions/S4_vs_S2/S4_vs_S2_phenotype_permutation.Gsea.1602035227339/HALLMARK_GLYCOLYSIS.tsv",sep = "\t",header = T, row.names = 2))
gs6<- data.frame(read.table(file = "CROP_R/2_Diffexpr_by_conditions/S4_vs_S2/S4_vs_S2_phenotype_permutation.Gsea.1602035227339/HALLMARK_OXIDATIVE_PHOSPHORYLATION.tsv",sep = "\t",header = T, row.names = 2))
gs7<- data.frame(read.table(file = "CROP_R/2_Diffexpr_by_conditions/S4_vs_S2/S4_vs_S2_phenotype_permutation.Gsea.1602035227339/HALLMARK_ADIPOGENESIS.tsv",sep = "\t",header = T, row.names = 2))
gs8<- data.frame(read.table(file = "CROP_R/2_Diffexpr_by_conditions/S4_vs_S2/S4_vs_S2_phenotype_permutation.Gsea.1602035227339/HALLMARK_FATTY_ACID_METABOLISM.tsv",sep = "\t",header = T, row.names = 2))

#Extract genes significally differentially expressed in DESea analysis and enriched in GSEA analysis
readM_gs1 = readM[intersect(comM2, rownames(gs1[gs1$CORE.ENRICHMENT == "Yes",])),]
readM_gs2 = readM[intersect(comM2, rownames(gs2[gs2$CORE.ENRICHMENT == "Yes",])),]
readM_gs3 = readM[intersect(comM2, rownames(gs3[gs3$CORE.ENRICHMENT == "Yes",])),]
readM_gs4 = readM[intersect(comM2, rownames(gs4[gs4$CORE.ENRICHMENT == "Yes",])),]
readM_gs5 = readM[intersect(comM2, rownames(gs5[gs5$CORE.ENRICHMENT == "Yes",])),]
readM_gs6 = readM[intersect(comM2, rownames(gs6[gs6$CORE.ENRICHMENT == "Yes",])),]
readM_gs7 = readM[intersect(comM2, rownames(gs7[gs7$CORE.ENRICHMENT == "Yes",])),]
readM_gs8 = readM[intersect(comM2, rownames(gs8[gs8$CORE.ENRICHMENT == "Yes",])),]

#Combine LFCs of these genes together
cut.pos = 0.3
cut.neg = -0.3
df = rbind(readM_gs1[readM_gs1$S4_vs_S2 > cut.pos,],
           readM_gs2[readM_gs2$S4_vs_S2 > cut.pos,],
           readM_gs4["NUPR1",],
           readM_gs3[readM_gs3$S4_vs_S2 < cut.neg,],
           readM_gs4[readM_gs4$S4_vs_S2 < cut.neg,],
           readM_gs5[readM_gs5$S4_vs_S2 < cut.neg,],
           readM_gs6[readM_gs6$S4_vs_S2 < cut.neg,],
           readM_gs7[readM_gs7$S4_vs_S2 < cut.neg,],
           readM_gs8[readM_gs8$S4_vs_S2 < cut.neg,]
)

df$pathways = rep(c("TNFα_signaling_via_NFκB", "Kras_signaling_DN", "Others", "Myc_targets", "mTORC1_signaling", "Glycolysis", "Oxidative_phosphorylation", "Adipogenesis","Fatty_acid_metabolism"),
                  c(length(which(readM_gs1$S4_vs_S2 > cut.pos)),
                    length(which(readM_gs2$S4_vs_S2 > cut.pos)),
                    1,
                    length(which(readM_gs3$S4_vs_S2 < cut.neg)),
                    length(which(readM_gs4$S4_vs_S2 < cut.neg)),
                    length(which(readM_gs5$S4_vs_S2 < cut.neg)),
                    length(which(readM_gs6$S4_vs_S2 < cut.neg)),
                    length(which(readM_gs7$S4_vs_S2 < cut.neg)),
                    length(which(readM_gs8$S4_vs_S2 < cut.neg))
                  ))

dim(df)

write.csv(df, "CROP_R/2_Diffexpr_by_conditions/S4_vs_S2/DEGs_in_different_pathways_induced_by_GF_removal_and_recovered_by_perturbations.csv")
#Parameters for heatmap
#1.annotation for columns
my_sample_col1 <- data.frame(Groups = rep(c("DEGs_induced_by_removal_of_GFs", "DEGs_recovered_by_perturbations"), c(1,14)))
my_sample_col2 <- data.frame(Clusters = rep(c("N", "Cluster_1", "Cluster_2", "Cluster_3", "Not_assigned", "Cluster_1"), c(1,1,1,7,2,3)))
my_sample_col = data.frame(Clusters = my_sample_col2, Groups =  my_sample_col1)
rownames(my_sample_col) <- colnames(readM)

#2.annotation for rows
my_gene_col <- data.frame(Pathways = df$pathways)
rownames(my_gene_col) <- rownames(df)

#3.colors for LFCs (heatmap)
#my.breaks <- c(seq(-1, 0, by=0.1), seq(0.1, 1, by=0.1)) 
my.breaks <- c(seq(-1, 1, by=0.1)) 
colfunc <- colorRampPalette(c("orange", "white", "purple"))
#my.colors <- c(colorRampPalette(colors = c("yellow", "white"))(length(seq(-1, 0, by=0.1))), colorRampPalette(colors = c("white", "purple"))(length(seq(0.1, 1, by=0.1))))
my.colors <- colfunc(length(my.breaks))

#4. annotation colors
pathways_col = colorRampPalette(c(c("red","yellow","springgreen","royalblue")))(50)

cairo_pdf(file = paste0("CROP_R/2_Diffexpr_by_conditions/S4_vs_S2/heatmap_DEGs_58genes_SKO_3DKO.pdf"),  width = 10.4, height =2.8)
pheatmap(as.matrix(t(df[,1:15])), angle_col = 90, fontsize_row =9, fontsize_col = 8.5,
         cluster_rows = FALSE, cluster_cols = FALSE, 
         color = my.colors, breaks = my.breaks, 
         border_color = NA,legend_labels = TRUE,
         annotation_row = my_sample_col, 
         annotation_col  = my_gene_col,
         annotation_colors = list(Groups = c(DEGs_induced_by_removal_of_GFs ="Red",DEGs_recovered_by_perturbations="blue"),
                                  Clusters = c( N = "white", Cluster_1 = "#a157db", Cluster_2 = "#db5f57", Cluster_3 = "#57d3db", Not_assigned = "#cccccc"),
                                  Pathways = c(TNFα_signaling_via_NFκB = pathways_col[6],
                                               Kras_signaling_DN = pathways_col[10],  
                                               Others = pathways_col[14], 
                                               Myc_targets = pathways_col[20],
                                               mTORC1_signaling = pathways_col[25],
                                               Glycolysis = pathways_col[36],
                                               Oxidative_phosphorylation = pathways_col[40], 
                                               Adipogenesis = pathways_col[45],Fatty_acid_metabolism = pathways_col[50]
                                  )))
dev.off()

#----------------------------------V. Upset diagrame showing the overlapping DEGs reversed by sgRNAs----------------------------------------
readM_S0 = readM_L[comM2,1:12]
readM_S1 = readM_S0[which(readM_S0$S4_vs_S2 < 0), 2:12]
dim(readM_S1)
readM_S2 = readM_S0[which(readM_S0$S4_vs_S2 > 0), 2:12]
dim(readM_S2)

readM_S = readM_S0[,2:12]
dim(readM_S)
#readM_S1_M = matrix(0, nrow(readM_S1), ncol(readM_S1))
#rownames(readM_S1_M) = rownames(readM_S1)
#rownames(readM_S1_M) = rownames(readM_S1)
readM_S[abs(readM_S) > 0.1 & readM_S*readM_S0[,1] < 0] = 1
readM_S[readM_S != 1] = 0
readM_S$sum = rowSums(readM_S)
#readM_S1_M[readM_S1_M >= -0.1] = 0

cairo_pdf("CROP_R/2_Diffexpr_by_conditions/S4_vs_S2/Intersected_DEGs_recovered_by_sgRNAs_S4.pdf", width =6.8, height =4.8)
upset(readM_S, nsets = 11, keep.order = T, nintersects = NA, matrix.color = "gray23", sets.bar.color = "#007AFF", sets.x.label = "No. of removal_of_GFs\n induced DEGs\nrecovered by sgRNAs",
      main.bar.color = "#007AFF", point.size = 1.8, line.size = 0.8, set_size.show = T, set_size.scale_max=340, text.scale = c(1.8,2,1.45,1.5,1.5,1.5), mb.ratio = c(0.5,0.5))
dev.off()


