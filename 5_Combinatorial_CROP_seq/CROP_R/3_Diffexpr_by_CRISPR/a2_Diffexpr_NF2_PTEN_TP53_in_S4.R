#Xiaoyu Zhao 
library(ggplot2)
library(data.table) ## for read.table("xxx.xls",sep = "\t",header = T) 
library(pheatmap) 
library(Cairo) #https://stackoverflow.com/questions/27690729/greek-letters-symbols-and-line-breaks-inside-a-ggplot-legend-label

setwd("~/Downloads/5_Combinatorial_CROP_seq/")
ls = c("NF2.PTEN", "NF2.TP53", "PTEN.TP53", "NF2.CTRL1", "PTEN.CTRL1", "TP53.CTRL1")
lfc =read.csv("CROP_Python/sequencing/outs_S4/results/lfc_Deseq2_S4.csv", row.names = 2)
data = lfc[, ls]

#---------------I. Myc_targets,Oxidative_phosphorylation, mTORC1_signaling, Adipogenesis, Glycolysis, Fatty_acid_metabolism---------------------
#Import gene sets from GSEA results
gs1<- data.frame(read.table(file = "CROP_R/3_Diffexpr_by_CRISPR/S4_GSEA/NF2-PTEN_vs_Control/NF2_PTEN.Gsea.1605292391784/HALLMARK_MYC_TARGETS_V1.tsv",sep = "\t",header = T, row.names = 2))
gs2<- data.frame(read.table(file = "CROP_R/3_Diffexpr_by_CRISPR/S4_GSEA/NF2-PTEN_vs_Control/NF2_PTEN.Gsea.1605292391784/HALLMARK_OXIDATIVE_PHOSPHORYLATION.tsv",sep = "\t",header = T, row.names = 2))
gs3<- data.frame(read.table(file = "CROP_R/3_Diffexpr_by_CRISPR/S4_GSEA/NF2-PTEN_vs_Control/NF2_PTEN.Gsea.1605292391784/HALLMARK_MTORC1_SIGNALING.tsv",sep = "\t",header = T, row.names = 2))
gs4<- data.frame(read.table(file = "CROP_R/3_Diffexpr_by_CRISPR/S4_GSEA/NF2-PTEN_vs_Control/NF2_PTEN.Gsea.1605292391784/HALLMARK_ADIPOGENESIS.tsv",sep = "\t",header = T, row.names = 2))
gs5<- data.frame(read.table(file = "CROP_R/3_Diffexpr_by_CRISPR/S4_GSEA/NF2-PTEN_vs_Control/NF2_PTEN.Gsea.1605292391784/HALLMARK_GLYCOLYSIS.tsv",sep = "\t",header = T, row.names = 2))
gs6<- data.frame(read.table(file = "CROP_R/3_Diffexpr_by_CRISPR/S4_GSEA/NF2-PTEN_vs_Control/NF2_PTEN.Gsea.1605292391784/HALLMARK_FATTY_ACID_METABOLISM.tsv",sep = "\t",header = T, row.names = 2))

#Myc_targets,Oxidative_phosphorylation, mTORC1_signaling, Adipogenesis, Glycolysis, Fatty_acid_metabolism
#--------------------t.test for dat1------------------------------
dat1 = data[intersect(rownames(gs1), rownames(lfc)), ] 
dim(dat1)
#NF2.PTEN vs max(NF2.CTRL1, PTEN.CTRL1)
t.test(dat1$NF2.PTEN, dat1$NF2.CTRL1, "greater") # p = 6.416e-06
#NF2.TP53 vs max(NF2.CTRL1, PTEN.CTRL1)
t.test(dat1$NF2.TP53, dat1$NF2.CTRL1, "greater") # p = 8.034e-05
#PTEN.TP53 vs max(PTEN.CTRL1, TP53.CTRL1)
t.test(dat1$PTEN.TP53, dat1$PTEN.CTRL1, "greater") # p = 4.779e-06

#--------------------t.test for dat2------------------------------
dat2 = data[intersect(rownames(gs2), rownames(lfc)), ]
dim(dat2)
#NF2.PTEN vs max(NF2.CTRL1, PTEN.CTRL1)
t.test(dat2$NF2.PTEN, dat2$PTEN.CTRL1, "greater") # p = 0.005994
#NF2.TP53 vs max(NF2.CTRL1, PTEN.CTRL1)
t.test(dat2$NF2.TP53, dat2$NF2.CTRL1, "greater") # p = 0.1426
#PTEN.TP53 vs max(PTEN.CTRL1, TP53.CTRL1)
t.test(dat2$PTEN.TP53, dat2$PTEN.CTRL1, "greater") # p = 1.091e-05

#--------------------t.test for dat3------------------------------
dat3 = data[intersect(rownames(gs3), rownames(lfc)), ]
dim(dat3)
#NF2.PTEN vs max(NF2.CTRL1, PTEN.CTRL1)
t.test(dat3$NF2.PTEN, dat3$NF2.CTRL1, "greater") # p = 0.004432
#NF2.TP53 vs max(NF2.CTRL1, PTEN.CTRL1)
t.test(dat3$NF2.TP53, dat3$NF2.CTRL1, "greater") # p = 0.484
#PTEN.TP53 vs max(PTEN.CTRL1, TP53.CTRL1)
t.test(dat3$PTEN.TP53, dat3$PTEN.CTRL1, "greater") # p = 0.0009195

#--------------------t.test for dat4------------------------------
dat4 = data[intersect(rownames(gs4), rownames(lfc)), ]
dim(dat4)
#NF2.PTEN vs max(NF2.CTRL1, PTEN.CTRL1)
t.test(dat4$NF2.PTEN, dat4$PTEN.CTRL1, "greater") # p = 0.005808
#NF2.TP53 vs max(NF2.CTRL1, PTEN.CTRL1)
t.test(dat4$NF2.TP53, dat4$NF2.CTRL1, "greater") # p = 0.03337
#PTEN.TP53 vs max(PTEN.CTRL1, TP53.CTRL1)
t.test(dat4$PTEN.TP53, dat4$PTEN.CTRL1, "greater") # p = 0.008549

#--------------------t.test for dat5------------------------------
dat5 = data[intersect(rownames(gs5), rownames(lfc)), ]
dim(dat5)
#NF2.PTEN vs max(NF2.CTRL1, PTEN.CTRL1)
t.test(dat5$NF2.PTEN, dat5$NF2.CTRL1, "greater") # p = 0.1061
#NF2.TP53 vs max(NF2.CTRL1, PTEN.CTRL1)
t.test(dat5$NF2.TP53, dat5$NF2.CTRL1, "greater") # p =  0.276
#PTEN.TP53 vs max(PTEN.CTRL1, TP53.CTRL1)
t.test(dat5$PTEN.TP53, dat5$PTEN.CTRL1, "greater") # p = 0.01856

#--------------------t.test for dat6------------------------------
dat6 = data[intersect(rownames(gs6), rownames(lfc)), ]
dim(dat6)
#NF2.PTEN vs max(NF2.CTRL1, PTEN.CTRL1)
t.test(dat6$NF2.PTEN, dat6$PTEN.CTRL1, "greater") # p = 0.101
#NF2.TP53 vs max(NF2.CTRL1, PTEN.CTRL1)
t.test(dat6$NF2.TP53, dat6$NF2.CTRL1, "greater") # p =  0.2956
#PTEN.TP53 vs max(PTEN.CTRL1, TP53.CTRL1)
t.test(dat6$PTEN.TP53, dat6$PTEN.CTRL1, "greater") # p = 0.1962


##Barplot of mean and SEM of Myc_targets,Oxidative_phosphorylation, mTORC1_signaling, Adipogenesis, Glycolysis, Fatty_acid_metabolism
df1 = cbind.data.frame(group =c("NF2-PTEN", "NF2-TP53", "PTEN-TP53", "NF2-CTRl1", "PTEN-CTRL1", "TP53-CTRL1"), 
                 hallmark = rep("Myc_targets", 6),
                 mean = apply(dat1, 2, function(x) mean(x)), 
                 sd = apply(dat1, 2, function(x) sd(x)),
                 se = apply(dat1, 2, function(x) sd(x)/sqrt(nrow(dat1)))
)


df2 = cbind.data.frame(group =c("NF2-PTEN", "NF2-TP53", "PTEN-TP53", "NF2-CTRl1", "PTEN-CTRL1", "TP53-CTRL1"), 
                       hallmark = rep("Oxidative_phosphorylation", 6),
                       mean = apply(dat2, 2, function(x) mean(x)), 
                       sd = apply(dat2, 2, function(x) sd(x)),
                       se = apply(dat2, 2, function(x) sd(x)/sqrt(nrow(dat2)))
)

df3 = cbind.data.frame(group =c("NF2-PTEN", "NF2-TP53", "PTEN-TP53", "NF2-CTRl1", "PTEN-CTRL1", "TP53-CTRL1"), 
                       hallmark = rep("mTORC1_signaling", 6),
                       mean = apply(dat3, 2, function(x) mean(x)), 
                       sd = apply(dat3, 2, function(x) sd(x)),
                       se = apply(dat3, 2, function(x) sd(x)/sqrt(nrow(dat3)))
)

df4 = cbind.data.frame(group =c("NF2-PTEN", "NF2-TP53", "PTEN-TP53", "NF2-CTRl1", "PTEN-CTRL1", "TP53-CTRL1"), 
                       hallmark = rep("Adipogenesis", 6),
                       mean = apply(dat4, 2, function(x) mean(x)), 
                       sd = apply(dat4, 2, function(x) sd(x)),
                       se = apply(dat4, 2, function(x) sd(x)/sqrt(nrow(dat4)))
)

df5 = cbind.data.frame(group =c("NF2-PTEN", "NF2-TP53", "PTEN-TP53", "NF2-CTRl1", "PTEN-CTRL1", "TP53-CTRL1"), 
                       hallmark = rep("Glycolysis", 6),
                       mean = apply(dat5, 2, function(x) mean(x)), 
                       sd = apply(dat5, 2, function(x) sd(x)),
                       se = apply(dat5, 2, function(x) sd(x)/sqrt(nrow(dat5)))
)

df6 = cbind.data.frame(group =c("NF2-PTEN", "NF2-TP53", "PTEN-TP53", "NF2-CTRl1", "PTEN-CTRL1", "TP53-CTRL1"), 
                       hallmark = rep("Fatty_acid_metabolism", 6),
                       mean = apply(dat6, 2, function(x) mean(x)), 
                       sd = apply(dat6, 2, function(x) sd(x)),
                       se = apply(dat6, 2, function(x) sd(x)/sqrt(nrow(dat6)))
)


d1 = rbind.data.frame(df1, df2, df3, df4, df5, df6)

groups = factor(d1$group, levels = c("NF2-PTEN", "NF2-TP53", "PTEN-TP53", "NF2-CTRl1", "PTEN-CTRL1", "TP53-CTRL1"))
ggplot(d1, aes(x=groups, y=mean, fill = hallmark)) + 
  geom_bar(position=position_dodge(), stat="identity",colour = "grey20", width = 0.85) +
  geom_errorbar(aes(ymin=mean, ymax=mean+se),
                width=.3,# Width of the error bars
                color = "grey20",
                position=position_dodge(0.85))+
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        axis.line = element_line(colour = "grey50"),
        # legend.position = c(2,0.9),
        #legend.position = "None",
        legend.title = element_text(size = 12, color = "black"),
        legend.text = element_text(size=12 , color = "black"),
        axis.text.y  =element_text(size=12, color = "black"),
        axis.text.x =element_text(size=12, angle = 45, color = "black", vjust = 0.8, hjust = 0.8, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title=element_text(size=12),
        plot.title = element_text(size=36,hjust=0.5)
        #strip.text.y  = element_text(size = 26, face="bold")
        ##strip.text is for facet
  )+
  #scale_fill_manual(values=c('#999999','#E69F00'))+
  scale_fill_brewer(direction = -1)+
  scale_y_continuous(limits = c(-0.03,0.12),breaks = seq(0, 0.12, by = 0.04))+
  labs(title="",x ="", y = "LFCs")
ggsave("CROP_R/3_Diffexpr_by_CRISPR/Myc_OXPHOS_mTORC1_Adipogenesis_Glycolysis_Fattyacid_in_3_DKOs.pdf", width = 6.6, height = 4)


#------------------------------------------II. E2F_targets, G2M_checkpoint------------------------------------------------
gs7<- data.frame(read.table(file = "CROP_R/3_Diffexpr_by_CRISPR/S4_GSEA/NF2-CTRL1_vs_Control/NF2_hallmarks.Gsea.1603490996038/HALLMARK_E2F_TARGETS.tsv",sep = "\t",header = T, row.names = 2))
gs8<- data.frame(read.table(file = "CROP_R/3_Diffexpr_by_CRISPR/S4_GSEA/NF2-CTRL1_vs_Control/NF2_hallmarks.Gsea.1603490996038/HALLMARK_G2M_CHECKPOINT.tsv",sep = "\t",header = T, row.names = 2))

#--------------------t.test for dat7------------------------------
dat7 = data[intersect(rownames(gs7), rownames(lfc)), ]
dim(dat7)
#NF2.PTEN vs max(NF2.CTRL1, PTEN.CTRL1)
t.test(dat7$NF2.PTEN, dat7$NF2.CTRL1, "greater") # p = 0.0001818
#NF2.TP53 vs max(NF2.CTRL1, PTEN.CTRL1)
t.test(dat7$NF2.TP53, dat7$NF2.CTRL1, "greater") # p =  0.02563
#PTEN.TP53 vs max(PTEN.CTRL1, TP53.CTRL1)
t.test(dat7$PTEN.TP53, dat7$PTEN.CTRL1, "greater") # p = 0.001827

#--------------------t.test for dat8------------------------------
dat8 = data[intersect(rownames(gs8), rownames(lfc)), ]
dim(dat8)
#NF2.PTEN vs max(NF2.CTRL1, PTEN.CTRL1)
t.test(dat8$NF2.PTEN, dat8$NF2.CTRL1, "greater") # p = 0.00293
#NF2.TP53 vs max(NF2.CTRL1, PTEN.CTRL1)
t.test(dat8$NF2.TP53, dat8$NF2.CTRL1, "greater") # p =  8.283e-05
#PTEN.TP53 vs max(PTEN.CTRL1, TP53.CTRL1)
t.test(dat8$PTEN.TP53, dat8$PTEN.CTRL1, "greater") # p = 0.0003278

df7 = cbind.data.frame(group =c("NF2-PTEN", "NF2-TP53", "PTEN-TP53", "NF2-CTRl1", "PTEN-CTRL1", "TP53-CTRL1"), 
                       hallmark = rep("E2F_targets", 6),
                       mean = apply(dat7, 2, function(x) mean(x)), 
                       sd = apply(dat7, 2, function(x) sd(x)),
                       se = apply(dat7, 2, function(x) sd(x)/sqrt(nrow(dat7)))
)

df8 = cbind.data.frame(group =c("NF2-PTEN", "NF2-TP53", "PTEN-TP53", "NF2-CTRl1", "PTEN-CTRL1", "TP53-CTRL1"), 
                       hallmark = rep("G2M_checkpoint", 6),
                       mean = apply(dat8, 2, function(x) mean(x)), 
                       sd = apply(dat8, 2, function(x) sd(x)),
                       se = apply(dat8, 2, function(x) sd(x)/sqrt(nrow(dat8)))
)

d2 = rbind.data.frame(df7, df8)
groups = factor(d2$group, levels = c("NF2-PTEN", "NF2-TP53", "PTEN-TP53", "NF2-CTRl1", "PTEN-CTRL1", "TP53-CTRL1"))
ggplot(d2, aes(x=groups, y=mean, fill = hallmark)) + 
  geom_bar(position=position_dodge(),width=.5, stat="identity", color = "grey20") +
  geom_errorbar(aes(ymin=mean, ymax=mean+se),
                width=.2,                    # Width of the error bars
                color = "grey20",
                position=position_dodge(.5))+
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        axis.line = element_line(colour = "grey60"),
        #legend.position = c(2,0.9),
        legend.position = "None",
        legend.title = element_blank(),
        legend.text = element_text(size=12 , color = "black"),
        axis.text.y  =element_text(size=12, color = "black"),
        axis.text.x =element_text(size=12, angle = 45, color = "black", vjust = 0.8, hjust = 0.8, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title=element_text(size=12),
        plot.title = element_text(size=36,hjust=0.5)
        #strip.text.y  = element_text(size = 26, face="bold")
        ##strip.text is for facet
  )+
  scale_fill_brewer(direction = -1)+
  scale_y_continuous(limits = c(-0.03,0.12),breaks = seq(0, 0.12, by = 0.04))+
  labs(title="",x ="", y = "LFCs")
ggsave("CROP_R/3_Diffexpr_by_CRISPR/E2F_and_G2M_in_3_DKOs.pdf", width =4, height = 4)





