#Xiaoyu Zhao
RNGkind(sample.kind = "Rounding")
# RNGkind is to get the results from set.seed() to match for different R versions
library(ggplot2)
library(pheatmap)
library(gplots)
library(RColorBrewer)
library(Cairo)
library(dplyr)
library(readxl)
setwd("~/Downloads/data/")
ref <- data.frame(read_excel("TableS1_Annotatation_of_asssyed_TSGs_or_TSG_combinations.xlsx", 5))
gene52 = as.character(unique(ref[,1]))
#----------------------------------------------@@@@@@@@@@@@@@@@@@@---------------------------------------------------------
#==============================PART 1: Comparison among results under three culture conditions=============================
#----------------------------------------------@@@@@@@@@@@@@@@@@@@---------------------------------------------------------
#-------------------------------------I. Distribution of sgRNA_level fitness in 3 conditions--------------------------------------------------------
##Full medium
r.f1 = read.csv("full/fitness/full.1.RPM_fitlm.csv", row.names =1, header = T)
dim(r.f1)
r.f2 = read.csv("full/fitness//full.2.RPM_fitlm.csv", row.names =1, header = T)
dim(r.f2)
#all(rownames(r.f1) == rownames(r.f2))
coM1 = intersect(rownames(r.f1), rownames(r.f2))
f.full= data.frame(f1 = r.f1[coM1,dim(r.f1)[2]], f2 = r.f2[coM1,dim(r.f2)[2]])
f.full$f.avg = apply(f.full, 1, function(x) mean(x[1:2]))
f.full$conditions = rep("Full", nrow(f.full))
##-------------------------------------------------------------------------------------------------------------------------
##Minimal medium
r.m1 = read.csv("minimal/fitness/min.1.RPM_fitlm.csv", row.names =1, header = T)
dim(r.m1)
r.m2 = read.csv("minimal/fitness/min.2.RPM_fitlm.csv", row.names =1, header = T)
dim(r.m2)
#all(rownames(r.m1) == rownames(r.m2))
coM2 = intersect(rownames(r.m1), rownames(r.m2))
f.min = data.frame(f1 = r.m1[coM2,dim(r.m1)[2]], f2 = r.m2[coM2,dim(r.m2)[2]])
f.min$f.avg = apply(f.min, 1, function(x) mean(x[1:2]))
f.min$conditions = rep("Minimal", nrow(f.min))
##-------------------------------------------------------------------------------------------------------------------------
##+ TGFbeta 1 medium
r.t1 = read.csv("tgf/fitness/tgf.1.RPM_fitlm.csv", row.names =1, header = T)
dim(r.t1)
r.t2 = read.csv("tgf/fitness/tgf.2.RPM_fitlm.csv", row.names =1, header = T)
dim(r.t2)
coM3 = intersect(rownames(r.t1), rownames(r.t2))
f.tgf = data.frame(f1 = r.t1[coM3,dim(r.t1)[2]], f2 = r.t2[coM3,dim(r.t2)[2]])
f.tgf$f.avg = apply(f.tgf, 1, function(x) mean(x[1:2]))
f.tgf$conditions = rep(paste0("TGF", "\u03B2" , "-1"), nrow(f.tgf))

m1 = median(f.full[,"f.avg"])
m2 = median(f.min[,"f.avg"])
m3 = median(f.tgf[,"f.avg"])
#commonS = intersect(rownames(r1),rownames(r2))
identical(names(f.full), names(f.min))
files = list(f.full, f.min, f.tgf)
dat1 = data.frame(do.call(rbind,files))
dim(dat1)

Conditions = factor(dat1$conditions, levels = c(paste0("TGF", "\u03B2" , "-1"),"Minimal", "Full"))
 p = ggplot(dat1, aes(x=f.avg, fill=Conditions)) +
  geom_histogram(binwidth=.004, alpha=.5, position="identity")+
  scale_fill_manual(values=c("blue", "green", "magenta"))+
  theme_set(theme_bw())+
  theme(panel.border = element_rect(fill=NA, size=1, color = "black"), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line = element_line(colour = "black"),
    legend.position = c(0.8,0.7),
    legend.title = element_text(size = 20),
    legend.text = element_text(size=20 ),
    axis.text.y  =element_text(size=24,color= "black"),
    axis.text.x =element_text(size=24, color = "black", margin = margin(c(t = 5, r = 0, b = 10, l = 0))),
    axis.title=element_text(size=30),#face="bold")
    plot.title = element_text(size=32,hjust=0.5)
  )+
  labs(title="",x ="Fitness", y = "Frequency")+
  geom_vline(xintercept = c(-0.1, 0, 0.1), linetype="dotted", color = c("grey","black","grey"), size=1.5)
p
ggsave("comparison_in_vitro/Distribution_sgRNA_level_fitness_3_conditions.pdf", width = 6, height = 6.5, device = cairo_pdf)

#------------------------------------II. Top 5 gene_level fitness in 3 conditions---------------------------------------------------
f1 = data.frame(read.csv("full/fitness/full_SKO_fitness_gene_level.csv", row.names = 1, header = T))
f2 = data.frame(read.csv("minimal/fitness/min_SKO_fitness_gene_level.csv", row.names = 1, header = T))
f3 = data.frame(read.csv("tgf/fitness/tgf_SKO_fitness_gene_level.csv", row.names = 1, header = T))

f1$conditions =  rep("Full", nrow(f1))
f2$conditions =  rep("Minimal", nrow(f2))
f3$conditions =  rep(paste0("TGF", "\u03B2" , "-1"), nrow(f3))

f1_order = f1[order(f1$f.avg, decreasing = T),]
f1_order$order = 1:53
f2_order = f2[order(f2$f.avg, decreasing = T),]
f2_order$order = 1:53
f3_order = f3[order(f3$f.avg, decreasing = T),]
f3_order$order = 1:53

dat2 = rbind.data.frame(f1_order[1:5,], f2_order[1:5,], f3_order[1:5,])
dat2$Conditions = factor(dat2$conditions, levels = c("Full", "Minimal", paste0("TGF", "\u03B2" , "-1")))

p <- ggplot(dat2,aes(x=order, y=f.avg))+
  geom_bar(position=position_dodge(), stat="identity", fill = "#5AC8FA", color ="#007AFF", width = 0.8)+
  geom_errorbar(aes(ymin=l, ymax=u), width=.2, color = "grey20")+
  scale_y_continuous(limits = c(-0.07, 0.25),
                     breaks = round(seq(-0.1, 0.25, by =0.05),2),
                     labels = round(seq(-0.1, 0.25, by = 0.05),2))+
  labs(title="",x ="TOP 5 SKOs", y = "Fitness")+
  scale_fill_manual(values = c("magenta", "green", "blue"))+
  theme(panel.border = element_rect(color ="black", fill = NA, size =0.5), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line = element_line(colour = "black"),
    legend.position = "none",
    #legend.title = element_text(size = 20, face= "bold"),
    #legend.text = element_text(size=18 ),
    axis.text.y = element_text(size=14, color = 'black'),
    axis.text.x =element_text(size=16, color = 'black',angle = 45, vjust =1, hjust = 1),
    axis.title =  element_text(size=14, color = 'black'),
    plot.title = element_text(size=18,hjust=0.5),
    strip.text  = element_text(size = 14, color = 'black'),
    strip.background = element_rect(color = "black", size = 0.5)
  )+
  facet_grid(.~Conditions)
  
p
ggsave("Comparison_in_vitro/TOP5_SKO_genes_3_conditions.pdf", device = cairo_pdf, width = 6.2, height= 3.5)
#------------------------------------------III. Proportion of 3 catergories_dual_sgRNAs---------------------------------------------------------
f1 = read.csv("full/fitness/full_relative_RPM_percentage_SKO_DKO_CTRL.csv", row.names = 1, header = T)
f2 = read.csv("minimal/fitness/min_relative_RPM_percentage_SKO_DKO_CTRL.csv", row.names = 1, header = T)
f3 = read.csv("tgf/fitness/tgf_relative_RPM_percentage_SKO_DKO_CTRL.csv", row.names = 1, header = T)

#Take the data at t0 and t6
full.t1 = data.frame(type = c("SKOs", "DKOs", "CTRLs"), mean = f1[1:3,"t1"], se = f1[4:6,"t1"], Conditions = rep("Full",3), t = rep("T1", 3))
full.t6 = data.frame(type = c("SKOs", "DKOs", "CTRLs"), mean = f1[1:3,"t6"], se = f1[4:6,"t6"], Conditions = rep("Full",3), t = rep("T6", 3) )
min.t1 = data.frame(type = c("SKOs", "DKOs", "CTRLs"), mean = f2[1:3,"t1"], se = f2[4:6,"t1"], Conditions = rep("Minimal",3), t = rep("T1", 3) )
min.t6 = data.frame(type = c("SKOs", "DKOs", "CTRLs"), mean = f2[1:3,"t6"], se = f2[4:6,"t6"], Conditions = rep("Minimal",3), t = rep("T6", 3) )
tgf.t1 = data.frame(type = c("SKOs", "DKOs", "CTRLs"), mean = f3[1:3,"t1"], se = f3[4:6,"t1"], Conditions = rep(paste0("+ TGF ", "\u03B2" , "1"),3), t = rep("T1", 3) )
tgf.t6 = data.frame(type = c("SKOs", "DKOs", "CTRLs"), mean = f3[1:3,"t6"], se = f3[4:6,"t6"], Conditions = rep(paste0("+ TGF ", "\u03B2" , "1"),3), t = rep("T6", 3) )

#Combine all the data into a single dataframe
data = rbind.data.frame(full.t1, full.t6, min.t1, min.t6, tgf.t1, tgf.t6)
data$l = data$mean - data$se
data$u = data$mean + data$se
data$l[data$type == "SKOs"] = with(data,l[type == "SKOs"] + mean[type == "CTRLs"])
data$u[data$type == "SKOs"] = with(data,u[type == "SKOs"] + mean[type == "CTRLs"])

data$l[data$type == "DKOs"] = with(data,l[type == "DKOs"] + mean[type == "SKOs"]+ mean[type == "CTRLs"])
data$u[data$type == "DKOs"] = with(data,u[type == "DKOs"] + mean[type == "SKOs"] + mean[type == "CTRLs"])

p <- ggplot(data, aes(x=t, y=mean, fill= factor(type, levels = c("DKOs","SKOs","CTRLs")))) +
  geom_bar(stat="identity", position = "stack", alpha = 0.8, width = 0.7)+
  geom_errorbar(aes(ymin=l, ymax=u), width=.2)+
  #geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, stat="identity",position = 'identity')+
  facet_grid(.~Conditions)+
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        axis.line = element_line(colour = "black"),
        #legend.position = c(1,0.8),
        legend.title = element_text(size = 17),
        legend.text = element_text(size=17),
        axis.text.y  =element_text(size=20,color = 'black'),
        axis.text.x =element_text(size=20,color = 'black',margin= margin(c(5,0,0,0))),
        axis.title=element_text(size=23),#face="bold")
        #plot.title = element_text(size=32,hjust=0.5),
        strip.text  = element_text(size = 21, color = "black")
  )+
  scale_y_continuous(limits = c(0,103), breaks = seq(0,100,25))+
  labs(title="",x ="Time points", y = paste("Relative percentage (%)"))+
  scale_fill_manual(values=c("magenta","#5AC8FA","grey"),
                    name="",
                    #breaks=c("Cells", "Tumors"),
                    labels=c("DKOs", "SKOs","CTRLs"))
p
ggsave('Comparison_in_vitro/Proportion of 3 catergories_dual_sgRNAs.pdf',device = cairo_pdf, width = 7, height = 4.5)

#----------------------------------------------@@@@@@@@@@@@@@@@@@@---------------------------------------------------------
#==================================PART 2: Comparison of results in vitro and in vivo======================================
#----------------------------------------------@@@@@@@@@@@@@@@@@@@---------------------------------------------------------
#------------------------------I. Heatmap comparing single gene effects in vitro and in vivo-------------------------------
#Import fitness results in vitro
f1 = data.frame(read.csv("full/fitness/full_SKO_fitness_gene_level.csv", row.names = 1, header = T))
f2 = data.frame(read.csv("minimal/fitness/min_SKO_fitness_gene_level.csv", row.names = 1, header = T))
f3 = data.frame(read.csv("tgf/fitness/tgf_SKO_fitness_gene_level.csv", row.names = 1, header = T))
#Import LFC results in vivo from TableS3
res.pten <- data.frame(read_excel("TableS3_Raw_and_processed_data_combinatorial_screening_in_vivo.xlsx", 4))
rownames(res.pten) = res.pten$gene_names
res.pik <- data.frame(read_excel("TableS3_Raw_and_processed_data_combinatorial_screening_in_vivo.xlsx", 6))
rownames(res.pik) = res.pik$gene_names
res.myc <- data.frame(read_excel("TableS3_Raw_and_processed_data_combinatorial_screening_in_vivo.xlsx", 8))
rownames(res.myc) = res.myc$gene_names

#Combine single-gene effects in vitro and in vivo
SKOs = paste0(gene52, "-CTRL")
res = cbind(f1[SKOs,"f.avg"], f2[SKOs,"f.avg"],  f3[SKOs,"f.avg"],
           res.pten[SKOs,"LFC"], res.pik[SKOs,"LFC"], res.myc[SKOs,"LFC"])
colnames(res) = c("Full", "Minimal", paste("TGF", "\u03B2" , "-1"), "PTEN-/-", "PIK3CA", "MYC")
rownames(res) = gene52
class(res)
sgMatrix = apply(res, 2, function(x) scale(x))
colnames(sgMatrix) = c("Full", "Minimal", paste("TGF", "\u03B2" , "-1"), "PTEN-/-", "PIK3CA", "MYC")
rownames(sgMatrix) = gene52

CairoPDF(file = paste0("Comparison_in_vitro_in_vivo/Single_gene_effects_invitro_invivo.pdf"), width = 3.0, height = 6.2)

heatmap.2(as.matrix(sgMatrix),
          Rowv = TRUE, 
          Colv = TRUE,
          na.rm = TRUE,
          hclustfun = function(x) hclust(x,method = "complete"),
          dendrogram="both",
          cexRow = 0.55, cexCol = 1,
          scale = "none", 
          key=TRUE, keysize=0.5, symkey=TRUE,
          density.info="none", trace="none",
          col=colorRampPalette(c("blue", "white", "magenta")),
          breaks = seq(-4,4,0.1),
          margins=c(8,8),
          colRow = rep("black",52)
)

dev.off()

#Scale color for fitness 
df <- reshape2::melt(outer(1:2,1:2), varnames = c("X1", "X2"))
ggplot(df, aes(X1, X2, fill = value))+ 
  geom_raster()+ 
  scale_fill_gradientn(name = "z.score",
                       colours=colorRampPalette(c("blue", "white", "magenta"))(length(seq(-4, 4, 0.1))), na.value = "transparent",
                       breaks=seq(-4,4,2),labels=seq(-4,4,2),
                       limits=c(-4,4))+
  guides(size = 3, fill = guide_colourbar(ticks = TRUE))+
  theme(legend.text = element_text(size=14),
        legend.title = element_text(size = 14))
ggsave("Comparison_in_vitro_in_vivo/scale_color_Single_gene_effects_invitro_invivo.pdf", width = 3, height =3)

#--------------------------------------------II. Overlapped_Nodes_invitro_and_invivo---------------------------------------------------------
#nodes in GI networks in vivo
cen_pten = read.csv("PTEN_tumors/Centrality_scores_PTEN.csv", row.names = 1, header = T)
dim(cen_pten)
cen_pik3 = read.csv("PIK3CA_tumors/Centrality_scores_PIK3CA.csv", row.names = 1, header = T)
dim(cen_pik3)
cen_myc = read.csv("HAMyc_tumors/Centrality_scores_HAMyc.csv", row.names = 1, header = T)
dim(cen_myc)

#nodes in GI networks in vitro
nodes_full = c("NF2", "PTEN", "TP53", "RB1", "MLL4", "GATA3", "KMT2C")    #7 nodes
nodes_min = c("NF2", "PTEN", "TP53", "RB1","CBFB","SMAD4","NF1", "RUNX1") #8 nodes

#How many nodes detected in full medium were found in GI networks in vivo
n1 = sum(nodes_full %in% rownames(cen_pten))/(length(nodes_full)-1)*100 #"PTEN" was removed in denominator when compared to PTEN-/- tumors
n2 = sum(nodes_full %in% rownames(cen_pik3))/(length(nodes_full))*100
n3 = sum(nodes_full %in% rownames(cen_myc))/(length(nodes_full))*100
#How many nodes detected in full medium were found in GI networks in vivo
n4 = sum(nodes_min %in% rownames(cen_pten))/(length(nodes_min)-1)*100 #PTEN" was removed in denominator when compared to PTEN-/- tumors
n5 = sum(nodes_min %in% rownames(cen_pik3))/(length(nodes_min))*100
n6 = sum(nodes_min %in% rownames(cen_myc))/(length(nodes_min))*100

nodes = data.frame(Groups = rep(c("PTEN-/-", "PIK3CA", "MYC"),2),
                   Conditions = c(rep("Full",3), rep("Minimal", 3)),
                   Fraction = c(n1, n2, n3, n4, n5, n6))

p <- ggplot(nodes,aes(x=factor(Conditions, levels = c("Minimal", "Full")), y= Fraction,fill = factor(Groups,levels = c("PTEN-/-", "PIK3CA", "MYC"))))+
  geom_bar(position=position_dodge(), stat="identity",  width = 0.8,color="grey10")+
  theme(#panel.border = element_blank(), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line = element_line(colour = "black"),
    #legend.position = "none",
    legend.title = element_text(size = 16,color = 'black'),
    legend.text = element_text(size=14 ),
    axis.text.y  =element_text(size=16, color = 'black'),
    axis.text.x =element_text(size=16, color = 'black'),
    axis.title=element_text(size=18, color = 'black'),#face="bold")
    plot.title = element_text(size=18,hjust=0.5),
    strip.text  = element_text(size = 16, color = 'black')
  )+  
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100, by =25))+
  scale_fill_manual(values=c("orchid1","#FFCC00","#5AC8FA"), 
                        name=expression(paste(italic("In vivo")," groups")),
                        labels=c("PTEN-/-", "PIK3CA", "MYC"))+
  labs(title="",x ="", y = expression(paste("% of nodes discovered ",italic("in vivo"))))
p
ggsave('Comparison_in_vitro_in_vivo/Overlapped_Nodes_invitro_and_invivo.pdf', width = 4.6, height = 3.5)

#--------------------------------------III. Overlapped_Edges_invitro_and_invivo----------------------------------------------
#Edges in GI networks in vivo
gi_pten = read.csv("PTEN_tumors/18_Tumor_promoting_GIs_PTEN.csv", row.names = 1, header = T)
dim(gi_pten)
gi_pik3 = read.csv("PIK3CA_tumors/20_Tumor_promoting_GIs_PIK3CA.csv", row.names = 1, header = T)
dim(gi_pik3)
gi_myc = read.csv("HAMyc_tumors/19_tumor_promoting_GIs_HAMyc.csv", row.names = 1, header = T)
dim(gi_myc)
#Edges in GI networks in vitro
gi_full = read.csv("full/GI/Oncogenic_GIs_cytoscape_full.csv" ,row.names = 1, header = T)
gi_min = read.csv("minimal/GI/Oncogenic_GIs_cytoscape_min.csv",row.names = 1, header = T)

#How many edges detected in full medium were found in GI networks in vivo
n1 = sum(rownames(gi_full) %in% rownames(gi_pten))/(nrow(gi_full)-2)*100 #"PTEN" linked GIs were removed in denominator when compared to PTEN-/- tumors
n2 = sum(rownames(gi_full) %in% rownames(gi_pik3))/(nrow(gi_full))*100
n3 = sum(rownames(gi_full) %in% rownames(gi_myc))/(nrow(gi_full))*100
#How many edges detected in full medium were found in GI networks in vivo
n4 = sum(rownames(gi_min) %in% rownames(gi_pten))/(nrow(gi_min)-4)*100 #"PTEN" linked GIs were removed in denominator when compared to PTEN-/- tumors
n5 = sum(rownames(gi_min)  %in% rownames(gi_pik3))/(nrow(gi_min))*100
n6 = sum(rownames(gi_min)  %in% rownames(gi_myc))/(nrow(gi_min))*100

edges = data.frame(Groups = rep(c("PTEN-/-", "PIK3CA", "MYC"),2),
                   Conditions = c(rep("Full",3), rep("Minimal", 3)),
                   Fraction = c(n1, n2, n3, n4, n5, n6))

p <- ggplot(edges,aes(x=factor(Conditions, levels = c("Minimal", "Full")), y= Fraction, fill = factor(Groups,levels = c("PTEN-/-", "PIK3CA", "MYC"))))+
  geom_bar(position=position_dodge(), stat="identity",  width = 0.8,color="grey10")+
  theme(#panel.border = element_blank(), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line = element_line(colour = "black"),
    #legend.position = "none",
    legend.title = element_text(size = 16,color = 'black'),
    legend.text = element_text(size=14 ),
    axis.text.y  =element_text(size=16, color = 'black'),
    axis.text.x =element_text(size=16, color = 'black'),
    axis.title=element_text(size=18, color = 'black'),#face="bold")
    plot.title = element_text(size=18,hjust=0.5),
    strip.text  = element_text(size = 16, color = 'black')
  )+  
  scale_y_continuous(limits = c(0, 80),breaks = seq(0,80, by =20))+
  scale_fill_manual(values=c("orchid1","#FFCC00","#5AC8FA"),
                    name=expression(paste(italic("In vivo")," groups")),
                    labels=c("PTEN-/-", "PIK3CA", "MYC"))+
  labs(title="",x ="", y = expression(paste("% of edges discovered ",italic("in vivo"))))
p
ggsave('Comparison_in_vitro_in_vivo/Overlapped_Edges_invitro_and_invivo.pdf', width = 4.6, height = 3.5)

#--------------------------------------IV. Overlapped_Nodes_and Edges_invitro_and_invivo----------------------------------------------
d1 = data.frame(Groups = nodes$Groups, num = nodes$Fraction, Conditions = nodes$Conditions, Category = "Node")
d2=  data.frame(Groups = edges$Groups, num = edges$Fraction, Conditions = edges$Conditions, Category = "Edge")
d = rbind.data.frame(d1, d2)

ggplot(d,aes(x=Category, y= num, fill= factor(Conditions,levels = c("Minimal", "Full"))))+
  geom_boxplot(position=position_dodge(1), width = 0.7)+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1, position=position_dodge(1))+
  theme(panel.border = element_blank(), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line = element_line(colour = "black"),
    #legend.position = "none",
    legend.title = element_text(size = 16, color = "black"),
    legend.text = element_text(size=16, color = "black"),
    axis.text.y  =element_text(size=16, color = 'black'),
    axis.text.x =element_text(size=16, color = 'black', margin = margin(c(5,0,0,0))),
    axis.title=element_text(size=20, color = 'black'),#face="bold")
    plot.title = element_text(size=18,hjust=0.5),
    strip.text  = element_text(size = 16, color = 'black')
  )+
  scale_y_continuous(limits = c(10,100),breaks = seq(0,100, by =20))+
  scale_fill_manual(values=c("magenta","#5AC8FA"),
                    name="Conditions",
                    labels=c("Minimal", "Full"))+
  labs(title="",x ="", y = expression(paste("% discovered ",italic("in vivo"))))

#statistical comparison 
t.test(d$num[1:3], d$num[4:6], alternative = "less")
t.test(d$num[7:9], d$num[10:12], alternative = "less")

ggsave('Comparison_in_vitro_in_vivo/Overlapped_Nodes_Edges_invitro_and_invivo.pdf', width = 4.4, height = 4.5)


                 
                              