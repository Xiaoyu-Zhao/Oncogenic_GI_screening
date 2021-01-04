#Xiaoyu Zhao 
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(VennDiagram)
library(tidyr)
setwd("~/Downloads/data/")

#----------------------I. Comparison of tumor promoting effects in three genetic contexts---------------------
#Tumor-promoting effects of single genes
ref <- data.frame(read_excel("TableS1_Annotatation_of_asssyed_TSGs_or_TSG_combinations.xlsx", 3))
gene52 = as.character(ref$gene_symbol)
#Import results from TableS3
res.pten <- data.frame(read_excel("TableS3_Raw_and_processed_data_combinatorial_screening_in_vivo.xlsx", 4))
rownames(res.pten) = res.pten$gene_names
res.pik <- data.frame(read_excel("TableS3_Raw_and_processed_data_combinatorial_screening_in_vivo.xlsx", 6))
rownames(res.pik) = res.pik$gene_names
res.myc <- data.frame(read_excel("TableS3_Raw_and_processed_data_combinatorial_screening_in_vivo.xlsx", 8))
rownames(res.myc) = res.myc$gene_names


singleEffect = matrix(0, nrow = length(gene52), ncol = 3)
rownames(singleEffect) = gene52
colnames(singleEffect) = c("PTEN-/-", "PIK3CA", "HAMyc")

nM.pten = rownames(res.pten[which(res.pten$LFC > 0 & res.pten$p_adj < 0.05 ),])
nM.pik3 = rownames(res.pik[which(res.pik$LFC > 0 & res.pik$p_adj < 0.05 ),])
nM.myc = rownames(res.myc[which(res.myc$LFC > 0 & res.myc$p_adj < 0.05 ),])

for(i in 1:nrow(singleEffect)){
  singleEffect[i,3] = sum(grepl(gene52[i], nM.pten))
  singleEffect[i,2] = sum(grepl(gene52[i], nM.pik3))
  singleEffect[i,1] = sum(grepl(gene52[i], nM.myc))
}

#sE1 = singleEffect[which(rowSums(singleEffect) > 1),]
sE = singleEffect[which(rowSums(singleEffect) > 1 
                        & (singleEffect[,1] >1 | singleEffect[,2] >1 | singleEffect[,3] >1)),]
dim(sE)
sE0 = sE[order(rowSums(sE), decreasing = T),]
write.csv(sE0, "Tumor_promoting_effects_52_single_genes_in_tumors.csv")
#----------------------------------------------------------------------------------------------------
res0 = read.csv("Effects_top8_single_gene_in_tumors.csv", row.names = 1, header = T)
colnames(res0) = c("HAMyc", "PIK3CA", "PTEN-/-")
#categories_col = brewer.pal(10, "Set3")
my.breaks <- c(seq(0, 40, by=1)) 
#colfunc <- colorRampPalette(c("white", "#007AFF","mediumblue", "royalblue4"))
colfunc <- colorRampPalette(c( "white","dodgerblue1", "dodgerblue2","dodgerblue3", "dodgerblue4"))
plot(rep(1,5),col=colfunc(5),pch=19,cex=3)
#my.colors <- c(colorRampPalette(colors = c("yellow", "white"))(length(seq(-1, 0, by=0.1))), colorRampPalette(colors = c("white", "purple"))(length(seq(0.1, 1, by=0.1))))
my.colors <- colfunc(length(my.breaks))
my_sample_col <- data.frame(Tumors = c("HAMyc", "PIK3CA", "PTEN-/-"))
rownames(my_sample_col) <- colnames(res0)
#categories_col = colorRampPalette(c("red","yellow","springgreen","royalblue"))(50)
categories_col = colorRampPalette(c("red","yellow","springgreen","purple"))(50)
colorRampPalette(c("red","yellow","springgreen","purple"))(50)
plot(rep(1,50),col = categories_col,pch=19,cex=3)
my_gene_col <- data.frame(Categories = c(rep("Background_indepedent",3), "Favored_in_PTEN-/-",
                                         "Favored_in_HAMyc", rep("No_effects_in_PTEN-/-",3)))
rownames(my_gene_col) <- rownames(res0)

#pathways_col = colorRampPalette(c("red","yellow","springgreen","royalblue"))(50)
pdf("Heatmap_top8_single_gene_in_tumors_20200908.pdf", width = 4.5, height = 3.8)
#par(mar=c(0,20,0,5))
pheatmap(res0, angle_col = 45, fontsize_row =14, fontsize_col = 14,
         mar=c(0,20,0,5),
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = my.colors, breaks = my.breaks, 
         border_color = "grey", legend_labels = TRUE,
         #annotation_col = my_sample_col, 
         annotation_row  = my_gene_col,
         annotation_names_row = FALSE)
dev.off()


#-------------------------II.Comparison of GI networks in three genetic contexts --------------------------
#Nodes: Centrality score of GI network in tumors
cen_pten = read.csv("PTEN_tumors/Centrality_scores_PTEN.csv", row.names = 1, header = T)
dim(cen_pten)
cen_pik3 = read.csv("PIK3CA_tumors/Centrality_scores_PIK3CA.csv", row.names = 1, header = T)
dim(cen_pik3)
cen_myc = read.csv("HAMyc_tumors/Centrality_scores_HAMyc.csv", row.names = 1, header = T)
dim(cen_myc)

cen = union(union(rownames(cen_pten), rownames(cen_pik3)), rownames(cen_myc))
unNION = data.frame(matrix(0, length(cen), 3)) 
colnames(unNION) = c("PTEN-/-", "PIK3CA", "HAMyc")
rownames(unNION) = cen
unNION[match(rownames(cen_pten), rownames(unNION)), "PTEN-/-"] = cen_pten[,"degree"]
unNION[match(rownames(cen_pik3), rownames(unNION)), "PIK3CA"] = cen_pik3[,"degree"]
unNION[match(rownames(cen_myc), rownames(unNION)), "HAMyc"] = cen_myc[,"degree"]
unNION
unNION$sum = rowSums(unNION[,1:3])


#com = apply(unNION, 1, function(x) sum(x!=0))
#df = unNION[which(com>1),]
unNION = unNION[order(unNION$sum, decreasing = T),]
df =unNION[,1:3]

df$sname = rownames(df)
#df$sum = rowSums(df[,1:3])
data = gather(df, condition, cen, -sname)
data$cen[data$cen == 0] = NA

snames = factor(data$sname, levels= rownames(unNION[order(-unNION$sum),]))
Conditions = factor(data$condition, levels = c("HAMyc", "PIK3CA", "PTEN-/-"))

ggplot(data, aes(x = Conditions, y = snames, size= cen**2 ))+ 
  geom_point(color = "#007AFF")+ scale_size_area(max_size = 12)+
  #scale_fill_gradient(color = "blue")+
  theme(panel.border = element_rect(fill=NA, size=0.2, color = "black"), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        axis.line = element_line(size=0.1, color = "black"),
        legend.position = "none",
        legend.title = element_text(size = 20, face= "bold"),
        legend.text = element_text(size=20 ),
        axis.text.y  = element_text(size=18,  color = 'black' ),
        axis.text.x =element_text(size=16, color = 'black',vjust=0, angle = 90),
        axis.title.x=element_text(size=16, color = 'black', vjust=-1.5, margin = margin(t =0, r =0, b =1, l=0, unit = "cm")),#face="bold")
        plot.title = element_text(size=26,hjust= 1, vjust = 3)
        #strip.text.y  = element_text(size = 26, face="bold")
        ##strip.text is for facet
  )+
  scale_y_discrete(position = "right")+
  labs(title="",x ="", y = "")+
  coord_flip()

ggsave("Centrality_score.pdf", width = 8.5, height = 4.8)

#Edges: Overlapped Oncogenic GIs
#-------------------------------------------III. Overlapped Oncogenic GIs -------------------------------------------------
gi_pten = read.csv("PTEN_tumors/18_Tumor_promoting_GIs_PTEN.csv", row.names = 1, header = T)
dim(gi_pten)
gi_pik3 = read.csv("PIK3CA_tumors/20_Tumor_promoting_GIs_PIK3CA.csv", row.names = 1, header = T)
dim(gi_pik3)
gi_myc = read.csv("HAMyc_tumors/19_tumor_promoting_GIs_HAMyc.csv", row.names = 1, header = T)
dim(gi_myc)

nGI = union(union(rownames(gi_pten), rownames(gi_pik3)), rownames(gi_myc))

unNION = data.frame(matrix(0, length(nGI), 3)) 
colnames(unNION) = c("PTEN-/-", "PIK3CA", "HAMyc")
rownames(unNION) = nGI

unNION[match(rownames(gi_pten), rownames(unNION)), "PTEN-/-"] = 1
unNION[match(rownames(gi_pik3), rownames(unNION)), "PIK3CA"] = 1
unNION[match(rownames(gi_myc), rownames(unNION)), "HAMyc"] = 1

overLap<- function(tumor.type) {
  ppl <- unNION
  
  for (i in 1:length(tumor.type)) {
    ppl <- subset(ppl, ppl[tumor.type[i]] == T)
  }
  nrow(ppl)
}

a = c("PTEN-/-", "PIK3CA", "HAMyc")

pdf("Overlapped_GIs.pdf", height = 6, width = 6)
grid.newpage()
overrideTriple=T  
##https://stackoverflow.com/questions/11727068/scaling-triple-venn-diagram-in-r-with-venndiagram-package
draw.triple.venn(overLap(a[1]), overLap(a[2]), overLap(a[3]), overLap(a[1:2]), 
                 overLap(a[2:3]), overLap(a[c(1,3)]), overLap(a),
                 #category = c("PTEN-/-", "PIK3CA", "HAMyc"), 
                 lty = "blank", 
                 fill = c("#5AC8FA", "pink1", "mediumorchid"),
                 scaled=TRUE,fontfamily ="Helvetica", cat.fontfamily="Helvetica",
                 cex = 2, cat.cex = 2)
#To change the main font and category font
#https://stackoverflow.com/questions/56583545/how-to-change-venn-diagram-font
dev.off()

unNION[which(rowSums(unNION) == 3),]
unNION[which(rowSums(unNION) == 1),]
unNION[which(rowSums(unNION) == 2),]










