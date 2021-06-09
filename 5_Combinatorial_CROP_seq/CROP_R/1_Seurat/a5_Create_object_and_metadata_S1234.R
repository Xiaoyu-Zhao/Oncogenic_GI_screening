##Xiaoyu Zhao  20200908
library(Seurat)
library(Matrix)
library(ggplot2)
setwd("~/Downloads/5_Combinatorial_CROP_seq/CROP_R/1_Seurat/")
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#--------------------- I. Merge four objects into S1234 ----------------------------
gex.S1 <- readRDS( file = "objects/gex.S1.rds") #dim(gex.S1) = (16735, 5864)
gex.S2 <- readRDS( file = "objects/gex.S2.rds") #dim(gex.S2) = (16575, 7844)
gex.S3 <- readRDS( file = "objects/gex.S3.rds") #dim(gex.S3) = (16975, 6511)
gex.S4 <- readRDS( file = "objects/gex.S4.rds") #dim(gex.S4) = (16915, 8364)
gex.S1234 <- merge(x = gex.S1, y = c(gex.S2, gex.S3, gex.S4), add.cell.ids = c("S1", "S2","S3","S4"), project = "S1234")
dim(gex.S1234) #(18556, 28583) 5864 + 7844 + 6511 + 8364 = 28583
head(gex.S1234@meta.data)
VlnPlot(object = gex.S1234, features = c("nCount_RNA","nFeature_RNA","percent_mito"), ncol = 3, pt.size = 0.1)
FeatureScatter(object = gex.S1234, feature1 = "nCount_RNA", feature2 = "percent_mito")
FeatureScatter(object = gex.S1234, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
write.csv(as.matrix(gex.S1234@meta.data),file="cell_identities_S1234.csv")
saveRDS(gex.S1234, file = "objects/gex.S1234_nofilter_28583cells.rds")

#-------------------------------------------------------------------------------------
gex.S1234 <- subset(x=gex.S1234, subset = nFeature_RNA > 800 & nFeature_RNA < 5000 & percent_mito < 0.1) #According to violin plot of nFeature_RNA
dim(gex.S1234)
gex.S1234 <- NormalizeData(gex.S1234,normalization.method = "LogNormalize", scale.factor = 1e4)
gex.S1234 <- FindVariableFeatures(gex.S1234, selection.method = 'vst', nfeatures = 3500)
length(x = VariableFeatures(object = gex.S1234))
VariableFeaturePlot(gex.S1234)
gex.S1234 <- ScaleData(object = gex.S1234)
gex.S1234 <- RunPCA(gex.S1234, features = VariableFeatures(object = gex.S1234), verbose = FALSE)
#print(x = gex.S1234[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)
#VizDimLoadings(object = gex.S1234, dims = 1:2)
ElbowPlot(gex.S1234,ndims = 50, reduction = "pca")
#-------------------------------------------------------------------------------------
#gex.S1234 <- FindNeighbors(gex.S1234, dims.use = 1:30, print.output = FALSE)
#gex.S1234 <- FindClusters(gex.S1234, resolution = 0.6)
#gex.S1234 <- RunTSNE(gex.S1234, dims = 1:20)
gex.S1234 <- RunUMAP(gex.S1234, dims = 1:20,  n.neighbors = 30L, min.dist = 2)
#RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
saveRDS(gex.S1234, file = "objects/gex.S1234_filter_26081cells.rds")
#-------------------------------------------------------------------------------------
pdf(file = "UMAP_S1234.pdf",10,8)
#DimPlot(object = gex.S1234, reduction = 'tsne', group.by = "Condition",label = FALSE,pt.size = 0.5, label.size = 10) + NoAxes()
DimPlot(object = gex.S1234, reduction = 'umap', group.by = "Condition",label = FALSE,pt.size = 0.5, label.size = 20) + NoAxes()
dev.off()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#--------------------- II. Summary information on S1234 ----------------------------
data <- read.csv("cell_identities_S1234.csv")
num <- as.numeric(data[,"num_guides"])
num[num > 6] <- ">6"
df<- cbind.data.frame(data,num)
#-------------------------------------------------------------------------------------
df1 <- data.frame(num_guides = as.character(df$num), percent.mito = df$percent_mito, identity=df$Condition)
ggplot(df1,aes(x=factor(num_guides,levels=c(0,1,2,3,4,5,6,">6")), y=percent.mito,fill=factor(num_guides,levels=c(0,1,2,3,4,5,6,">6")))) + 
  #geom_boxplot()
  geom_violin(trim = FALSE)+
  #geom_point(binaxis='y', stackdir='center',dotsize = 0.01,binwidth = 1.5,position = "dodge")+
  #geom_jitter(shape=20, position=position_jitter(0.2))
  geom_boxplot(width=0.1,fill="white")+
  labs(title="",x="No. of sgRNAs per cell", y = "Percent.mito")+
  geom_hline(yintercept=0.1, linetype="dashed", color = "red")+
  scale_fill_brewer(palette="Dark2") + 
  facet_grid(identity ~. )+
  #theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.text=element_text(size=32, colour = "black"),
        axis.title=element_text(size=34, colour = "black"),#face="bold")
        plot.title = element_text(face="bold",size=30,hjust=0.5),
        strip.text.y  = element_text(size = 25)
        ##strip.text is for facet
        )#+coord_fixed(ratio = 2)
ggsave("Percent.mito_S1234.pdf",width = 9.8, height = 11.2)
#-------------------------------------------------------------------------------------
df2 <- data.frame(num=as.character(df$num),nFeature_RNA=df$nFeature_RNA,identity=df$Condition)
ggplot(df2,aes(x=factor(num,levels=c(0,1,2,3,4,5,6,">6")), y=nFeature_RNA,fill=factor(num,levels=c(0,1,2,3,4,5,6,">6")))) + 
  #geom_boxplot()
  geom_violin(trim = FALSE)+
  #geom_point(binaxis='y', stackdir='center',dotsize = 0.01,binwidth = 1.5,position = "dodge")+
  #geom_jitter(shape=20, position=position_jitter(0.2))
  geom_boxplot(width=0.1,fill="white")+
  labs(title="",x="No. of sgRNAs per cell", y = "nFeature_RNA")+
  geom_hline(yintercept=3000, linetype="dashed", color = "black")+
  scale_fill_brewer(palette="Dark2") + 
  facet_grid(identity ~. )+
  #theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.text.y  =element_text(size=30, colour = "black"),
        axis.text.x =element_text(size=32, colour = "black"),
        axis.title=element_text(size=34, colour = "black"),#face="bold")
        plot.title = element_text(face="bold",size=30,hjust=0.5),
        strip.text.y  = element_text(size = 25)
        ##strip.text is for facet
  )
#+scale_y_continuous(limits = c(0,6000),breaks = round(seq(0, 5000, by = 2500),0), labels = c("0", "2.5K", "5K"))
ggsave("nFeature_RNA.S1234.pdf",width = 9.8, height = 11.2)
#-------------------------------------------------------------------------------------
df3 <- data.frame(num=as.character(df$num),nCount_RNA=df$nCount_RNA,identity=df$Condition)
ggplot(df3,aes(x=factor(num,levels=c(0,1,2,3,4,5,6,">6")), y=nCount_RNA,fill=factor(num,levels=c(0,1,2,3,4,5,6,">6")))) + 
  #geom_boxplot()
  geom_violin(trim = FALSE)+
  #geom_point(binaxis='y', stackdir='center',dotsize = 0.01,binwidth = 1.5,position = "dodge")+
  #geom_jitter(shape=20, position=position_jitter(0.2))
  geom_boxplot(width=0.1,fill="white")+
  labs(title="",x="No. of sgRNAs per cell", y = "nCount_RNA")+
  geom_hline(yintercept=15000, linetype="dashed", color = "black")+
  scale_fill_brewer(palette="Dark2") + 
  facet_grid(identity ~. )+
  #theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.text.y  =element_text(size=30, colour = "black"),
        axis.text.x =element_text(size=32, colour = "black"),
        axis.title=element_text(size=34, colour = "black"),#face="bold")
        plot.title = element_text(face="bold",size=30,hjust=0.5),
        strip.text.y  = element_text(size = 25)
        ##strip.text is for facet
  )
ggsave("nCount_RNA.S1234.pdf",width = 10, height = 11.2)
#-------------------------------------------------------------------------------------
##Number of sgRNAs per cell
df4 <- df[c("Condition","num")]
d1 <- df4$num[df$Condition=="S1_D0_full"]
d2 <- df4$num[df$Condition=="S2_D6_full"]
d3 <- df4$num[df$Condition=="S3_D6_tgf"]
d4 <- df4$num[df$Condition=="S4_D6_min"]

###Combine the 
d_1234<- as.matrix(t(rbind(table(d1),table(d2),table(d3),table(d4))))
per <- apply(d_1234,2,function(x){x/sum(x)*100})
std <- apply(per,1,sd)
avg <- apply(per,1,mean)
d_Summary <- data.frame(rownames(per),per,as.numeric(avg),as.numeric(std))
colnames(d_Summary) <- c("num","S1","S2","S3","S4","avg","std")
# data <- setNames(melt(per), c('num', 'samples', 'values'))
# std
#data_stat <- summarySE(data, measurevar="values", groupvars=c("num","samples"))
ggplot(d_Summary,aes(x=factor(num,levels=c(0,1,2,3,4,5,6,">6")), y=avg,fill=factor(num,levels=c(0,1,2,3,4,5,6,">6")))) +
      geom_bar(position=position_dodge(), stat="identity",width = 0.8)+
      geom_errorbar(aes(ymin=avg-std, ymax=avg+std),
                    width=.2,                    # Width of the error bars
                    position=position_dodge(.9))+
  scale_fill_brewer(palette="Dark2") + 
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.text.y  =element_text(size=21, colour = "black"),
        axis.text.x =element_text(size=22, colour = "black"),
        axis.title.x=element_text(size=24, colour = "black"),#face="bold")
        axis.title.y=element_text(size=24, colour = "black", vjust = 2),
        plot.title = element_text(face="bold",size=30,hjust=0.5)
  )+
  scale_y_continuous(limits = c(0,60), breaks = seq(0,60,10))+
  labs(title="",x="No. of sgRNAs per cell", y = "Percentage (%)")
ggsave("Number_of_sgs_per_cell.S1234.pdf",width = 6.6, height = 4.6)

