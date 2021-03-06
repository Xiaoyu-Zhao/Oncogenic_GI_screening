##Xiaoyu Zhao  
library(Seurat)
library(Matrix)
library(ggplot2)
library(data.table) 
setwd("~/Downloads/5_Combinatorial_CROP_seq/CROP_R/")
#======================I: Get data and sample info for GSEA==========================================
#---------------------Control cells in Condition B Vs Control cells Condition A--------------------------
#-------------------------------@@@@@@@@@@@@@@@@@@@@@@@@@@@@---------------------------------------------
#Import Seurat objects
gex.S1 <- readRDS( file = "1_Seurat/objects/gex.S1.rds") #dim(gex.S1) = (16735, 5864)
gex.S1 <- subset(gex.S1, subset = knock_outs %in% c("control", "SKO", "DKO") & num_guides %in% c(1,2) & nFeature_RNA > 800 & nFeature_RNA < 5000 & percent_mito < 0.1 )
dim(gex.S1)

gex.S2 <- readRDS( file = "1_Seurat/objects/gex.S2.rds") #dim(gex.S2) = (16575, 7844)
gex.S2 <- subset(gex.S2, subset = knock_outs %in% c("control", "SKO", "DKO") & num_guides %in% c(1,2) & nFeature_RNA > 800 & nFeature_RNA < 5000 & percent_mito < 0.1 )
dim(gex.S2)

gex.S3 <- readRDS( file = "1_Seurat/objects/gex.S3.rds") #dim(gex.S3) = (16975, 6511)
gex.S3 <- subset(gex.S3, subset = knock_outs %in% c("control", "SKO", "DKO") & num_guides %in% c(1,2) & nFeature_RNA > 800 & nFeature_RNA < 5000 & percent_mito < 0.1 )
dim(gex.S3)

gex.S4 <- readRDS( file = "1_Seurat/objects/gex.S4.rds") #dim(gex.S4) = (16915, 8364)
gex.S4 <- subset(gex.S4, subset = knock_outs %in% c("control", "SKO", "DKO") & num_guides %in% c(1,2) & nFeature_RNA > 800 & nFeature_RNA < 5000 & percent_mito < 0.1 )
dim(gex.S4)

#Function "get_data_for_GSEA" is to to get data and sample info of Control cells 
setwd("~/Downloads/5_Combinatorial_CROP_seq/CROP_R/2_Diffexpr_by_conditions/")
get_data_for_GSEA = function(gexB = gex.S4, gexA = gex.S2, sampleB = "S4", sampleA = "S2"){
  #Output expression data for GSEA
  symb = read.csv("symbol_entrezid_CROPseq.csv")
  file_name = paste0(sampleB, "_vs_", sampleA)
  path = paste0(getwd(),"/",file_name,"/")
  dir.create(path)
  gexB <- subset(gexB, perturbed == "CTRL-CTRL")
  gexA <- subset(gexA, perturbed == "CTRL-CTRL")
  d1 <- as.data.frame(as.matrix(GetAssayData(gexB, slot = "counts")))
  d2 <- as.data.frame(as.matrix(GetAssayData(gexA, slot = "counts")))
  commS = intersect(rownames(d1),rownames(d2))
  data <- data.frame(d1[commS,], d2[commS,])
  names <- symb[match(rownames(data),symb$V1),]$V2
  df <- data.frame(names,rownames(d1[commS,]),d1[commS,],d2[commS,])
  colnames(df)[1:2] <- c("NAME","DESCRIPTION")
  fwrite(df,row.names = FALSE, file = paste0(path, file_name, ".txt"), sep="\t")
  #Output sample information for GSEA
  n1 <- ncol(d1)
  n2 <- ncol(d2)
  nn <- sum(n1,n2)
  info1 <- NULL
  info2 <- NULL
  sink(paste0(path, "sampleinfo_", file_name, ".cls"))
  cat(paste(nn,2,1,sep="\t"))
  cat("\n")
  cat(paste("#",sampleB,sampleA,sep="\t"))
  cat("\n")
  for(i in 1:n1){
    info1 <- paste(sampleB,info1,sep="\t")
  }
  cat(info1)
  for(j in 1:n2){
    info2 <- paste(sampleA,info2,sep="\t")
  }
  cat(info2)
  sink()
  
}

#Exexute function "get_data_for_GSEA" 
get_data_for_GSEA(gexB = gex.S3, gexA = gex.S2, sampleB = 'S3', sampleA = "S2")
get_data_for_GSEA(gexB = gex.S4, gexA = gex.S2, sampleB = 'S4', sampleA = "S2")

