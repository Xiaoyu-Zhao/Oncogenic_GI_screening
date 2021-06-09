##Xiaoyu Zhao  20200908
library(Seurat)
library(Matrix)
library(ggplot2)
setwd("~/Downloads/5_Combinatorial_CROP_seq/CROP_R/")
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#--------------------- I. Isolate Gene Expression Matrix-----------------------------
#----------------(Whole matrix data = Gene Expression + CRISPR)----------------------
matrix_dir.S3 ="outs/S3_tgf_D6/S3_filtered_feature_bc_matrix/"
barcode.path.S3 <- paste0(matrix_dir.S3, "barcodes.tsv.gz")
features.path.S3 <- paste0(matrix_dir.S3, "features.tsv.gz")
matrix.path.S3 <- paste0(matrix_dir.S3, "matrix.mtx.gz")
mat.S3 <- readMM(file = matrix.path.S3)

feature.names.S3 = read.delim(features.path.S3, 
                              header = FALSE,
                              stringsAsFactors = FALSE)
barcode.names.S3 = read.delim(barcode.path.S3, 
                              header = FALSE,
                              stringsAsFactors = FALSE)

#Feature names for Gene expression and CRISPR
GEX.names.S3 <- feature.names.S3[grepl("Gene Expression",feature.names.S3[,3]),]
nrow(GEX.names.S3)
CRISPR.names.S3 <- feature.names.S3[grepl("CRISPR",feature.names.S3[,3]),]
nrow(CRISPR.names.S3)

#Matrices for Gene expression and CRISPR
mat.GEX.S3 <- mat.S3[grepl("Gene Expression",feature.names.S3[,3]),]
dim(mat.GEX.S3)
rownames(mat.GEX.S3) <- GEX.names.S3$V1
colnames(mat.GEX.S3) <- barcode.names.S3$V1
mat.GEX.S3[1:3,1:3]

mat.CRISPR.S3 <- mat.S3[grepl("CRISPR",feature.names.S3[,3]),]
dim(mat.CRISPR.S3)
rownames(mat.CRISPR.S3) <- paste0("CRIPSR_",CRISPR.names.S3$V1)
colnames(mat.CRISPR.S3) <- barcode.names.S3$V1
mat.CRISPR.S3[1:3,1:3]
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#-------------------- II. Get CRISPR information from Cellranger--------------------
ident.S3<- as.matrix(read.csv("outs/S3_tgf_D6/S3_protospacer_calls_per_cell.csv",row.names = 1))
ident.S3 <- ident.S3[,1:3]
dim(ident.S3)  #sgRNA/sgRNAs were dectedd in 5937 out of 6511 cells
ident.S3[1:3,]

num.S3 <-sort(unique(ident.S3[,1]))
num.S3
table(ident.S3[,1])
# 1    2    3    4    5    6 
#1922 3755  176   77    6    1 

#---------------II.a Transformations of sgRNA/sgRNA names include:-----------------
sg.comb.S3 <- ident.S3[,2]
# @1. "CTRL000" --> "CTRL"
temp1.S3 <- gsub("CTRL000*","CTRL",sg.comb.S3)
# @2. "|" --> "-"
temp2.S3 <- gsub("[|]","-",temp1.S3)    # Use [|] or "\\|" for "|" 

#-------------------II.b Merge sgRNAs targeting the same genes:---------------------
# @1. Remove no. in sgRNAs, example: "NF2_2" --> "NF2"
temp3.S3 <- gsub("_[1-9]","",temp2.S3)
length(unique(temp3.S3))  #227 unique sgRNA combs (gene level)
sort(unique(temp3.S3))
table(temp3.S3)

#-------------------II.c single sgRNAs  ---> SKO dual sgRNAs: ---------------------
#example: "CBFB_2" --> "CBFB_2-CTRL1"
temp4.S3 = temp3.S3
m = temp3.S3[!grepl("[-]", temp3.S3)]
sort(unique(m)) 
table(m) 
#There are 14 different single sgRNAs in S3
#CASP8(87), CBFB(86), CDH1(69), CTRL1(23), CTRL2(4), hRosa26(3), NF1(123), NF2(130),
#PTEN(174), RB1(103), SMAD4(956), TBX3(20), TP53(132), USP9X(12)
temp4.S3[!grepl("[-]", temp4.S3)] = paste0(m,"-CTRL1")
length(unique(temp4.S3))   #215 unique sgRNA combs (gene level), combing different sgRNAs
sort(unique(temp4.S3))
table(temp4.S3)

#-----------------------II.c Combine all the info above----------------------------
ident2.S3 <- cbind(ident.S3, temp2.S3, temp3.S3, temp4.S3 )
dim(ident2.S3)
#Make the matrix for gex@meta.data
sgrna.ident.S3 <- matrix("Unassigned", nrow = ncol(mat.GEX.S3), 6)
dim(sgrna.ident.S3)
rownames(sgrna.ident.S3) <- barcode.names.S3$V1
overlap.S3 <- intersect(rownames(ident2.S3),rownames(sgrna.ident.S3))
length(overlap.S3)
sgrna.ident.S3[overlap.S3,1:6] <- ident2.S3[overlap.S3,1:6]
colnames(sgrna.ident.S3) <- c("num","ori_guides", "umi","sgRNA_combs", "targets", "targets_modi")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#------------------------- III. Create Seurat object--------------------------------
#gex.S3 <- CreateSeuratObject(counts = mat.GEX.S3, project = "S3_D6_1per", min.cells = 3, min.features = 200)
gex.S3 <- CreateSeuratObject(counts = mat.GEX.S3, project = "S3_D6_tgf", min.cells = 3)
dim(gex.S3)
comb.S3 = sgrna.ident.S3[rownames(gex.S3@meta.data),]
all(rownames(sgrna.ident.S3) == rownames(gex.S3@meta.data))
dim(comb.S3)
comb.S3[comb.S3[,1]=="Unassigned",1] <- 0
comb.S3 <- cbind.data.frame(as.numeric(as.character(comb.S3[,1])),comb.S3[,2:6])

#Add CRISPR info to gex@metadata
gex.S3 <- AddMetaData(gex.S3, metadata = data.frame(comb.S3))
colnames(gex.S3@meta.data) <- c("Condition", "nCount_RNA", "nFeature_RNA", "num_guides","guide_identity","guide_UMI_count", "perturbed_1", "perturbed_2","perturbed_3")
gex.S3@meta.data[1:5,]
#Determine how many genes were targeted by sgRNAs in each cell
num_CTRL <- apply(gex.S3@meta.data, 1, function(x) sum(unlist(gregexpr("CTRL", x[5])) != -1))
gex.S3[["num_targets"]] = gex.S3@meta.data["num_guides"] - num_CTRL
#Determine each cell as SKO, DKO, MKO or Broken
gex.S3[["knock_outs"]] <- rep("MKO", nrow(gex.S3@meta.data))
control <- c("CTRL2-CTRL1", "CTRL1-CTRL1", "hRosa26-CTRL1")
gex.S3@meta.data[gex.S3@meta.data$perturbed_3 == "Unassigned", "knock_outs"] = 'Broken'
gex.S3@meta.data[gex.S3@meta.data$num_targets == 1, "knock_outs"] = 'SKO'
gex.S3@meta.data[gex.S3@meta.data$num_targets == 2, "knock_outs"] = 'DKO'
gex.S3@meta.data[gex.S3@meta.data$perturbed_3 %in% control, "knock_outs"] = 'control'

#Reorder sgRNAs in alphabetical order
reorder = function(char){
  if(grepl("CTRL", char)){
    return(char)
  }
  else{
    return(paste(sort(strsplit(char,split = "-")[[1]]), collapse = "-"))
  }
}

gex.S3[["perturbed_1"]]= sapply(as.character(gex.S3@meta.data$perturbed_1), function(x) reorder(x))
gex.S3[["perturbed_2"]]= sapply(as.character(gex.S3@meta.data$perturbed_2), function(x) reorder(x))
gex.S3[["perturbed_3"]]= sapply(as.character(gex.S3@meta.data$perturbed_3), function(x) reorder(x))
gex.S3@meta.data[gex.S3@meta.data$perturbed_1 == "CTRL1-RB1_3","perturbed_1"] = "RB1_3-CTRL1"
gex.S3@meta.data[gex.S3@meta.data$perturbed_2 == "CTRL1-RB1","perturbed_2"] = "RB1-CTRL1"
gex.S3@meta.data[gex.S3@meta.data$perturbed_3 == "CTRL1-RB1","perturbed_3"] = "RB1-CTRL1"
#Final info of perturbations in gex[["perturbed"]]
perturbed  = as.vector(gex.S3@meta.data[["perturbed_3"]])
perturbed[perturbed %in% control] = 'CTRL-CTRL'
gex.S3[["perturbed"]] = perturbed

#Add mito.features info to gex@metadata
mito.features <- GEX.names.S3$V1[grep(pattern = "^MT-", GEX.names.S3$V2)]
percent.mito <- Matrix::colSums(x = GetAssayData(object = gex.S3, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = gex.S3, slot = 'counts'))
gex.S3[["percent_mito"]] <- percent.mito

#Over!
saveRDS(gex.S3, file = "1_Seurat/objects/gex.S3.rds")
write.csv(as.matrix(gex.S3@meta.data),file="1_Seurat/cell_identities_S3.csv")

plot1 <- FeatureScatter(gex.S3, feature1 = "nCount_RNA", feature2 = "percent_mito")+ NoLegend()
plot2 <- FeatureScatter(gex.S3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ NoLegend()
plot1 + plot2
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#--------------------- IV. Summary of perturbations in S3-----------------------------
#All perturbations 
length(unique(gex.S3@meta.data$perturbed_3))
meta = as.data.frame(gex.S3@meta.data)
levels(meta$perturbed_3)
table(meta$perturbed_3)
meta.df = as.data.frame(table(meta$perturbed_3))
meta.df$num_targets = meta[match(meta.df$Var1, meta$perturbed_3), 'num_targets']
num_targets_Category = meta.df$num_targets
num_targets_Category[num_targets_Category > 6] = '> 6'
meta.df$num_targets_Category = num_targets_Category
meta.df$knock_outs = meta[match(meta.df$Var1, meta$perturbed_3), 'knock_outs']
colnames(meta.df)[1] = "perturbed"
dim(meta.df)
write.csv(meta.df, "1_Seurat/S3_all_215_perturbations.csv")

#All SKO, DKO and control perturbations
meta_sel = subset(meta, meta$num_guides == 1 | meta$num_guides == 2)
#a = droplevels(meta_sel$perturbed)
#b = factor(meta_sel$perturbed)
#all(a == b)
table(meta_sel$perturbed)
meta_sel.df = as.data.frame(table(meta_sel$perturbed))
meta_sel.df$num_targets = meta_sel[match(meta_sel.df$Var1, meta_sel$perturbed), 'num_targets']
meta_sel.df$knock_outs = meta_sel[match(meta_sel.df$Var1, meta_sel$perturbed), 'knock_outs']
colnames(meta_sel.df)[1] = "perturbed"
dim(meta_sel.df)
write.csv(meta_sel.df, "1_Seurat/S3_all_66_SKO_DKO_perturbations.csv")






