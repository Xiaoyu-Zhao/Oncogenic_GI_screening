##Xiaoyu Zhao 
library(Seurat)
library(Matrix)
library(ggplot2)
setwd("~/Downloads/5_Combinatorial_CROP_seq/CROP_R/")
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#--------------------- I. Isolate Gene Expression Matrix-----------------------------
#----------------(Whole matrix data = Gene Expression + CRISPR)----------------------
matrix_dir.S1 ="outs/S1_full_D0/S1_filtered_feature_bc_matrix/"
barcode.path.S1 <- paste0(matrix_dir.S1, "barcodes.tsv.gz")
features.path.S1 <- paste0(matrix_dir.S1, "features.tsv.gz")
matrix.path.S1 <- paste0(matrix_dir.S1, "matrix.mtx.gz")
mat.S1 <- readMM(file = matrix.path.S1)

feature.names.S1 = read.delim(features.path.S1, 
                              header = FALSE,
                              stringsAsFactors = FALSE)
barcode.names.S1 = read.delim(barcode.path.S1, 
                              header = FALSE,
                              stringsAsFactors = FALSE)

#Feature names for Gene expression and CRISPR
GEX.names.S1 <- feature.names.S1[grepl("Gene Expression",feature.names.S1[,3]),]
nrow(GEX.names.S1)
CRISPR.names.S1 <- feature.names.S1[grepl("CRISPR",feature.names.S1[,3]),]
nrow(CRISPR.names.S1)

#Matrices for Gene expression and CRISPR
mat.GEX.S1 <- mat.S1[grepl("Gene Expression",feature.names.S1[,3]),]
dim(mat.GEX.S1)
rownames(mat.GEX.S1) <- GEX.names.S1$V1
colnames(mat.GEX.S1) <- barcode.names.S1$V1
mat.GEX.S1[1:3,1:3]

mat.CRISPR.S1 <- mat.S1[grepl("CRISPR",feature.names.S1[,3]),]
dim(mat.CRISPR.S1)
rownames(mat.CRISPR.S1) <- paste0("CRIPSR_",CRISPR.names.S1$V1)
colnames(mat.CRISPR.S1) <- barcode.names.S1$V1
mat.CRISPR.S1[1:3,1:3]
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#-------------------- II. Get CRISPR information from Cellranger--------------------
ident.S1<- as.matrix(read.csv("outs/S1_full_D0/S1_protospacer_calls_per_cell.csv",row.names = 1))
ident.S1 <- ident.S1[,1:3]
dim(ident.S1)  #sgRNA/sgRNAs were dectedd in 5296 out of 5864 cells
ident.S1[1:3,]

num.S1 <-sort(unique(ident.S1[,1]))
num.S1
table(ident.S1[,1])
#  1    2    3    4    5    6    7 
#1854 3183  178   73    6    1    1

#---------------II.a Transformations of sgRNA/sgRNA names include:-----------------
sg.comb.S1 <- ident.S1[,2]
# @1. "CTRL000" --> "CTRL"
temp1.S1 <- gsub("CTRL000*","CTRL",sg.comb.S1)
# @2. "|" --> "-"
temp2.S1 <- gsub("[|]","-",temp1.S1)    # Use [|] or "\\|" for "|" 

#-------------------II.b Merge sgRNAs targeting the same genes:---------------------
# @1. Remove no. in sgRNAs, example: "NF2_2" --> "NF2"
temp3.S1 <- gsub("_[1-9]","",temp2.S1)
length(unique(temp3.S1))  #263 unique sgRNA combs (gene level)
sort(unique(temp3.S1))
table(temp3.S1)

#-------------------II.c single sgRNAs  ---> SKO dual sgRNAs: ---------------------
#example: "CBFB_2" --> "CBFB_2-CTRL1"
temp4.S1 = temp3.S1
m = temp3.S1[!grepl("[-]", temp3.S1)]
sort(unique(m)) 
table(m) 
#There are 14 different single sgRNAs in S1
#CASP8(138), CBFB(173), CDH1(113), CTRL1(20), CTRL2(11), hRosa26(9), NF1(242), NF2(93), 
#PTEN(277), RB1(143), SMAD4(388), TBX3(23), TP53(216), USP9X(8)
temp4.S1[!grepl("[-]", temp4.S1)] = paste0(m,"-CTRL1")
length(unique(temp4.S1))   #251 unique sgRNA combs (gene level), combining different sgRNAs
sort(unique(temp4.S1))
table(temp4.S1)

#-----------------------II.d Combine all the info above----------------------------
ident2.S1 <- cbind(ident.S1, temp2.S1, temp3.S1, temp4.S1 )
dim(ident2.S1)
#Make the matrix for gex@meta.data
sgrna.ident.S1 <- matrix("Unassigned", nrow = ncol(mat.GEX.S1), 6)
dim(sgrna.ident.S1)
rownames(sgrna.ident.S1) <- barcode.names.S1$V1
overlap.S1 <- intersect(rownames(ident2.S1),rownames(sgrna.ident.S1))
length(overlap.S1)
sgrna.ident.S1[overlap.S1,1:6] <- ident2.S1[overlap.S1,1:6]
colnames(sgrna.ident.S1) <- c("num","ori_guides", "umi","sgRNA_combs", "targets","targets_modi")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#------------------------- III. Create Seurat object--------------------------------
#gex.S1 <- CreateSeuratObject(counts = mat.GEX.S1, project = "S1_D6_1per", min.cells = 3, min.features = 200)
gex.S1 <- CreateSeuratObject(counts = mat.GEX.S1, project = "S1_D0_full", min.cells = 3)
dim(gex.S1)
comb.S1 = sgrna.ident.S1[rownames(gex.S1@meta.data),]
all(rownames(sgrna.ident.S1) == rownames(gex.S1@meta.data))
dim(comb.S1)
comb.S1[comb.S1[,1]=="Unassigned",1] <- 0
comb.S1 <- cbind.data.frame(as.numeric(as.character(comb.S1[,1])),comb.S1[,2:6])

#Add CRISPR info to gex@metadata
gex.S1 <- AddMetaData(gex.S1, metadata = data.frame(comb.S1))
colnames(gex.S1@meta.data) <- c("Condition", "nCount_RNA", "nFeature_RNA", "num_guides","guide_identity","guide_UMI_count", "perturbed_1", "perturbed_2","perturbed_3")
gex.S1@meta.data[1:5,]
#Determine how many genes were targeted by sgRNAs in each cell
num_CTRL <- apply(gex.S1@meta.data, 1, function(x) sum(unlist(gregexpr("CTRL", x[5])) != -1))
gex.S1[["num_targets"]] = gex.S1@meta.data["num_guides"] - num_CTRL
#Determine each cell as SKO, DKO, MKO or Broken
gex.S1[["knock_outs"]] <- rep("MKO", nrow(gex.S1@meta.data))
control <- c("CTRL2-CTRL1", "CTRL1-CTRL1", "hRosa26-CTRL1")
gex.S1@meta.data[gex.S1@meta.data$perturbed_3 == "Unassigned", "knock_outs"] = 'Broken'
gex.S1@meta.data[gex.S1@meta.data$num_targets == 1, "knock_outs"] = 'SKO'
gex.S1@meta.data[gex.S1@meta.data$num_targets == 2, "knock_outs"] = 'DKO'
gex.S1@meta.data[gex.S1@meta.data$perturbed_3 %in% control, "knock_outs"] = 'control'

#Reorder sgRNAs in alphabetical order
reorder = function(char){
  if(grepl("CTRL", char)){
    return(char)
  }
  else{
    return(paste(sort(strsplit(char,split = "-")[[1]]), collapse = "-"))
  }
}
gex.S1[["perturbed_1"]]= sapply(as.character(gex.S1@meta.data$perturbed_1), function(x) reorder(x))
gex.S1[["perturbed_2"]]= sapply(as.character(gex.S1@meta.data$perturbed_2), function(x) reorder(x))
gex.S1[["perturbed_3"]]= sapply(as.character(gex.S1@meta.data$perturbed_3), function(x) reorder(x))
gex.S1@meta.data[gex.S1@meta.data$perturbed_1 == "CTRL1-RB1_3","perturbed_1"] = "RB1_3-CTRL1"
gex.S1@meta.data[gex.S1@meta.data$perturbed_2 == "CTRL1-RB1","perturbed_2"] = "RB1-CTRL1"
gex.S1@meta.data[gex.S1@meta.data$perturbed_3 == "CTRL1-RB1","perturbed_3"] = "RB1-CTRL1"
#Final info of perturbations in gex[["perturbed"]]
perturbed  = as.vector(gex.S1@meta.data[["perturbed_3"]])
perturbed[perturbed %in% control] = 'CTRL-CTRL'
gex.S1[["perturbed"]] = perturbed

#Add mito.features info to gex@metadata
mito.features <- GEX.names.S1$V1[grep(pattern = "^MT-", GEX.names.S1$V2)]
percent.mito <- Matrix::colSums(x = GetAssayData(object = gex.S1, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = gex.S1, slot = 'counts'))
gex.S1[["percent_mito"]] <- percent.mito

#Over!
saveRDS(gex.S1, file = "1_Seurat/objects/gex.S1.rds")
write.csv(as.matrix(gex.S1@meta.data),file="1_Seurat/cell_identities_S1.csv")

plot1 <- FeatureScatter(gex.S1, feature1 = "nCount_RNA", feature2 = "percent_mito")+ NoLegend()
plot2 <- FeatureScatter(gex.S1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ NoLegend()
plot1 + plot2
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#--------------------- IV. Summary of perturbations in S1-----------------------------
#All perturbations 
length(unique(gex.S1@meta.data$perturbed_3))
meta = as.data.frame(gex.S1@meta.data)
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
write.csv(meta.df, "1_Seurat/S1_all_251_perturbations.csv")

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
write.csv(meta_sel.df, "1_Seurat/S1_all_67_SKO_DKO_perturbations.csv")






