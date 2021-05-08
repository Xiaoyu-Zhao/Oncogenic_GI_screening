##Xiaoyu Zhao  20200908
library(Seurat)
library(Matrix)
library(ggplot2)
setwd("~/Downloads/5_Combinatorial_CROP_seq/CROP_R/")
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#--------------------- I. Isolate Gene Expression Matrix-----------------------------
#----------------(Whole matrix data = Gene Expression + CRISPR)----------------------
matrix_dir.S4 ="outs/S4_min_D6/S4_filtered_feature_bc_matrix/"
barcode.path.S4 <- paste0(matrix_dir.S4, "barcodes.tsv.gz")
features.path.S4 <- paste0(matrix_dir.S4, "features.tsv.gz")
matrix.path.S4 <- paste0(matrix_dir.S4, "matrix.mtx.gz")
mat.S4 <- readMM(file = matrix.path.S4)

feature.names.S4 = read.delim(features.path.S4, 
                              header = FALSE,
                              stringsAsFactors = FALSE)
barcode.names.S4 = read.delim(barcode.path.S4, 
                              header = FALSE,
                              stringsAsFactors = FALSE)

#Feature names for Gene expression and CRISPR
GEX.names.S4 <- feature.names.S4[grepl("Gene Expression",feature.names.S4[,3]),]
nrow(GEX.names.S4)
CRISPR.names.S4 <- feature.names.S4[grepl("CRISPR",feature.names.S4[,3]),]
nrow(CRISPR.names.S4)

#Matrices for Gene expression and CRISPR
mat.GEX.S4 <- mat.S4[grepl("Gene Expression",feature.names.S4[,3]),]
dim(mat.GEX.S4)
rownames(mat.GEX.S4) <- GEX.names.S4$V1
colnames(mat.GEX.S4) <- barcode.names.S4$V1
mat.GEX.S4[1:3,1:3]

mat.CRISPR.S4 <- mat.S4[grepl("CRISPR",feature.names.S4[,3]),]
dim(mat.CRISPR.S4)
rownames(mat.CRISPR.S4) <- paste0("CRIPSR_",CRISPR.names.S4$V1)
colnames(mat.CRISPR.S4) <- barcode.names.S4$V1
mat.CRISPR.S4[1:3,1:3]
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#-------------------- II. Get CRISPR information from Cellranger--------------------
ident.S4<- as.matrix(read.csv("outs/S4_min_D6/S4_protospacer_calls_per_cell.csv",row.names = 1))
ident.S4 <- ident.S4[,1:3]
dim(ident.S4)  #sgRNA/sgRNAs were dected in 7699 out of 8364 cells
ident.S4[1:3,]

num.S4 <-sort(unique(ident.S4[,1]))
num.S4
table(ident.S4[,1])
# 1    2    3    4    5    6   12 
#2981 4176  348  160   24    9    1 

#---------------II.a Transformations of sgRNA/sgRNA names include:-----------------
sg.comb.S4 <- ident.S4[,2]
# @1. "CTRL000" --> "CTRL"
temp1.S4 <- gsub("CTRL000*","CTRL",sg.comb.S4)
# @2. "|" --> "-"
temp2.S4 <- gsub("[|]","-",temp1.S4)    # Use [|] or "\\|" for "|" 

#-------------------II.b Merge sgRNAs targeting the same genes:---------------------
# @1. Remove no. in sgRNAs, example: "NF2_2" --> "NF2"
temp3.S4 <- gsub("_[1-9]","",temp2.S4)
length(unique(temp3.S4))  #397 unique sgRNA combs (gene level)
sort(unique(temp3.S4))
table(temp3.S4)

#-------------------II.c single sgRNAs  ---> SKO dual sgRNAs: ---------------------
#example: "CBFB_2" --> "CBFB_2-CTRL1"
temp4.S4 = temp3.S4
m = temp3.S4[!grepl("[-]", temp3.S4)]
sort(unique(m)) 
table(m) 
#There are 14 different single sgRNAs in S4
#CASP8(175), CBFB(202), CDH1(155), CTRL1(40), CTRL2(9), hRosa26(10), NF1(313), 
#NF2(150),PTEN(631), RB1(223), SMAD4(677), TBX3(56), TP53(321), USP9X_3(19)
temp4.S4[!grepl("[-]", temp4.S4)] = paste0(m,"-CTRL1")
length(unique(temp4.S4))   #385 unique sgRNA combs (gene level), combing different sgRNAs
sort(unique(temp4.S4))
table(temp4.S4)

#-----------------------II.d Combine all the info above----------------------------
ident2.S4 <- cbind(ident.S4, temp2.S4, temp3.S4, temp4.S4 )
dim(ident2.S4)
#Make the matrix for gex@meta.data
sgrna.ident.S4 <- matrix("Unassigned", nrow = ncol(mat.GEX.S4), 6)
dim(sgrna.ident.S4)
rownames(sgrna.ident.S4) <- barcode.names.S4$V1
overlap.S4 <- intersect(rownames(ident2.S4),rownames(sgrna.ident.S4))
length(overlap.S4)
sgrna.ident.S4[overlap.S4,1:6] <- ident2.S4[overlap.S4,1:6]
colnames(sgrna.ident.S4) <- c("num","ori_guides", "umi","sgRNA_combs", "targets","targets_modi")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#------------------------- III. Create Seurat object--------------------------------
#gex.S4 <- CreateSeuratObject(counts = mat.GEX.S4, project = "S4_D6_1per", min.cells = 3, min.features = 200)
gex.S4 <- CreateSeuratObject(counts = mat.GEX.S4, project = "S4_D6_min", min.cells = 3)
dim(gex.S4)
comb.S4 = sgrna.ident.S4[rownames(gex.S4@meta.data),]
all(rownames(sgrna.ident.S4) == rownames(gex.S4@meta.data))
dim(comb.S4)
comb.S4[comb.S4[,1]=="Unassigned",1] <- 0
comb.S4 <- cbind.data.frame(as.numeric(as.character(comb.S4[,1])),comb.S4[,2:6])

#Add CRISPR info to gex@metadata
gex.S4 <- AddMetaData(gex.S4, metadata = data.frame(comb.S4))
colnames(gex.S4@meta.data) <- c("Condition", "nCount_RNA", "nFeature_RNA", "num_guides","guide_identity","guide_UMI_count", "perturbed_1", "perturbed_2","perturbed_3")
gex.S4@meta.data[1:5,]
#Determine how many genes were targeted by sgRNAs in each cell
num_CTRL <- apply(gex.S4@meta.data, 1, function(x) sum(unlist(gregexpr("CTRL", x[5])) != -1))
gex.S4[["num_targets"]] = gex.S4@meta.data["num_guides"] - num_CTRL
#Determine each cell as SKO, DKO, MKO or Broken
gex.S4[["knock_outs"]] <- rep("MKO", nrow(gex.S4@meta.data))
control <- c("CTRL2-CTRL1", "CTRL1-CTRL1","hRosa26-CTRL1")
gex.S4@meta.data[gex.S4@meta.data$perturbed_3 == "Unassigned", "knock_outs"] = 'Broken'
gex.S4@meta.data[gex.S4@meta.data$num_targets == 1, "knock_outs"] = 'SKO'
gex.S4@meta.data[gex.S4@meta.data$num_targets == 2, "knock_outs"] = 'DKO'
gex.S4@meta.data[gex.S4@meta.data$perturbed_3 %in% control, "knock_outs"] = 'control'

#Reorder sgRNAs in alphabetical order
reorder = function(char){
  if(grepl("CTRL", char)){
    return(char)
  }
  else{
    return(paste(sort(strsplit(char,split = "-")[[1]]), collapse = "-"))
  }
}
gex.S4[["perturbed_1"]]= sapply(as.character(gex.S4@meta.data$perturbed_1), function(x) reorder(x))
gex.S4[["perturbed_2"]]= sapply(as.character(gex.S4@meta.data$perturbed_2), function(x) reorder(x))
gex.S4[["perturbed_3"]]= sapply(as.character(gex.S4@meta.data$perturbed_3), function(x) reorder(x))
gex.S4@meta.data[gex.S4@meta.data$perturbed_1 == "CTRL1-RB1_3","perturbed_1"] = "RB1_3-CTRL1"
gex.S4@meta.data[gex.S4@meta.data$perturbed_2 == "CTRL1-RB1","perturbed_2"] = "RB1-CTRL1"
gex.S4@meta.data[gex.S4@meta.data$perturbed_3 == "CTRL1-RB1","perturbed_3"] = "RB1-CTRL1"
#Final info of perturbations in gex[["perturbed"]]
perturbed  = as.vector(gex.S4@meta.data[["perturbed_3"]])
perturbed[perturbed %in% control] = 'CTRL-CTRL'
gex.S4[["perturbed"]] = perturbed

#Add mito.features info to gex@metadata
mito.features <- GEX.names.S4$V1[grep(pattern = "^MT-", GEX.names.S4$V2)]
percent.mito <- Matrix::colSums(x = GetAssayData(object = gex.S4, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = gex.S4, slot = 'counts'))
gex.S4[["percent_mito"]] <- percent.mito

#Over!
saveRDS(gex.S4, file = "1_Seurat/objects/gex.S4.rds")
write.csv(as.matrix(gex.S4@meta.data),file="1_Seurat/cell_identities_S4.csv")

plot1 <- FeatureScatter(gex.S4, feature1 = "nCount_RNA", feature2 = "percent_mito")+ NoLegend()
plot2 <- FeatureScatter(gex.S4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ NoLegend()
plot1 + plot2
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#--------------------- IV. Summary of perturbations in S4-----------------------------
#All perturbations 
length(unique(gex.S4@meta.data$perturbed_3))
meta = as.data.frame(gex.S4@meta.data)
levels(meta$perturbed_3)
table(meta$perturbed_3)
meta.df = as.data.frame(table(meta$perturbed_3))
meta.df$num_targets = meta[match(meta.df$Var1, meta$perturbed_3), 'num_targets']
num_targets_Category = meta.df$num_targets
num_targets_Category[num_targets_Category > 6] = '> 6'
meta.df$num_targets_Category = num_targets_Category
meta.df$knock_outs = meta[match(meta.df$Var1, meta$perturbed_2), 'knock_outs']
colnames(meta.df)[1] = "perturbed"
dim(meta.df)
write.csv(meta.df, "1_Seurat/S4_all_385_perturbations_S4.csv")

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
write.csv(meta_sel.df, "1_Seurat/S4_all_68_SKO_DKO_perturbations.csv")






