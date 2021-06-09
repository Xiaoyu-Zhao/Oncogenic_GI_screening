##Xiaoyu Zhao  
library(Seurat)
library(Matrix)
library(ggplot2)
setwd("~/Downloads/5_Combinatorial_CROP_seq/CROP_R/")
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#--------------------- I. Isolate Gene Expression Matrix-----------------------------
#----------------(Whole matrix data = Gene Expression + CRISPR)----------------------
matrix_dir.S2 ="outs/S2_full_D6/S2_filtered_feature_bc_matrix/"
barcode.path.S2 <- paste0(matrix_dir.S2, "barcodes.tsv.gz")
features.path.S2 <- paste0(matrix_dir.S2, "features.tsv.gz")
matrix.path.S2 <- paste0(matrix_dir.S2, "matrix.mtx.gz")
mat.S2 <- readMM(file = matrix.path.S2)

feature.names.S2 = read.delim(features.path.S2, 
                              header = FALSE,
                              stringsAsFactors = FALSE)
barcode.names.S2 = read.delim(barcode.path.S2, 
                              header = FALSE,
                              stringsAsFactors = FALSE)

#Feature names for Gene expression and CRISPR
GEX.names.S2 <- feature.names.S2[grepl("Gene Expression",feature.names.S2[,3]),]
nrow(GEX.names.S2)
CRISPR.names.S2 <- feature.names.S2[grepl("CRISPR",feature.names.S2[,3]),]
nrow(CRISPR.names.S2)

#Matrices for Gene expression and CRISPR
mat.GEX.S2 <- mat.S2[grepl("Gene Expression",feature.names.S2[,3]),]
dim(mat.GEX.S2)
rownames(mat.GEX.S2) <- GEX.names.S2$V1
colnames(mat.GEX.S2) <- barcode.names.S2$V1
mat.GEX.S2[1:3,1:3]

mat.CRISPR.S2 <- mat.S2[grepl("CRISPR",feature.names.S2[,3]),]
dim(mat.CRISPR.S2)
rownames(mat.CRISPR.S2) <- paste0("CRIPSR_",CRISPR.names.S2$V1)
colnames(mat.CRISPR.S2) <- barcode.names.S2$V1
mat.CRISPR.S2[1:3,1:3]
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#-------------------- II. Get CRISPR information from Cellranger--------------------
ident.S2<- as.matrix(read.csv("outs/S2_full_D6/S2_protospacer_calls_per_cell.csv",row.names = 1))
ident.S2 <- ident.S2[,1:3]
dim(ident.S2)  #sgRNA/sgRNAs were dectedd in 6862 out of 7844 cells
ident.S2[1:3,]

num.S2 <-sort(unique(ident.S2[,1]))
num.S2
table(ident.S2[,1])
#  1    2    3    4    5    6    7    9 
#2395 4092  250  114    6    3    1    1

#---------------II.a Transformations of sgRNA/sgRNA names include:-----------------
sg.comb.S2 <- ident.S2[,2]
# @1. "CTRL000" --> "CTRL"
temp1.S2 <- gsub("CTRL000*","CTRL",sg.comb.S2)
# @2. "|" --> "-"
temp2.S2 <- gsub("[|]","-",temp1.S2)    # Use [|] or "\\|" for "|" 

#-------------------II.b Merge sgRNAs targeting the same genes:---------------------
# @1. Remove no. in sgRNAs, example: "NF2_2" --> "NF2"
temp3.S2 <- gsub("_[1-9]","",temp2.S2)
length(unique(temp3.S2))  #301 unique sgRNA combs (gene level)
sort(unique(temp3.S2))
table(temp3.S2)

#-------------------II.c single sgRNAs  ---> SKO dual sgRNAs: ---------------------
temp4.S2 = temp3.S2
#example: "CBFB_2" --> "CBFB_2-CTRL1"
m = temp3.S2[!grepl("[-]", temp3.S2)]
sort(unique(m)) 
table(m) 
#There are 14 different single sgRNAs in S2
#CASP8(166), CBFB(171), CDH1(124), CTRL1(22), CTRL2(5), hRosa26(8), NF1(214), NF2(103)
#PTEN(314), RB1(147), SMAD4(870), TBX3(32), TP53(211),  USP9X(8)
temp4.S2[!grepl("[-]", temp4.S2)] = paste0(m,"-CTRL1")
length(unique(temp4.S2))   #301 unique sgRNA combs (gene level),combing different sgRNAs
sort(unique(temp4.S2))
table(temp4.S2)

#-----------------------II.d Combine all the info above----------------------------
ident2.S2 <- cbind(ident.S2, temp2.S2, temp3.S2, temp4.S2 )
dim(ident2.S2)
#Make the matrix for gex@meta.data
sgrna.ident.S2 <- matrix("Unassigned", nrow = ncol(mat.GEX.S2), 6)
dim(sgrna.ident.S2)
rownames(sgrna.ident.S2) <- barcode.names.S2$V1
overlap.S2 <- intersect(rownames(ident2.S2),rownames(sgrna.ident.S2))
length(overlap.S2)
sgrna.ident.S2[overlap.S2,1:6] <- ident2.S2[overlap.S2,1:6]
colnames(sgrna.ident.S2) <- c("num","ori_guides", "umi","sgRNA_combs", "targets","targets_modi")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#------------------------- III. Create Seurat object--------------------------------
#gex.S2 <- CreateSeuratObject(counts = mat.GEX.S2, project = "S2_D6_1per", min.cells = 3, min.features = 200)
gex.S2 <- CreateSeuratObject(counts = mat.GEX.S2, project = "S2_D6_full", min.cells = 3)
dim(gex.S2)
comb.S2 = sgrna.ident.S2[rownames(gex.S2@meta.data),]
all(rownames(sgrna.ident.S2) == rownames(gex.S2@meta.data))
dim(comb.S2)
comb.S2[comb.S2[,1]=="Unassigned",1] <- 0
comb.S2 <- cbind.data.frame(as.numeric(as.character(comb.S2[,1])),comb.S2[,2:6])

#Add CRISPR info to gex@metadata
gex.S2 <- AddMetaData(gex.S2, metadata = data.frame(comb.S2))
colnames(gex.S2@meta.data) <- c("Condition", "nCount_RNA", "nFeature_RNA", "num_guides","guide_identity","guide_UMI_count", "perturbed_1", "perturbed_2","perturbed_3")
gex.S2@meta.data[1:5,]
#Determine how many genes were targeted by sgRNAs in each cell
num_CTRL <- apply(gex.S2@meta.data, 1, function(x) sum(unlist(gregexpr("CTRL", x[5])) != -1))
gex.S2[["num_targets"]] = gex.S2@meta.data["num_guides"] - num_CTRL
#Determine each cell as SKO, DKO, MKO or Broken
gex.S2[["knock_outs"]] <- rep("MKO", nrow(gex.S2@meta.data))
control <- c("CTRL2-CTRL1", "CTRL1-CTRL1", "hRosa26-CTRL1")
gex.S2@meta.data[gex.S2@meta.data$perturbed_3 == "Unassigned", "knock_outs"] = 'Broken'
gex.S2@meta.data[gex.S2@meta.data$num_targets == 1, "knock_outs"] = 'SKO'
gex.S2@meta.data[gex.S2@meta.data$num_targets == 2, "knock_outs"] = 'DKO'
gex.S2@meta.data[gex.S2@meta.data$perturbed_3 %in% control, "knock_outs"] = 'control'

#Reorder sgRNAs in alphabetical order
reorder = function(char){
  if(grepl("CTRL", char)){
    return(char)
  }
  else{
    return(paste(sort(strsplit(char,split = "-")[[1]]), collapse = "-"))
  }
}
gex.S2[["perturbed_1"]]= sapply(as.character(gex.S2@meta.data$perturbed_1), function(x) reorder(x))
gex.S2[["perturbed_2"]]= sapply(as.character(gex.S2@meta.data$perturbed_2), function(x) reorder(x))
gex.S2[["perturbed_3"]]= sapply(as.character(gex.S2@meta.data$perturbed_3), function(x) reorder(x))
gex.S2@meta.data[gex.S2@meta.data$perturbed_1 == "CTRL1-RB1_3","perturbed_1"] = "RB1_3-CTRL1"
gex.S2@meta.data[gex.S2@meta.data$perturbed_2 == "CTRL1-RB1","perturbed_2"] = "RB1-CTRL1"
gex.S2@meta.data[gex.S2@meta.data$perturbed_3 == "CTRL1-RB1","perturbed_3"] = "RB1-CTRL1"
#Final info of perturbations in gex[["perturbed"]]
perturbed  = as.vector(gex.S2@meta.data[["perturbed_3"]])
perturbed[perturbed %in% control] = 'CTRL-CTRL'
gex.S2[["perturbed"]] = perturbed

#Add mito.features info to gex@metadata
mito.features <- GEX.names.S2$V1[grep(pattern = "^MT-", GEX.names.S2$V2)]
percent.mito <- Matrix::colSums(x = GetAssayData(object = gex.S2, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = gex.S2, slot = 'counts'))
gex.S2[["percent_mito"]] <- percent.mito

#Over!
saveRDS(gex.S2, file = "1_Seurat/objects/gex.S2.rds")
write.csv(as.matrix(gex.S2@meta.data),file="1_Seurat/cell_identities_S2.csv")

plot1 <- FeatureScatter(gex.S2, feature1 = "nCount_RNA", feature2 = "percent_mito")+ NoLegend()
plot2 <- FeatureScatter(gex.S2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ NoLegend()
plot1 + plot2
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#--------------------- IV. Summary of perturbations in S2-----------------------------
#All perturbations 
length(unique(gex.S2@meta.data$perturbed_3))
meta = as.data.frame(gex.S2@meta.data)
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
write.csv(meta.df, "1_Seurat/S2_all_301_perturbations.csv")

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
write.csv(meta_sel.df, "1_Seurat/S2_all_67_SKO_DKO_perturbations.csv")






