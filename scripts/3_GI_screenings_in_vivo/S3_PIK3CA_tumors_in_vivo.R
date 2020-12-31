#Xiaoyu Zhao 
RNGkind(sample.kind = "Rounding")
# RNGkind is to get the results from set.seed() to match for different R versions
library(ggplot2)
library(dplyr)
library(scales)
library(car) 
library(readxl)
library(RColorBrewer)
setwd("~/Downloads/data/")
#=====================Part 1: Calculate LFC based on 70th percentile==========================================
#-------------------------------@@@@@@@@@@@@@@@@@@@@@@@@@@@@-------------------------------------------------
#---------------------------I. Intra tumor normalization: RPM (reads per million) ---------------------------
ct.pik <- data.frame(read_excel("TableS3_Raw_and_processed_data_combinatorial_screening_in_vivo.xlsx", 2))
ref <- data.frame(read_excel("TableS1_Annotatation_of_asssyed_TSGs_or_TSG_combinations.xlsx", 1))

path = paste0(getwd(), "/", "PIK3CA_tumors/")
dir.create(path) 
setwd(path)
#------------------------------------------------------------------------------------------------------------
rownames(ct.pik) = ct.pik$gene_names
#Remove rows with "CTRL0003","PTEN_1","PTEN_2","PTEN_3" sgRNAs
bad_candidate = c("CTRL0003","PTEN_1","PTEN_2","PTEN_3")
filter_flag_1 = sapply(bad_candidate, function(x) {grepl(x, rownames(ct.pik))})
filter_flag_2 = apply(filter_flag_1, 1, sum)
filter_flag = ifelse(filter_flag_2==0, TRUE, FALSE)
counts = ct.pik[filter_flag, 3:ncol(ct.pik) ]
dim(counts)
write.csv(counts, "cell_tumor_counts_PIK3CA.csv")
#------------------------------------------------------------------------------------------------------------
counts = read.csv('cell_tumor_counts_PIK3CA.csv', row.names = 1, header = T)
keep <- rowSums(counts[,1:2]) >= 2
cts = counts[keep,]
RPM = as.matrix(apply(cts[,1:ncol(cts)]+0.1, 2, function(x) {x/sum(x)*1000000})) 

#reference gene_comb names for PTEN-/- tumors
sname <- ref$gene_comb
write.csv(sname, "gene_comb_names_PIK3CA.csv")

#----------------------------II. Calculate RPM at gene level: 70th percentile --------------------------------
##Qct function is to calculate the 70th percentile RPM of each gene pair of each sample
Qct  = function(count, qct = 0.7){
  QctM = NULL
  for(i in 1:length(sname)){
    if(sname[i] != "CTRL-CTRL"){
      g1 <- unlist(strsplit(sname[i], split = "-")[1])[1]
      g2 <- unlist(strsplit(sname[i], split = "-")[1])[2]
      gg <- which(grepl(g1, rownames(count)) & grepl(g2, rownames(count)) == 1)
      temp = apply(count, 2, function(x){quantile(as.numeric(x[gg]),qct)})
    }
    else {
      gg = which(grepl("^CTRL", rownames(count)) & grepl("-CTRL", rownames(count)))
      temp = apply(count, 2, function(x){quantile(as.numeric(x[gg]),qct)})
    }
    QctM = rbind(QctM, temp)
  } 
  rownames(QctM) = sname
  return (QctM)
}
QctM = Qct(RPM, qct = 0.7)
QctM = QctM[rowSums(is.na(QctM))==0,]
dim(QctM)

## Normalize gene-level RPM within each tumor sample
QctMM = apply(QctM, 2, function(x) {x/sum(x)*1000000})
#----------------------------III. Determine gene level effect size (LFC) --------------------------------
#Determine effect size lfc:
#Formula : lfc = log2(mean(Qct_tumor))/mean(Qct_cell))
t1 = apply(QctMM[,1:2], 1, function(x) mean(x))
range(t1)
#sum(t1 == 0)
t2 = apply(QctMM[,3:ncol(QctMM)], 1, function(x) mean(x))
range(t2)
fc = t2/t1
range(fc)
#hist(fc, breaks = 100)
lfc <- log2(fc)
LFC = cbind.data.frame(Qct_cell = t1, Qct_tumor = t2, FC = fc, LFC = lfc)
#----------------------------IV. Determine gene level effect size (LFC) --------------------------------
#Determine significance of effect size (lfc) and confidence interval (CI) using Bootstrap Scott proposed
set.seed(42)
LFC.Boot = NULL
Boot = NULL
nBoot = 1000
for(i in 1:nBoot){
  idx.cell = sample(1:2, size = 2, replace = T)
  idx.tumor = sample(3:ncol(QctMM), size = ncol(QctMM)-2, replace = T) 
  temp = apply(QctMM, 1, function(x){log2(mean(c(x[idx.tumor]))/mean(c(x[idx.cell])))})
  Boot = cbind(Boot, temp)
}

LFC.Boot = cbind(lfc, Boot)
p_value = apply(LFC.Boot, 1, function(x) {if(x[1] > 0 ) {sum(x[2:(nBoot+1)] <= 0)/nBoot} else {sum(x[2:(nBoot+1)] >= 0)/nBoot}})
p_adj = p.adjust(p_value, method = 'BH')

l <- apply(Boot, 1, function(x){quantile(x, 0.025)})
u <- apply(Boot, 1, function(x){quantile(x, 0.975)})
LFC.zscore = scale(LFC$LFC)
res.Boot = cbind.data.frame(QctMM, LFC,  l, u, LFC.zscore, p_value, p_adj)
write.csv(res.Boot, 'res_LFC_PIK3CA.csv')

#----------------------------------@@@@@@@@@@@@@@@@@@@@@@@@@@@@-------------------------------------------------
#========================================Part 2: Outlier Test===================================================
#----------------------------------@@@@@@@@@@@@@@@@@@@@@@@@@@@@-------------------------------------------------
#---------------------------------I. Calculate the mean log2RPM at gene level-----------------------------------
counts = read.csv('cell_tumor_counts_PIK3CA.csv', row.names = 1, header = T)
ref <- read.csv("gene_comb_names_PIK3CA.csv",row.names = 1)
sname <- as.character(ref[,1])
#-------------------------------------------------------------------------------------------------------------
#Calculate log2RPM 
keep <- rowSums(counts[,1:2]) >= 2
cts = counts[keep,]
RPM = as.matrix(apply(cts[,1:ncol(cts)], 2, function(x) {x/sum(x)*1000000})) 
dim(RPM)
log2rpm = as.data.frame(apply(RPM[,1:ncol(RPM)], 2, function(x) {log2(x+1)})) ###0.01 or 0.001
write.csv(log2rpm,"log2RPM_PIK3CA.csv")
#-------------------------------------------------------------------------------------------------------------
#Calculate the mean log2RPM at gene level
c.gg <- c() ##ave of log-relative-freq of each gene pair in cells
t.gg <- c()##ave of log-relative-freq of each gene pair in tumors

for(i in 1:length(sname)){
  g1 <- paste0("^", unlist(strsplit(sname[i], split = "-")[1])[1])
  g2 <- paste0("-", unlist(strsplit(sname[i], split = "-")[1])[2])
  gg <- rownames(log2rpm)[grepl(g1, rownames(log2rpm)) & grepl(g2, rownames(log2rpm))]
  c.gg[i] <- mean(as.numeric(as.matrix(log2rpm[gg,1:2])))
  t.gg[i] <- mean(as.numeric(as.matrix(log2rpm[gg,3:ncol(log2rpm)])))
}
lf <- cbind(c.gg, t.gg)
rownames(lf) <- sname
df <- data.frame(sname = sname, lf)
fit = lm(df$t.gg ~ df$c.gg)
#plot(fit)
#------------------------------------II. Outlier test on Avg.log2RPM--------------------------------------------
##"outlierTest" reports the Bonferroni p-values for testing each observation in turn to be a mean-shift outlier, 
#based Studentized residuals in linear (t-tests) and generalized linear models (normal tests).
tt = outlierTest(fit, cutoff = Inf, n.max = nrow(df)) 
d1 = df[as.numeric(names(tt[[1]])), ]
d2 = do.call(cbind.data.frame, tt)
p.adj.FDR = p.adjust(d2$p, method = "BH")
m = cbind(d1, d2[,1:3], p.adj.FDR)
write.csv(m, "Outlier_gene_level_PIK3CA.csv", row.names = F)

#----------------------------------@@@@@@@@@@@@@@@@@@@@@@@@@@@@-------------------------------------------------
#========================Part 3: How many dual-sgRNAs were enriched for each gene pair==========================
#----------------------------------@@@@@@@@@@@@@@@@@@@@@@@@@@@@-------------------------------------------------
#-----------------------I. Calculate the mean log2RPM from all cell/tumor samples-----------------------------
log2rpm = read.csv("log2RPM_PIK3CA.csv", row.names = 1)
mean1 = apply(log2rpm[,1:2], 1, function(x) mean(x))
mean2 = apply(log2rpm[, 3:ncol(log2rpm)], 1, function(x) mean(x))
dat0 = data.frame(mean1 = mean1, mean2 = mean2)###mean1, mean2
fit0 = lm(dat0$mean2 ~ dat0$mean1)
#plot(fit0)
#------------------------------------II. Outlier test on Avg.log2RPM--------------------------------------------
##"outlierTest" reports the Bonferroni p-values for testing each observation in turn to be a mean-shift outliner, 
#based Studentized residuals in linear (t-tests) and generalized linear models (normal tests).
tt0 = outlierTest(fit0, cutoff = Inf, n.max = nrow(dat0)) 
d1 = dat0[as.numeric(names(tt0[[1]])), ]
d2 = do.call(cbind.data.frame, tt0)
p.adj.FDR = p.adjust(d2$p, method = "BH")
m0 = cbind(d1, d2[,1:3], p.adj.FDR)
write.csv(m0, "Outlier_sgRNA_level_PIK3CA.csv")

#----------------------III. Optional: Scatter plot of Avg.log2RPM in tumors vs in cells (sgRNA level)-------------
ggplot(dat0, aes(x = mean1, y = mean2))+
  geom_point()+
  geom_smooth(method=lm,linetype = "dashed", color= "red")+ # Add linear regression line #  (by default includes 95% confidence region) 
  scale_size_manual(values = c(1,4)) +
  theme(panel.border = element_rect(colour = "black",size = 0.6,fill=NA), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.text.y  =element_text(size=26),
        axis.text.x =element_text(size=26,margin = margin(t = 1, r = 1, b = 1, l = 1)),
        axis.title=element_text(size=24),#face="bold")
        plot.title = element_text(face="bold",size=30,hjust=0.5)
        #strip.text.y  = element_text(size = 12, face="bold")
        ##strip.text is for facet
  )+
  labs(title = "",x = expression(paste(log[2]~RPM, " in cells")), y = expression(paste(log[2]~RPM, " in tumors")))+
  scale_x_continuous(limits = c(0,9), breaks = round(seq(0, 10, by = 1))) +
  scale_y_continuous(limits = c(-2,15),breaks = round(seq(-5, 15, by = 5),0))+
  geom_hline(yintercept=-4, linetype="dashed", color = "grey")
#ggsave("logRPMtumor.vs.cell_PIK3CA_20200908_Bootstrap_permutation.pdf", width = 7, height =8)

#---------------------------------IV. Examine how many dual-sgRNAs enriched--------------------------------------
rP = read.csv("Outlier_sgRNA_level_PIK3CA.csv", row.names = 1)
rAvg <- read.csv("Outlier_gene_level_PIK3CA.csv",row.names = 1)
sname = as.character(rownames(rAvg))
cfP = 0.1 ##cutoff for p value
rAvg$num.total.sg = 0
rAvg$num.sg.outlier = 0
rAvg$num.total.sg1 = 0
rAvg$num.sg1.outlier = 0
rAvg$num.total.sg2 = 0
rAvg$num.sg2.outlier = 0

for(i in 1:nrow(rAvg)){
  g1 <- paste0("^", unlist(strsplit(sname[i], split = "-")[1])[1])
  g2 <- paste0("-", unlist(strsplit(sname[i], split = "-")[1])[2])
  gg <- rownames(rP)[grepl(g1, rownames(rP)) & grepl(g2, rownames(rP))]
  rAvg$num.total.sg[i] <- length(gg)
  rAvg$num.sg.outlier[i] <- sum(as.numeric(rP[gg,"p.adj.FDR"]) < cfP)
  
  gg_outlier = gg[which(as.numeric(rP[gg,"p.adj.FDR"]) < cfP)]
  rAvg$num.total.sg1[i] = length(unique(unlist(lapply(gg, function(x) {strsplit(x, split = "-")[[1]][1]}))))
  rAvg$num.sg1.outlier[i]= length(unique(unlist(lapply(gg_outlier, function(x) {strsplit(x, split = "-")[[1]][1]}))))
  
  rAvg$num.total.sg2[i] = length(unique(unlist(lapply(gg, function(x) {strsplit(x, split = "-")[[1]][2]}))))
  rAvg$num.sg2.outlier[i]= length(unique(unlist(lapply(gg_outlier, function(x) {strsplit(x, split = "-")[[1]][2]}))))
  
}
rAvg$perP = rAvg$num.sg.outlier/rAvg$num.total.sg*100
write.csv(rAvg, "Number_outlier_sgRNAs_PIK3CA.csv")

#----------------------------------@@@@@@@@@@@@@@@@@@@@@@@@@@@@-------------------------------------------------
#============================Part 4: Calculate GI score based on LFC============================================
#----------------------------------@@@@@@@@@@@@@@@@@@@@@@@@@@@@-------------------------------------------------
#-----------------------------I. Determine GI score using Max, Min，Additive model------------------------------
ref <- read.csv("gene_comb_names_PIK3CA.csv",row.names = 1)
sname <- as.character(ref[,1])

read = read.csv('res_LFC_PIK3CA.csv', row.names = 1, header = T)
lfc1 = t(apply(read, 1, function(x){(x[3:5])/(x[1])}))
lfc2 = t(apply(read, 1, function(x){(x[3:5])/(x[2])}))
readM = cbind(lfc1, lfc2, read[,"LFC"]) # dim(readM): 1378, 19

SKO = sname[grepl('CTRL', sname)]
SKO = SKO[SKO != 'CTRL-CTRL']
DKO = sname[!grepl('CTRL', sname)]

readS = readM[SKO,]
readD = readM[DKO,]

GGIcompute = function(doubleK, singleM = readS) {
  sname = as.character(doubleK[1])
  lfc.gg = as.numeric(doubleK[length(doubleK)])
  g1 = strsplit(sname, split = '-')[[1]][1]
  g2 = strsplit(sname, split = '-')[[1]][2]
  lfc.g1 = as.numeric(singleM[grepl(g1,rownames(singleM)),ncol(singleM)][!is.na(singleM[grepl(g1,rownames(singleM)),ncol(singleM)])])
  lfc.g2 = as.numeric(singleM[grepl(g2,rownames(singleM)),ncol(singleM)][!is.na(singleM[grepl(g2,rownames(singleM)),ncol(singleM)])])
  ggi = differentModels(sname, lfc.gg, lfc.g1, lfc.g2)
  return (data.frame(ggi))
}

differentModels = function(str_name, num_d, num_s1, num_s2) {
  n1 = length(num_s1)
  n2 = length(num_s2)
  #totalName = NULL
  #totalLine = NULL
  if ((n1==0) | (n2==0)) {
    # NA for all single knockout
    naout = matrix(c(NA, NA, NA, NA, NA, NA, NA), nrow=1)
    colnames(naout) = c("LFC.DKO", "LFC.SKO1", "LFC.SKO2", "Add", "Max", "Min")
    return(naout)
  } else {
    #name_line = paste(str_name, names(num_s1), names(num_s2), sep="-") 
    add_line = num_d - (num_s1 + num_s2) 
    max_line = num_d - max(num_s1, num_s2)
    min_line = num_d - min(num_s1, num_s2)
    all_line = matrix(c(num_d, num_s1, num_s2, add_line, max_line, min_line), nrow = 1)
    #totalName = c(totalName, str_name)
    #totalLine = rbind(totalLine, all_line)
    rownames(all_line) = str_name
    colnames(all_line) = c("LFC.DKO", "LFC.SKO1", "LFC.SKO2", "Add", "Max", "Min")
  }
  return(all_line)
}
GGIvalues.0 = apply(cbind(rownames(readD), readD), 1, GGIcompute)
GGIvalues = do.call(rbind, GGIvalues.0)

#--------------------------II. Determine the statistical significance of GI scores------------------------------
##Construct the DKO and SKO matrices (ncol = 6 + 6 = 12, nrow = 1325 gene pairs)
LFC.max = function(doubleK, singleM = readS){
  sname = as.character(doubleK[1])
  g1 = strsplit(sname, split = '-')[[1]][1]
  g2 = strsplit(sname, split = '-')[[1]][2]
  lfc.g1 = as.numeric(singleM[grepl(g1,rownames(singleM)),7][!is.na(singleM[grepl(g1,rownames(singleM)),7])])
  lfc.g2 = as.numeric(singleM[grepl(g2,rownames(singleM)),7][!is.na(singleM[grepl(g2,rownames(singleM)),7])])
  gmax.sname = ifelse(lfc.g1 > lfc.g2, g1, g2)
  lfc.gg = as.numeric(doubleK[2:7])
  lfc.gmax = as.numeric(singleM[grepl(gmax.sname,rownames(singleM)),1:6])
  lfcM = matrix(c(lfc.gg, lfc.gmax), nrow = 1)
  colnames(lfcM) = c(paste0('DKO_',names(doubleK[2:7])), paste0('Max_SKO_',names(doubleK[2:7])))
  rownames(lfcM) = sname
  return(data.frame(lfcM))
}

LFC.Add = function(doubleK, singleM = readS){
  sname = as.character(doubleK[1])
  g1 = strsplit(sname, split = '-')[[1]][1]
  g2 = strsplit(sname, split = '-')[[1]][2]
  
  lfc.gg = as.numeric(doubleK[2:7])
  lfc.add = as.numeric(singleM[grepl(g1,rownames(singleM)),1:6]) + as.numeric(singleM[grepl(g2,rownames(singleM)),1:6])
  lfcM = matrix(c(lfc.gg, lfc.add), nrow = 1)
  colnames(lfcM) = c(paste0('DKO_',names(doubleK[2:7])), paste0('Add_SKO_',names(doubleK[2:7])))
  rownames(lfcM) = sname
  return(data.frame(lfcM))
}

lfcM.max.0 = apply(cbind(rownames(readD), readD), 1, LFC.max)
lfcM.max = do.call(rbind, lfcM.max.0)

lfcM.add.0 = apply(cbind(rownames(readD), readD), 1, LFC.Add)
lfcM.add = do.call(rbind, lfcM.add.0)

#Determine the significance of GIs by comparing 6 ratios from 3 tumor samples and 2 cell samples using Bootstrap_permutation methods
set.seed(42)
max.Boot = NULL
add.Boot = NULL
nBoot = 1000
for(i in 1:nBoot){
  idx = sample(1:ncol(lfcM.max), size = 3*ncol(lfcM.max), replace = T)
  idxD = idx[1:(length(idx)/2)]
  idxS = idx[(length(idx)/2+1) : length(idx)]
  temp1 = apply(lfcM.max, 1 , function(x) {log2(mean(x[idxD])/mean(x[idxS]))})
  max.Boot = cbind(max.Boot, temp1)
  
  temp2 = apply(lfcM.add, 1, function(x) {log2(mean(x[idxD])/mean(x[idxS]))})
  add.Boot = cbind(add.Boot, temp2)
}

#Two sided test
diff.max = apply(lfcM.max, 1, function(x){log2(mean(x[1:(ncol(lfcM.max)/2)])/mean(x[(ncol(lfcM.max)/2+1) : ncol(lfcM.max)]))})
diff.max.Boot = cbind(diff.max, max.Boot)
p.max = apply(diff.max.Boot, 1, function(x) {return (sum( abs(x[2:(nBoot+1)]) >= abs(x[1]))/nBoot) }   )
p.max.adj = p.adjust(p.max, method = 'BH')

#Two sided test
diff.add = apply(lfcM.add, 1, function(x){log2(mean(x[1:(ncol(lfcM.add)/2)])/mean(x[(ncol(lfcM.add)/2+1) : ncol(lfcM.add)]))})
diff.add.Boot = cbind(diff.add, add.Boot)
p.add = apply(diff.add.Boot, 1,  function(x){return (sum( abs(x[2:(nBoot+1)]) >= abs(x[1]))/nBoot)}) 
p.add.adj = p.adjust(p.add, method = 'BH')

zscore.Add = scale(GGIvalues$Add)
zscore.Max = scale(GGIvalues$Max)

GGI = cbind(GGIvalues, p.max, p.max.adj, zscore.Max, p.add, p.add.adj,zscore.Add)
write.csv(GGI, 'GI_score_PIK3CA.csv')

#----------------------------------@@@@@@@@@@@@@@@@@@@@@@@@@@@@-------------------------------------------------
#==========================Part 5: Combine outlier results with LFC and GI score=================================
#----------------------------------@@@@@@@@@@@@@@@@@@@@@@@@@@@@-------------------------------------------------
Outlier = read.csv("Number_outlier_sgRNAs_PIK3CA.csv", row.names = 1 )
dim(Outlier)
cN = paste0("OT_",colnames(Outlier))
colnames(Outlier) = cN
#--------------------------I. Combine LFC file and number of outliers file-------------------------------------
res = read.csv("res_LFC_PIK3CA.csv", row.names =1)
sum(!(rownames(res) %in% rownames(Outlier)))
gP1 = cbind.data.frame(res,Outlier[rownames(res),])
dim(gP1)
write.csv(gP1,"Res_LFC_Outlier_PIK3CA_1377_combs.csv" )
#--------------------------II. Combine GI file and number of outliers file-------------------------------------
GGi = read.csv("GI_score_PIK3CA.csv", row.names =1)
sum(!(rownames(GGi) %in% rownames(Outlier)))
dim(GGi)
gP2 = cbind.data.frame(gP1[rownames(GGi),], GGi)
dim(gP2)
write.csv(gP2,"Final_LFC_Outlier_GI_score_PIK3CA_1325_DKOs.csv")

#----------------------------------@@@@@@@@@@@@@@@@@@@@@@@@@@@@-------------------------------------------------
#=====================Part 6: Analysis or plotting after determination of LFC and GI score======================
#----------------------------------@@@@@@@@@@@@@@@@@@@@@@@@@@@@-------------------------------------------------

#----------------------------------------I. LFC_based_on_gene_rank----------------------------------------------
lfc_pik3 = as.data.frame(read.csv("Res_LFC_Outlier_PIK3CA_1377_combs.csv", header = T, row.names = 1))
dim(lfc_pik3)
lfc_pik3$sname = rownames(lfc_pik3)
lfc_pik3$hL = "not.hL"
lfc_pik3[which(lfc_pik3$LFC > 0 & lfc_pik3$p_adj < 0.05),"hL" ] = "hL"
lfc_pik3$Order[order(lfc_pik3$LFC, decreasing = T)] = 1:nrow(lfc_pik3)
data = lfc_pik3[,c("Order", "LFC", "p_adj","OT_perP", "hL")] 
sum(data$hL == "hL") #58

p <- ggplot(data, aes(x = Order, y = LFC))+
  geom_point(color= "grey60", size = 1.5)+
  geom_point(data = subset(data, hL == "hL"), color = "#FF3B30", size = 3, alpha = 1)+
  theme(#panel.border = element_rect(colour = "black",size = 0.6,fill=NA),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line = element_line(colour = "black"),
    legend.position = "none",
    #legend.title=element_text(size=18),
    legend.text=element_text(size=18),
    legend.title.align = 0.2,
    #legend.spacing.x = unit(2, 'cm'),
    legend.key.size = unit(1.5, "cm"),
    axis.text.y  =element_text(size=24,color = 'black'),
    axis.text.x =element_text(size=22, hjust = 0.8, margin = margin(c(5,0,0,0)),color = 'black'),
    axis.title=element_text(size=24),
    plot.title = element_text(size=22,hjust=0.5, vjust = -0.5),
    strip.text  = element_text(size = 16,color = 'black'))+
  ##strip.text is for facet)+
  labs(title = "PIK3CA",x = "Rank of gene pairs", y = expression(paste(log[2]~fold~change, " (LFC)")))+
  #scale_color_identity(name = "", labels = ""), paste0("")), guide = "legend")+
  scale_y_continuous(limits = c(-6,10),breaks = round(seq(-5, 10, by = 5),0)) +
  scale_x_continuous(limits = c(0,1500),breaks = round(seq(0, 1500, by = 500),0), labels = c("0", "0.5K", "1K", "1.5K"))+
  geom_hline(yintercept= 0, linetype="dashed", color = "grey60", size = 1)
p  
ggsave("LFC_based_on_gene_rank_PIK3CA.pdf", width = 5.2, height = 5.8)

#--------------------------- II. Scatter plot of Avg.log2RPM in tumors vs in cells--------------------------------------------
f1 = lfc_pik3[which(lfc_pik3$LFC > 0 & lfc_pik3$p_adj < 0.05 ), ]
hl.names = as.character(f1[order(f1[,"LFC"], decreasing = T), "sname"][1:5])
lfc_pik3$text.or.not = "not.texted"
lfc_pik3[hl.names, "text.or.not"] = "texted"
#Scatterplot of log2RPM in tumors versus log2RPM in cells
p <- ggplot(lfc_pik3, aes(x=OT_c.gg, y=OT_t.gg)) +
  geom_point( aes(size=hL, color=hL, alpha = text.or.not)) +    # Use hollow circles
  scale_alpha_manual(values= c(0.4,1),guide=F)+
  scale_size_manual(values = c(3,1.5)) +
  scale_color_manual(values= c("red","grey60"))+
  geom_smooth(method=lm,linetype = "dashed", color= "black")+ # Add linear regression line #  (by default includes 95% confidence region) 
  theme(panel.border = element_rect(colour = "black",size = 0.6,fill=NA), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.text.y  =element_text(size=24),
        axis.text.x =element_text(size=24),
        axis.title=element_text(size=24),#face="bold")
        plot.title = element_text(size=21,hjust=0.5)
        #strip.text.y  = element_text(size = 12, face="bold")
        ##strip.text is for facet
  )+
  labs(title = "PIK3CA",x = expression(paste(Avg.~log[2]~RPM, " in cells")), y = expression(paste(Avg.~log[2]~RPM, " in tumors")))+
  scale_x_continuous(limits = c(1,8),breaks = round(seq(1,8, by = 1),0)) +
  scale_y_continuous(limits = c(-1.5,5),breaks = round(seq(-1, 5, by = 1),0))+
  geom_hline(yintercept=1, linetype="dashed", color = "grey")
p
require("ggrepel")
set.seed(68)
p+ geom_text_repel(
  data = subset(lfc_pik3, text.or.not == "texted"),
  aes(label = sname),
  size = 6,
  color = "darkblue",
  box.padding = unit(0.3, "lines"),
  point.padding = unit(0.2, "lines"),
  #label.padding = unit(0.2, "lines"),
  #segment.color = 'grey50',
  nudge_y = 0.4,
  nudge_x = -0.15
)
ggsave("AvglogRPMtumor.vs.cell_PIK3CA.pdf", width = 5.2, height = 5.8)

#-------------------------------III. Top 10 tumor-promoting gene pairs----------------------------------------------
res = lfc_pik3[which(lfc_pik3$hL == "hL"),]
res = res[with(res, order(LFC, decreasing = T)),]
res0 = head(res, 10)
nM.SKO = unique(unlist(strsplit(rownames(res0), split = '-')))

Res = rbind.data.frame(res0, lfc_pik3[paste0(nM.SKO, '-CTRL'),])
Res = Res[order(Res$LFC, decreasing = T),]
dat = cbind.data.frame(sg = rownames(Res), Res)

dat$group = "DKO"
dat[grepl("CTRL", dat$sg),"group"] = "SKO"

dim(Res)
sgs = factor(dat$sg, levels = rev(dat$sg))
p <- ggplot(dat, aes(x = sgs, y=LFC,fill = group)) + #fill = group
  geom_bar(stat="identity",position=position_dodge(),width = 0.6,  color = "#007AFF")+ #fill = "#5AC8FA",
  geom_errorbar(aes(ymin=l, ymax=u), width=.2,position=position_dodge(0.9),color = "grey20")+
  scale_fill_manual(values = c("#5AC8FA","grey90"))+
  #facet_grid(~ tumor)+
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none',
        #legend.position = c(0.75,0.9),
        #legend.title = element_text(size = 20, face= "bold"),
        legend.text = element_text(size=20),
        axis.text.y  =element_text(size=16,color = 'black'),
        axis.text.x =element_text(size=20,color = 'black', margin = margin(5,0,0,0)),
        axis.title=element_text(size=22,margin = margin(10,0,0,0)),#face="bold")
        plot.title = element_text(size=23,hjust= 1, margin = margin(0,0,10,0))
        #strip.text  = element_text(size = 16, face="bold")
        ##strip.text is for facet
  )+
  scale_y_continuous(limits = c(-5,10), breaks = seq(-4,10,2))+
  #labs(title="",x ="", y = "Percentage(%) of \n log2 RPM > 0 ")+
  labs(title="Top 10 tumor-promoting DKOs",x ="", y = "LFC (PIK3CA)")+
  coord_flip()+
  geom_hline(yintercept=0, color = "grey5")

p
ggsave('Top10_tumor_promoting_DKOs_PIK3CA.pdf', width = 5.6, height = 6.2)

#-------------------------------IV. Correlation of LFC and outlier_residual----------------------------------------------
library("ggpubr")
ggscatter(lfc_pik3, x = "LFC" , y = "OT_rstudent", size = 3,alpha = 0.5,
          add = "reg.line", conf.int = FALSE, 
          add.params = list(color = "#007AFF", fill = "lightgray"), # Customize reg. line
          cor.coef = TRUE,rug = TRUE,
          cor.coeff.args = list(method = "pearson", label.x = -3.5,label.y = 13, label.sep = "\n", size = 7))+
  theme(#panel.border = element_blank(), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    #axis.ticks.y = element_line(colour = "black"),
    legend.position = 'none',
    #legend.position = c(0.75,0.9),
    #legend.title = element_text(size = 20, face= "bold"),
    legend.text = element_text(size=20),
    axis.text.y  =element_text(size=25,color = 'black'),
    axis.text.x =element_text(size=25,color = 'black', margin = margin(5,0,0,0)),
    axis.title=element_text(size=24,margin = margin(10,0,0,0)),#face="bold")
    plot.title = element_text(size=23,hjust= 1, margin = margin(0,0,10,0))
    #strip.text  = element_text(size = 16, face="bold")
    ##strip.text is for facet
  )+
  scale_x_continuous(limits = c(-4,8), breaks = seq(-5,7.5,2.5))+
  scale_y_continuous(limits = c(-2,15), breaks = seq(-2,16,2))+
  labs(title="",x ="LFC (PIK3CA)", y = "Studentized residual")

ggsave("Correlation_outlier_analysis_LFC_PIK3CA.pdf", width = 6, height =5)

#-------------------------------V. Tumor-promoting GGIs will be selected and plotted----------------------------------------------
ggi = read.csv("Final_LFC_Outlier_GI_score_PIK3CA_1325_DKOs.csv", row.names = 1, header = T)
ggi = as.data.frame(ggi)
ggi$sname = rownames(ggi)
ggi$oncogenic = "Not_oncogenic"

##How many oncogenic gene pairs
ggi[which(ggi$LFC > 0 & ggi$p_adj < 0.05),"oncogenic"] = "Yes_oncogenic"   ##55
nrow(ggi[which(ggi$oncogenic == "Yes_oncogenic"),]) #33

##How many tumor promoting positive interactions using "Add model"
ggi.add = ggi[which((ggi$oncogenic == "Yes_oncogenic") & (ggi$Add > 0) & (ggi$p.add.adj < 0.05) ),]
nrow(ggi.add) #27

##How many tumor promoting positive interactions using "Max model"
ggi$hL = 0
ggi[which((ggi$LFC > 1.7 & ggi$p_adj < 0.05) & (ggi$Max > 1.5) & (ggi$p.max.adj < 0.05)),"hL"] = 1 ##15
nrow(ggi[which(ggi$hL == 1),])
d0 = ggi[which(ggi$hL == 1),]
write.csv(d0, "20_Tumor_promoting_GIs_PIK3CA.csv")

p = ggplot(ggi, aes(x = Max, y = Min))+
  geom_point(aes(size=factor(hL), color=factor(oncogenic), alpha = factor(hL))) +    # Use hollow circles
  scale_alpha_manual(values= c(0.2,1), guide=F)+
  scale_size_manual(values = c(2,4)) +
  scale_color_manual(values= c("grey","darkorchid1"))+
  theme(panel.border = element_rect(colour = "black",size = 0.6,fill=NA), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.text.y  =element_text(size=24, color = "black"),
        axis.text.x =element_text(size=24, color = "black"),
        axis.title=element_text(size=24, color = "black"),#face="bold")
        plot.title = element_text(size=24,hjust=0.5, color = "black")
        #strip.text.y  = element_text(size = 12, face="bold")
        ##strip.text is for facet
  )+
  labs(title="PIK3CA",x =expression(paste("Diff: ",LFC[(DKO)]," - ", LFC[(SKO1)])), y = expression(paste("Diff: ",LFC[(DKO)]," - ", LFC[(SKO2)])))+
  scale_x_continuous(limits =c(-6,10.5), breaks = seq(-5,10,2.5)) +
  scale_y_continuous(limits= c(-4,10),breaks = seq(-2.5,10,2.5))+
  geom_hline(yintercept= 0, linetype="dashed", color = "grey")+
  geom_vline(xintercept= 0, linetype="dashed", color = "grey")
p
require("ggrepel")
set.seed(42)
data = subset(ggi, hL == 1)
d0 = data[with(data, order(Max, decreasing = T)),]
p+ geom_text_repel(
  data = head(d0,6), 
  aes(label = sname),
  size = 6,
  color = "darkblue",
  box.padding = unit(0.6, "lines"),
  point.padding = unit(0.28, "lines"),
  nudge_y = 0.4,
  nudge_x = 2.2
) 
ggsave("Max_Min_GGi_PIK3CA_tumors.pdf", width=6.2, height =6.8)

#--------------------------VI.Tumor-promoting GIs are selected for cytoscape (GI network map)---------------------------------------------
library(igraph)
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
read0 = read.csv("20_Tumor_promoting_GIs_PIK3CA.csv", row.names = 1, header = T)
dim(read0)
gene1 <- sapply(rownames(read0), function(x){unlist(strsplit(x,"-")[[1]][1])} )
gene2 <- sapply(rownames(read0), function(x){unlist(strsplit(x,"-")[[1]][2])} )

PPI_net =  cbind.data.frame(gene1, gene2, ggi = read0[,"Max"])
write.csv(PPI_net, "GI.score_PIK3CA_cytoscape.csv", quote = FALSE)

net = graph_from_data_frame(d = PPI_net, directed = F)
deg = igraph::degree(net)
deg_order = deg[order(deg, decreasing = F)]
V(net)$size =deg*3# Node size stands for the degree of protein
V(net)$color = apple_colors[6]
#V(net)$label = NA
E(net)$width <- E(net)$ggi*2
E(net)$color[E(net)$ggi !=0 ] = apple_colors[8]
#E(net)$color[E(net)$Significant == "Not Sig" ] = apple_colors[10]
l=layout_in_circle(net)
#l=layout_as_star(net)
#l = layout_nicely(net)
#l = layout_on_sphere
#l = layout_with_dh(net)
#l = norm_coords
net_clean <- delete.edges(net, which(E(net)$PTEN == 0))
net_clean <- delete.edges(net, which(E(net)$PTEN== 0))
plot(net_clean, layout = l,  vertex.frame.color=NA, vertex.shape="circle",vertex.label.color = apple_colors[11],
     vertex.label.cex = 0.8, margin = c(0,0,0,0),main = "PIK3CA", vertex.label.family= "Arial")

#pdf("22_oncogenic_GIs_PIK3CA_20200908.pdf", height =6, width = 6)
#plot(net_clean, layout = l,  vertex.frame.color=NA, vertex.shape="circle",vertex.label.color = apple_colors[11],
      #vertex.label.cex = 0.78, margin = c(0,0,0,0),main = "PIK3CA", vertex.label.family= "Helvetica")
#dev.off()

#Centrality score of nodes
g = net
cent_scores <- data.frame(
  degree = degree(g),
  betweenness = round(betweenness(g),4),
  closeness = round(closeness(g),4),
  eigenvector = round(eigen_centrality(g)$vector,4),
  subgraph = round(subgraph_centrality(g),4))
write.csv(cent_scores, "Centrality_scores_PIK3CA.csv")

#-------------------------------@@@@@@@@@@@@@@@@@@@@@@@@@@@@----------------------------------------------
#=====================================Additional analysis==============================================
#-------------------------------@@@@@@@@@@@@@@@@@@@@@@@@@@@@----------------------------------------------
#-------------------VII. Percentage of dual-sgRNAs enriched of tumor-promoting gene pairs-----------------------
library(dplyr)
library(scales)
sum = 11+8+1   #20
dat = data.frame(group = c("g1", "g2", "g3"),
                 num = c(11,8,1),
                 #num = c(7,13,17,8,1),
                 per = round(c(11,8,1)/20*100, 1)
                 #lab.ypos = cumsum(dat$per) - 0.8*dat$per
)

#lab.ypos = cumsum(dat$per) - 0.5*dat$per
#colfunc <- colorRampPalette(c("#EFEFF4","#007AFF"))
#colfunc(40)
#plot(rep(1,40),col=colfunc(40),pch=19,cex=3)
#mycols <- c("magenta", "blue","red","yellow", "green")
#mycols <- colfunc(40)
ggplot(dat, aes(x = 2, y = per, fill = group))+
  geom_bar(width = 1,stat = "identity", color = "white")+
  coord_polar("y", start = 0)+
  geom_text(aes(label = num), color = "black", size = 8, position = position_stack(vjust = 0.5))+
  xlim(1.5, 2.5)+
  theme_void()+
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        #axis.line = element_line(colour = "black"),
        legend.position = c(1.2,0.6),
        legend.title = element_text(size = 20, face= "bold"),
        legend.text = element_text(size=18 ),
        #axis.text.y  =element_text(size=24,face= "bold"),
        #axis.text.x =element_text(size=24, face = "bold"),#, margin = margin(t =20, r = 10, b = 10, l = 10)),
        #axis.title=element_text(size=30),#face="bold")
        plot.title = element_text(size=22,hjust=0.5)
        #strip.text.y  = element_text(size = 26, face="bold")
        ##strip.text is for facet
  )+
  scale_fill_brewer(palette =  "Pastel2",
                    name="",
                    #breaks=c("Cells", "Tumors"),
                    labels=c("20-40%", "40-60%", "80-100%"))+
  labs(title="PIK3CA",x ="", y = "")
ggsave("Percentage_of_dual-sgRNAs_enriched_in_20_oncogenic_GIs_PIK3CA.pdf", width=10, height = 5)

#-------------------------VIII. Deseq & tumor sample distance for heatmap--------------------------------------
library(DESeq2)
library(gplots)
library(pheatmap)
counts = read.csv("cell_tumor_counts_PIK3CA.csv", header=T, row.names=1)
condition <- c(rep("Cell",2),rep("Tumor",(ncol(counts)-2)))
condition
coldata <- data.frame(condition,row.names = colnames(counts))
dds<- DESeqDataSetFromMatrix(countData = counts,
                             colData = coldata,
                             design = ~ condition)
keep <- rowSums(counts(dds[,1:2])) >= 1
dds <- dds[keep,]
dds$condition <- factor(dds$condition, levels = c("Cell","Tumor"))
#dds <- DESeq(dds)
rld <- rlog( dds )
sampleDists <- dist( t( assay(rld) ) )
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )

pdf("Heatmap_sample_distances_PIK3CA.pdf",width = 8, height=8)
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(128)
heatmap.2(sampleDistMatrix, trace="none", col=colours, margins=c(10,10),cexRow=3.4,cexCol=3.4)
#pheatmap(sampleDistMatrix, color = colours, border_color = "NA", fontsize = 12)
dev.off()

#--------------------IX. Dotplot showing DKO dual-sgRNAs dominate the cell population in tumors-----------------------------
library(tidyr)
raw_cts = read.csv("cell_tumor_counts_PIK3CA.csv", row.names = 1) ###without two columns of mean_log2RPM
keep <- rowSums(raw_cts[,1:2]) > 0
cts = raw_cts[keep,]
RPM = as.matrix(apply(cts[,1:ncol(cts)], 2, function(x) {x/sum(x)*1000000})) 
dim(RPM)
write.csv(RPM, "RPM_PIK3CA.csv")
rpmM = as.data.frame(apply(RPM[,1:ncol(RPM)], 2, function(x) {log2(x+1)})) ###0.01 or 0.001
apply(rpmM, 2, function(x) sum(x>10))

messy = cbind.data.frame(dual.sgRNAs = rownames(rpmM), rpmM)
messy$sg_type = "NTC"

for(row in 1:nrow(messy)){
  splitName =  strsplit(as.character(messy[row,"dual.sgRNAs"]), split = "-")[[1]]
  if(sum(grepl("CTRL", splitName)) ==0){
    messy[row,"sg_type"]= "DKO"
  } else if(sum(grepl("CTRL", splitName)) == 1){
    messy[row,"sg_type"] = "SKO"
  }
}
#which(messy1$sg_type == "NTC")
#messy1[35472,]
apply(messy[,2:12], 2, function(x) {sum(x > 10)})
dim(messy) ##32883     7

df = gather(messy, samples, log2RPM, -dual.sgRNAs, -sg_type)
dim(df) #164415  4  (32883 * 5 samples)
sg_type1 = factor(df$sg_type, levels = c("DKO", "SKO", "NTC"))
sample = factor(df$samples, levels = c("Cell.1", "Cell.2", "F995.1", "F995.2","F996.1"))
ggplot(df, aes(x = sample, y = log2RPM, color = sg_type1))+
  geom_jitter(aes(color = sg_type1, alpha = sg_type1), position=position_jitter(0.25), size = 1)+
  scale_color_manual(values = c("orchid1", "blue", "green"))+
  scale_alpha_manual(values = c(0.3, 1,0.6))+
  theme(#panel.border = element_blank(), 
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_blank(),
    #panel.grid= element_line(),
    axis.line = element_line(colour = "black"),
    #legend.position = c(0.11,0.87),
    legend.title = element_text(size = 20, face= "bold"),
    legend.text = element_text(size=18),
    axis.text.y  =element_text(size=20, color = "black"),
    axis.text.x =element_text(size=20,color = "black",angle = 90),
    axis.title=element_text(size=22,color = "black"),#face="bold")
    plot.title = element_text(size=32,hjust=0.5)
    #strip.text.y  = element_text(size = 26, face="bold")
    ##strip.text is for facet
  )+
  scale_y_continuous(limits = c(-2,20), breaks = seq(0,20,5) )+
  labs(title="",x ="", y = expression(log[2]~RPM), color = "") +
  geom_hline(yintercept = c(0,10), linetype = "dashed", color = "darkslategray")+
  guides(alpha = FALSE) ##not show legend for alpha
#scale_fil8
ggsave("SKO_DKO_NTC_individual_tumors_PIK3CA.pdf", width = 4.8, height = 5.8)

