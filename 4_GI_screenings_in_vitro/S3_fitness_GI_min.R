#Xiaoyu Zhao
RNGkind(sample.kind = "Rounding")
# RNGkind is to get the results from set.seed() to match for different R versions
library(readxl)
library(gplots)
library(ggplot2)
library(ggpubr)
library(reshape)
library(RColorBrewer)
library(scales)  #For trans_breaks
setwd("~/Downloads/data/")
ref <- data.frame(read_excel("TableS1_Annotatation_of_asssyed_TSGs_or_TSG_combinations.xlsx", 5))
gene52 = as.character(unique(res[,1]))

#----------------------------------------------@@@@@@@@@@@@@@@@@@@---------------------------------------------------------
#==============================================PART 1: Fitness measurement=================================================
#----------------------------------------------@@@@@@@@@@@@@@@@@@@---------------------------------------------------------

#-----------------------------------------I. Measure fitness = lm(log2RPM ~ time)------------------------------------------
#Input: sheet3 in TableS5 ---> min.1; sheet4 in TableS5 --->min.2
cond = "min.2"
read = data.frame(read_excel("TableS5_Raw_and_processed_data_combinatorial_screening_in_vitro.xlsx",4))
rownames(read) = read[,"gene_names"]
dat00 = as.matrix(read[,3:ncol(read)])
# 1. filter some bad sgRNAs
filtered_list_others = c("ARID5B_2", "ARID5B_3", "PTEN_1", "PTEN_2", "PTEN_3")
filter_flag_1 = sapply(filtered_list_others, function(x) {grepl(x, rownames(dat00))})
filter_flag_2 = apply(filter_flag_1, 1, sum)
filter_flag = ifelse(filter_flag_2==0, TRUE, FALSE)
# 2. filter rows with low counts
dat = dat00[filter_flag, ]
t1 <- dat[, 1]>9
t2 <-apply(dat,1, function(x){ sum(x > 0) >=2})     
dat0 = dat[t1&t2,]
# 3. calculate RPM and log2RPM
rpm = apply(dat0+1, 2, function(x) {x/sum(x)*1000000})
logrpm = apply(dat0+1, 2, function(x) {log2(x/sum(x)*1000000)})
# 4. fitness = lm(log2RPM ~ time)
x = c(0, 3, 6, 9, 12, 15)
fitgrowth_lm = function(y) {
  outlm = coef(lm(y ~ x))
}
fitlm = t(apply(logrpm, 1, fitgrowth_lm))
#Save results
write.csv(cbind(rpm, fitlm), file=paste("minimal/fitness/",cond, ".RPM_fitlm.csv", sep=""))

##----------------------------------------I. lineage trajectories --------------------------------------------------
#------------------------------------------II.Plot lineage trajectories--------------------------------------------------
##function tranforming the matrix into a dataframe for plotting
traj_Transformation <- function(x, sname_pos, D0_pos, fitlm_pos){
  time_points = rep(c(0,3,6,9,12,15), nrow(x))
  traj = c(t(x[, D0_pos : (D0_pos+5)]))
  fitlm = rep(0, length(time_points))
  sname = rep("0", length(time_points))
  
  for(i in 1:nrow(x)){
    sname[(i*6-5): (i*6)] = as.character(rep(x[i,sname_pos], 6))
    fitlm[(i*6-5): (i*6)] = rep(x[i,fitlm_pos], 6)
  }
  trajM = data.frame(sname, fitness = fitlm, time_points, traj)
  return (trajM)
}
##function for plotting linear trajectories
traj_Plot = function(trajM, output){
  ggplot(data = trajM, aes(x = time_points, y= traj, group = sname, color = fitness)) +
    ggtitle("Minimal")+
    geom_line(alpha = 1, size = 0.3) +
    scale_color_gradientn(colors = c("#5AC8FA", "grey80", "magenta"), limits = c(-0.5, 0.5)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.text.x = element_text(size = 18, color = "black"), 
          axis.text.y = element_text(size = 18, color = "black"))+
    theme(text = element_text(size=20))+
    scale_y_continuous(name = "RPM",
                       limits = c(5e-2, 1e+4),
                       trans = log10_trans(),
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x))
    )+
    annotation_logticks(sides = "l") +  ## sides = "l" , annotation_logticks only shown on y axis
    scale_x_continuous(name = "Days",breaks = seq(0,15, by =3),labels = seq(0,15, by= 3))
  
  ggsave(output,  width = 6, height = 4)     
}

set.seed(250)
cond = "min.1"
fitlm = read.csv(file=paste("minimal/fitness/",cond, ".RPM_fitlm.csv", sep=""), header =T)
range(fitlm[,9])
dim(fitlm)
#Remove some noise, and show representative lineages
rd1 = fitlm[which(abs(fitlm[,9]) > 0.02),]
dim(rd1)
rd2 = rd1[which((rd1[,9]>0 & rd1[,3] > rd1[,2] & rd1[,4] > rd1[,3] & rd1[,5] > rd1[,4] & rd1[,6] > rd1[,5]) | (rd1[,9]<0 & rd1[,3] < rd1[,2] & rd1[,4] < rd1[,3] & rd1[,5] < rd1[,4] )) ,]
dim(rd2)
ns<-sample(nrow(rd2),5000)
read0<-rd2[ns,]
trajM =  traj_Transformation(read0, sname_pos = 1, D0_pos = 2, fitlm_pos = 9)
traj_Plot(trajM, paste("minimal/fitness/",cond, ".lineage_trajectories.pdf", sep=""))

#------------------------------III.Calculate percentage of RPM of DKO, SKO and CTRL at each timepoint--------------------------------------
#Function "perM_cal" is to calculate the relative percentage of RPM of DKOs, SKOs and CTRL in each replicate
perM_cal = function(cond){
  rpm = read.csv(file=paste("minimal/fitness/",cond, ".RPM_fitlm.csv", sep=""), header =T, row.names = 1)
  double.CTRL = grepl("^CTRL", row.names(rpm)) & grepl("-CTRL", row.names(rpm))
  sum(double.CTRL)
  rpm_KO = rpm[!double.CTRL,]
  
  single <- apply(rpm_KO[grepl("CTRL",rownames(rpm_KO)),1:6], 2, function(x) {sum(x)/sum(grepl("CTRL",row.names(rpm_KO)))})
  double <- apply(rpm_KO[!grepl("CTRL",rownames(rpm_KO)),1:6], 2, function(x) {sum(x)/sum(!grepl("CTRL",row.names(rpm_KO)))})
  ctrl <- apply(rpm[double.CTRL,1:6], 2, function(x) {sum(x)/sum(double.CTRL)})
  
  totalM = rbind(single, double, ctrl, single+double+ctrl)
  SKO = apply(totalM, 2, function(x){x[1]/x[4]*100})
  DKO = apply(totalM, 2, function(x){x[2]/x[4]*100})
  CTRL = apply(totalM, 2, function(x){x[3]/x[4]*100})
  
  perM = rbind(SKO, DKO, CTRL)
  return(perM)
}
perM.1 = perM_cal(cond = "min.1")
perM.2 = perM_cal(cond = "min.2")

perM.mean  = matrix(mapply(function(x,y) mean(c(x,y)), perM.1, perM.2), ncol=ncol(perM.1))
colnames(perM.mean) = c("t1", "t2" ,"t3", "t4", "t5", "t6")
rownames(perM.mean) = c("SKO.mean", "DKO.mean", "CTRL.mean")

perM.se  = matrix(mapply(function(x,y) sd(c(x,y))/sqrt(2), perM.1, perM.2), ncol=ncol(perM.1))
colnames(perM.se) = c("t1", "t2" ,"t3", "t4", "t5", "t6")
rownames(perM.se) = c("SKO.se", "DKO.se", "CTRL.se")

perM = rbind.data.frame(perM.mean, perM.se)
write.csv(perM, "minimal/fitness/min_relative_RPM_percentage_SKO_DKO_CTRL.csv")

#------------------------------V.Measure mean fitness of 1379 DKO and SKO, CTRL gene pairs--------------------------------------
#Organize fitness of gene-pair level at each row
set.seed(42)
ref <- data.frame(read_excel("TableS1_Annotatation_of_asssyed_TSGs_or_TSG_combinations.xlsx", 1))
sname = as.character(ref[,1])

totalcond = c("min.1", "min.2")
for (cond in totalcond) {
  read0 = read.csv(paste("minimal/fitness/",cond, ".RPM_fitlm.csv", sep=""), header=T, row.names=1)
  
  fitM= matrix(NA, nrow = length(sname),  ncol = 50 )
  rownames(fitM) = sname
  colnames(fitM) = rep("fitlm", 50)
  
  for(i in 1:length(sname)){
    if(sname[i] != "CTRL-CTRL"){
      g1 = strsplit(sname[i], split = '-')[[1]][1]
      g2 = strsplit(sname[i], split = '-')[[1]][2]
      
      target = read0[grepl(g1, rownames(read0)) & grepl(g2, rownames(read0)),8]
      fitM[i,1:length(target)] = target
    }
    else {
      target = read0[grepl("^CTRL", rownames(read0)) & grepl("-CTRL", rownames(read0)),8]
      fitM[i,1:length(target)] = target
    }
    
  }
  write.csv(fitM, file=paste("minimal/fitness/", cond, ".fitlm.Matrix.1379_gene_pairs.csv", sep=""))
}

#Construct fitness matrix for bootstrapping or permutation
#---------------------------------------------------------------------------------------------------------------------------------
cond = "min.1"
fitM0.1 = read.csv(file=paste("minimal/fitness/", cond, ".fitlm.Matrix.1379_gene_pairs.csv", sep=""), row.names = 1, header = T)
colSums(!is.na(fitM0.1))
fitM1.1 = as.matrix(fitM0.1[,1:49])
dim(fitM1.1)
ctrl.1 = t(fitM0.1[which(rownames(fitM0.1) == "CTRL-CTRL"), 1:15]) #Discard NAs
CTRL.1 = apply(ctrl.1, 1, function(x){return (rep(x, nrow(fitM0.1)))})
dim(CTRL.1)

cond = "min.2"
fitM0.2 = read.csv(file=paste("minimal/fitness/", cond, ".fitlm.Matrix.1379_gene_pairs.csv", sep=""), row.names = 1, header = T)
colSums(!is.na(fitM0.2))
fitM1.2 = as.matrix(fitM0.2[,1:49])
ctrl.2 = t(fitM0.2[which(rownames(fitM0.2) == "CTRL-CTRL"),1:17]) #Discard NAs
CTRL.2 = apply(ctrl.2, 1, function(x){return (rep(x, nrow(fitM0.2)))})

fitM = cbind(fitM1.1, fitM1.2[rownames(fitM1.1),], CTRL.1, CTRL.2)
dim(fitM)
write.csv(fitM, file="minimal/fitness/min.combined.fitMatrix_1379_gene_pairs.csv")

#---------------------------------------------------------------------------------------------------------------------------------
nBoot = 1000
cond = "min.1"
#Bootstrapping to calculate the 95% CI, and significance of mean fitness (compared with 0) for min.1
Boot1 = NULL
for(i in 1:nBoot){
  idx = sample(1:49, size = 49, replace = T)
  temp = apply(fitM, 1, function(x) mean(x[idx][!is.na(x[idx])]))
  Boot1 = cbind(Boot1, temp)
}

dim(Boot1)
fitlm.1 = apply(fitM, 1, function(x) mean(x[1:49][!is.na(x[1:49])]))
se.1 = apply(fitM, 1, function(x) sd(x[1:49][!is.na(x[1:49])])/(sqrt(length(x[1:49][!is.na(x[1:49])]))))

l1 <- apply(Boot1, 1, function(x){quantile(x[!is.na(x)], 0.025)})
u1<- apply(Boot1, 1, function(x){quantile(x[!is.na(x)], 0.975)})

f.Boot1 = cbind(fitlm.1, Boot1)
p_f.avg_vs_0_1 = apply(f.Boot1, 1, function(x) {if(x[1] > 0) {sum(x[2:(nBoot+1)] <= 0)/nBoot} else {sum(x[2:(nBoot+1)] >= 0)/nBoot}})
p.adj_1 = p.adjust(p_f.avg_vs_0_1, method = 'BH')
res1 = cbind.data.frame(f.avg = fitlm.1, f.se = se.1, l= l1, u =u1, p_f.avg_vs_0 = p_f.avg_vs_0_1, p.adj_I = p.adj_1)

#Bootstrapping to calculate the 95% CI, and significance of mean fitness (compared with fitness of CTRL) for min.1
idxF1 = 1:49
idxR1 = 99:(98+nrow(ctrl.1))
diff1 = apply(fitM, 1, function(x) mean(x[idxF1][!is.na(x[idxF1])]) - mean(x[idxR1][!is.na(x[idxR1])]) )
ctrl1.mean = apply(fitM, 1, function(x)  mean(x[idxR1][!is.na(x[idxR1])]) )

diff.Boot1 = NULL
nBoot = 1000
for(i in 1:nBoot){
  pool = c(1:49,99:(98+nrow(ctrl.1)))
  idx = sample(pool, size = 3*length(pool), replace = T)
  idxF = idx[1:(49*3)]
  idxR = idx[(49*3+1) : length(idx)]
  temp = apply(fitM, 1 , function(x) {mean(x[idxF][!is.na(x[idxF])]) - mean(x[idxR][!is.na(x[idxR])])})
  diff.Boot1 = cbind(diff.Boot1, temp)
}

diff.Matrix1 = cbind(diff1, diff.Boot1)
p.Boot1 = apply(diff.Matrix1, 1, function(x) {return (sum( abs(x[2:(nBoot+1)]) >= abs(x[1]))/nBoot) }   )
p.Boot1.adj = p.adjust(p.Boot1, method = 'BH')

RES1 = cbind.data.frame(res1, f.ctrl = ctrl1.mean, diff = diff1, p_f.avg_vs_CTRL = p.Boot1, p.adj_II = p.Boot1.adj)
write.csv(RES1, file=paste("minimal/fitness/", cond, ".fitlm_se_CI_diff_p_1379_gene_pairs.csv", sep=""))

#---------------------------------------------------------------------------------------------------------------------------------
cond = "min.2"
#Bootstrapping to calculate the 95% CI, and significance of mean fitness (compared with 0) for min.2
Boot2 = NULL
for(i in 1:nBoot){
  idx = sample(50:98, size = 49, replace = T)
  temp = apply(fitM, 1, function(x) mean(x[idx][!is.na(x[idx])]))
  Boot2 = cbind(Boot2, temp)
}

dim(Boot2)
fitlm.2 = apply(fitM, 1, function(x) mean(x[50:98][!is.na(x[50:98])]))
se.2 = apply(fitM, 1, function(x) sd(x[50:98][!is.na(x[50:98])])/(sqrt(length(x[50:98][!is.na(x[50:98])]))))

l2 <- apply(Boot2, 1, function(x){quantile(x[!is.na(x)], 0.025)})
u2<- apply(Boot2, 1, function(x){quantile(x[!is.na(x)], 0.975)})

f.Boot2 = cbind(fitlm.2, Boot2)
p_f.avg_vs_0_2 = apply(f.Boot2, 1, function(x) {if(x[1] > 0 ) {sum(x[2:(nBoot+1)] <= 0)/nBoot} else {sum(x[2:(nBoot+1)] >= 0)/nBoot}})
p.adj_2 = p.adjust(p_f.avg_vs_0_2, method = 'BH')
res2 = cbind.data.frame(f.avg = fitlm.2, f.se = se.2, l= l2, u =u2, p_f.avg_vs_0 = p_f.avg_vs_0_2, p.adj_I = p.adj_2)

#Bootstrapping to calculate the 95% CI, and significance of mean fitness (compared with fitness of CTRL) for min.2
idxF2 = 50:98
idxR2 = (99+nrow(ctrl.1)):ncol(fitM)
diff2 = apply(fitM, 1, function(x) mean(x[idxF2][!is.na(x[idxF2])]) - mean(x[idxR2][!is.na(x[idxR2])]) )
ctrl2.mean = apply(fitM, 1, function(x)  mean(x[idxR2][!is.na(x[idxR2])]) )

diff.Boot2 = NULL
nBoot = 1000
for(i in 1:nBoot){
  #pool = c(50:98,114:(113+nrow(ctrl.2)))
  pool = c(50:98,(99+nrow(ctrl.1)):ncol(fitM))
  idx = sample(pool, size = 3*length(pool), replace = T)
  idxF = idx[1:(49*3)]
  idxR = idx[(49*3+1) : length(idx)]
  temp = apply(fitM, 1 , function(x) {mean(x[idxF][!is.na(x[idxF])])- mean(x[idxR][!is.na(x[idxR])])})
  diff.Boot2 = cbind(diff.Boot2, temp)
}

diff.Matrix2 = cbind(diff2, diff.Boot2)
p.Boot2 = apply(diff.Matrix2, 1, function(x) {return (sum( abs(x[2:(nBoot+1)]) >= abs(x[1]))/nBoot) }   )

p.Boot2.adj = p.adjust(p.Boot2, method = 'BH')

RES2 = cbind.data.frame(res2, f.ctrl = ctrl2.mean, diff = diff2, p_f.avg_vs_CTRL = p.Boot2, p.adj_II = p.Boot2.adj)
write.csv(RES2, file=paste("minimal/fitness/", cond, ".fitlm_se_CI_diff_p_1379_gene_pairs.csv", sep=""))

#---------------------------------------------------------------------------------------------------------------------------------
cond = "min.combined"

#Bootstrapping to calculate the 95% CI, and significance of mean fitness (compared with 0) for min.1 and min.2
Boot3 = NULL
for(i in 1:nBoot){
  idx = sample(1:98, size = 98, replace = T)
  temp = apply(fitM, 1, function(x) mean(x[idx][!is.na(x[idx])]))
  Boot3 = cbind(Boot3, temp)
}

dim(Boot3)
fitlm.3 = apply(fitM, 1, function(x) mean(x[1:98][!is.na(x[1:98])]))
se.3 = apply(fitM, 1, function(x) sd(x[1:98][!is.na(x[1:98])])/(sqrt(length(x[1:98][!is.na(x[1:98])]))))

l3 <- apply(Boot3, 1, function(x){quantile(x[!is.na(x)], 0.025)})
u3<- apply(Boot3, 1, function(x){quantile(x[!is.na(x)], 0.975)})

f.Boot3 = cbind(fitlm.3, Boot3)
p_f.avg_vs_0_3 = apply(f.Boot3, 1, function(x) {if(x[1] > 0 ) {sum(x[2:(nBoot+1)] <= 0)/nBoot} else {sum(x[2:(nBoot+1)] >= 0)/nBoot}})
p.adj_3 = p.adjust(p_f.avg_vs_0_3, method = 'BH')
res3 = cbind.data.frame(f.avg = fitlm.3, f.se = se.3, l= l3, u =u3, p_f.avg_vs_0 = p_f.avg_vs_0_3, p.adj_I = p.adj_3)

#Bootstrapping to calculate the 95% CI, and significance of mean fitness (compared with fitness of CTRL) for min.1 and min.2
idxF3 = 1:98
idxR3 = 99:ncol(fitM)
diff3 = apply(fitM, 1, function(x) mean(x[idxF3][!is.na(x[idxF3])]) - mean(x[idxR3][!is.na(x[idxR3])]) )
ctrl3.mean = apply(fitM, 1, function(x)  mean(x[idxR3][!is.na(x[idxR3])]) )

diff.Boot3 = NULL
nBoot = 1000
for(i in 1:nBoot){
  idx = sample(1:ncol(fitM), size = 3*ncol(fitM), replace = T)
  idxF = idx[1:(98*3)]
  idxR = idx[(98*3+1) : length(idx)]
  temp = apply(fitM, 1 , function(x) {mean(x[idxF][!is.na(x[idxF])])- mean(x[idxR][!is.na(x[idxR])])})
  diff.Boot3 = cbind(diff.Boot3, temp)
}

diff.Matrix3 = cbind(diff3, diff.Boot3)
p.Boot3 = apply(diff.Matrix3, 1, function(x) {return (sum( abs(x[2:(nBoot+1)]) >= abs(x[1]))/nBoot) }   )
p.Boot3.adj = p.adjust(p.Boot3, method = 'BH')

RES3 = cbind.data.frame(res3, f.ctrl = ctrl3.mean, diff = diff3, p_f.avg_vs_CTRL = p.Boot3, p.adj_II = p.Boot3.adj)
write.csv(RES3, file=paste("minimal/fitness/", cond, ".fitlm_se_CI_diff_p_1379_gene_pairs.csv", sep=""))

#----------------------------------------------@@@@@@@@@@@@@@@@@@@---------------------------------------------------------
#==============================================PART 2: GI scoring=======================================================
#----------------------------------------------@@@@@@@@@@@@@@@@@@@---------------------------------------------------------
#--------------------------------------------I. Construct fitness matrix of SKO sgRNAs-------------------------------------
CTRL_CTRL = list()
totalcond = c("min.1", "min.2")
for (cond in totalcond) {
  read0 = read.csv(paste("minimal/fitness/",cond, ".RPM_fitlm.csv", sep=""), header=T, row.names=1)
  spli = strsplit(rownames(read0), split="-")
  # this is sgRNA_label1
  spli0 = unique(unlist(spli))  
  # this is gene only
  spli1 = unique(unlist(lapply(strsplit(spli0, split="_"), function(x) {x[[1]][1]}))) 
  #CTRL sgRNA
  spli0_CTRL = spli0[grepl("^CTRL", spli0)]
  #Non-CTRL sgRNA
  spli0_noCTRL = spli0[!grepl("^CTRL", spli0)]
  #CTRL gene
  spli1_CTRL = spli1[grepl("^CTRL", spli1)]
  #Non-CTRL gene
  spli1_noCTRL = spli1[!grepl("^CTRL", spli1)]
  
  allsgRNAs = sort(spli0_noCTRL)
  allCTRLs = sort(spli0_CTRL)
  wholeM = matrix(data=NA, nrow=length(allsgRNAs), ncol=length(allCTRLs))
  rownames(wholeM) = allsgRNAs
  colnames(wholeM) = allCTRLs
  for (sname in rownames(read0)) {
    splitName = strsplit(sname, split="-")[[1]]
    if ((sum(splitName %in% allsgRNAs)==1) & (sum(splitName %in% allCTRLs)==1)) {
      if (substr(splitName[1], 1, 4) == "CTRL") {
        wholeM[splitName[2], splitName[1]] = read0[sname, 8]
      } else {
        wholeM[splitName[1], splitName[2]] = read0[sname, 8]
      }
    }
  }
  write.csv(wholeM, file=paste("minimal/fitness/", cond, "_fitMatrix_SKO_sgRNAs.csv", sep=""))
}

#--------------------------------------------II. Functions for GI scoring-------------------------------------
# Function1: "differentModels" will use f_xy-function(f_x, f_y) to determine GGI score for each DKOs. 
differentModels = function(str_name, num_d, num_s1, num_s2) {
  n1 = length(num_s1)
  n2 = length(num_s2)
  totalName = NULL
  totalLine = NULL
  if ((n1==0) | (n2==0)) {
    # NA for all single knockout
    naout = matrix(c(NA, NA, NA, NA, NA, NA), nrow=1)
    colnames(naout) = c("Add", "Max", "Min", "f.DKO", "f.SKO1", "f.SKO2")
    return(naout)
  } else {
    for (i1 in 1:n1) {
      for (i2 in 1:n2) {
        name_line = paste(str_name, names(num_s1)[i1], names(num_s2)[i2], sep="-") 
        add_line = num_d - (num_s1[i1] + num_s2[i2]) 
        #multiple_line = num_d/num_ctrl - (num_s1[i1]/num_ctrl * num_s2[i2]/num_ctrl)
        max_line = num_d - max(num_s1[i1], num_s2[i2])
        min_line = num_d - min(num_s1[i1], num_s2[i2])
        all_line = c(add_line, max_line, min_line, num_d, num_s1[i1], num_s2[i2])
        totalName = c(totalName, name_line)
        totalLine = rbind(totalLine, all_line)
      }
    }
    rownames(totalLine) = totalName
    colnames(totalLine) = c("Add", "Max", "Min", "f.DKO", "f.SKO1", "f.SKO2")
    return(totalLine)
  }
}

#Function2: "computeGGI" organizes the inputs(sname, f.DKO, f.SKO1, f.SKO2) for "differentModels"
computeGGI = function(doubleK, singlerate=wholeM) {
  # doubleK will be one line of double-knockout construct in read0
  sname = as.character(as.numeric(doubleK[1]))
  sname = as.character(doubleK[1])
  sname = as.character(unlist(doubleK[1]))
  # print(sname)
  fitvalue = as.numeric(doubleK[length(doubleK)])
  sgname1 = strsplit(sname, split="-")[[1]][1]
  sgname2 = strsplit(sname, split="-")[[1]][2]
  sgname1_single = singlerate[sgname1, ][!is.na(singlerate[sgname1, ])]
  sgname2_single = singlerate[sgname2, ][!is.na(singlerate[sgname2, ])]
  ggiscore = differentModels(sname, fitvalue, sgname1_single, sgname2_single)
  return(data.frame(ggiscore))
}

#-------------------III. Execute the GI calculation for the two respective replicates-------------------------------------
# Inputs for sgRNA-level GI score calculation: fitness of DKO and SKO sgRNAs, respectively
cond = 'min.2'
read0 = read.csv(file = paste0("minimal/fitness/",cond,".RPM_fitlm.csv"), header=T, row.names=1)
wholeM = as.matrix(read.csv(file = paste0("minimal/fitness/",cond,"_fitMatrix_SKO_sgRNAs.csv"), header=T, row.names=1))

read0_noctrl = read0[!grepl("CTRL", rownames(read0)), ]
GGIvalues.0 = apply(cbind(rownames(read0_noctrl), read0_noctrl), 1, computeGGI)
GGIvalues.1 = do.call(rbind, GGIvalues.0)
GGIvalues = GGIvalues.1[!is.na(rowSums(GGIvalues.1)), ]
rownames(GGIvalues) = as.character(sapply(rownames(GGIvalues), function(x) as.character(strsplit(x, split="\\.")[[1]][length(strsplit(x, split="\\.")[[1]])])))
write.csv(GGIvalues, file=paste("minimal/GI/", cond, "_GI_sgRNA_level.csv",sep=""))

#---------IV. Calculate the gene-level GI scores (mean, sd, p.value and z.score) from the sgRNA-level GI scores-----------
# Function "MeanSD"
MeanSD = function(xx) {
  # this function will calculate mean across all CTRLs and sgRNAs, p-values of double knockouts vs.
  #single knockout+CTRL, and averaged double knockout fitness and single knockout+CTRL fitness.
  xx_mean = apply(xx[, 1:3], 2, mean, na.rm=T)
  xx_sd = apply(xx[, 1:3], 2, sd, na.rm=T)
  
  xx_calculateP = function(xx_row) {
    xx_double_value = xx_row[4]
    # add
    xx_add_value = xx_row[5] + xx_row[6]
    # max
    xx_max_value = max(xx_row[5], xx_row[6])
    # min
    xx_min_value = min(xx_row[5], xx_row[6])
    c(xx_double_value, xx_add_value, xx_max_value, xx_min_value)
  }
  
  xx_ttest = t(apply(xx, 1, xx_calculateP))
  xx_p_add = t.test(xx_ttest[, 1], xx_ttest[, 2], alternative = "greater")$p.value
  xx_p_max = t.test(xx_ttest[, 1], xx_ttest[, 3],alternative = "greater")$p.value
  xx_p_min = t.test(xx_ttest[, 1], xx_ttest[, 4],alternative = "less")$p.value
  
  xx_avg_double = mean(unique(xx[, 4]), na.rm=T)
  xx_avg_single1 = mean(unique(xx[, 5]), na.rm=T)
  xx_avg_single2 = mean(unique(xx[, 6]), na.rm=T)
  
  c(rbind(xx_mean, xx_sd), xx_p_add, xx_p_max, xx_p_min, xx_avg_double, xx_avg_single1, xx_avg_single2)
}

# Execute the Function "MeanSD", and the input 
totalGenes = gene52
Num = length(gene52)
total_outline = NULL
total_combo = NULL

for (idx1 in 1:(Num-1)) {
  for (idx2 in (idx1+1):Num) {
    g1 = totalGenes[idx1]
    g2 = totalGenes[idx2]
    # cat("working on", g1, g2, "\n")
    geneflag = grepl(g1, names(GGIvalues.0)) & grepl(g2, names(GGIvalues.0))
    sub_GGIvalues.0 = GGIvalues.0[geneflag]
    sub_GGIvalues.3 = do.call(rbind, sub_GGIvalues.0)
    outline = MeanSD(sub_GGIvalues.3)
    total_outline = rbind(total_outline, outline)
    total_combo = c(total_combo, paste(g1, g2, sep="-"))
  }
}
colnames(total_outline) = c("Add.mean", "Add.sd", "Max.mean", "Max.sd", "Min.mean", "Min.sd", "p.value_Add", "p.value_Max", "p.value_Min", "avg.f.DKO", "avg.f.SKO1", "avg.f.SKO2")
rownames(total_outline) = total_combo
z_Add <- (total_outline[,1] - mean(total_outline[,7]))/sd(total_outline[,7])
z_Max <- (total_outline[,3] - mean(total_outline[,8]))/sd(total_outline[,8])
z_Min <- (total_outline[,5] - mean(total_outline[,9]))/sd(total_outline[,9])
p.adj_Add <- p.adjust(total_outline[,7],method = "BH")
p.adj_Max <- p.adjust(total_outline[,8],method = "BH")
p.adj_Min <- p.adjust(total_outline[,9],method = "BH")
total <- cbind(total_outline, z_Add, z_Max, z_Min, p.adj_Add, p.adj_Max, p.adj_Min)
write.csv(total, file=paste("minimal/GI/", cond, "_GI_gene_level.csv",sep=""))

#----------------------------------------------@@@@@@@@@@@@@@@@@@@---------------------------------------------------------
#================================PART 3: Combine fitness and GI results of two replicates=================================================
#----------------------------------------------@@@@@@@@@@@@@@@@@@@---------------------------------------------------------
#Function to alphabetically order the two gene names in each gene pair
Path = "~/Downloads/data/"
Order_rname = function(file_name){
  file = read.csv(paste0(Path,file_name),  header = T, row.names = 1)
  
  f1 = unlist(lapply(rownames(file), function(x) {strsplit(x, split = "-")[[1]][1]}))
  f2 = unlist(lapply(rownames(file), function(x) {strsplit(x, split = "-")[[1]][2]}))
  f = cbind(f1, f2)
  for(i in 1:nrow(f)){
    if( sum(f[i,] %in% "CTRL")==0 ){
      if(f[i,1] > f[i,2]){
        p = f[i,1]
        f[i,1] = f[i,2]
        f[i,2] = p
      }
    }
    else if(sum(f[i,] %in% "CTRL") == 1){
      f[i,1] = f[i][which(f[i,] != "CTRL" )]
      f[i,2] = "CTRL"
    }
  }
  new.name = paste(f[,1], f[,2], sep = "-")
  rownames(file) = new.name
  return (file)
}

#fitness file: Change rownames by alphabetic order
Fit = Order_rname(file_name = "minimal/fitness/min.combined.fitlm_se_CI_diff_p_1379_gene_pairs.csv")
dim(Fit)                                   
write.csv(Fit, "minimal/fitness/final_min_combined_fitlm_se_CI_diff_p_1379_gene_pairs_ordered.csv")

#GI score file: Change rownames by alphabetic order
GI_min.1 = Order_rname(file_name = "minimal/GI/min.1_GI_gene_level.csv")
colnm = colnames(GI_min.1)
colnames(GI_min.1) = paste0("Rep1_", colnm)
dim(GI_min.1)
GI_min.1 = GI_min.1[order(rownames(GI_min.1)),]
write.csv(GI_min.1, "minimal/GI/final_min.1_GI_gene_level_ordered.csv")

GI_min.2 = Order_rname(file_name = "minimal/GI/min.2_GI_gene_level.csv")
colnm = colnames(GI_min.2)
colnames(GI_min.2) = paste0("Rep2_", colnm)
dim(GI_min.2)
GI_min.2 = GI_min.2[order(rownames(GI_min.2)),]
write.csv(GI_min.2, "minimal/GI/final_min.2_GI.gene_level_ordered.csv")

#Combine fitness and GI score
all(rownames(GI_min.1) == rownames(GI_min.2))
Res = cbind.data.frame(Fit[rownames(GI_min.1),],GI_min.1, GI_min.2)

##
N.add1 = which(colnames(Res) == "Rep1_Add.mean") # N.add1 = 11
N.add2 = which(colnames(Res) == "Rep2_Add.mean") # N.add2 = 29
N.max1 = which(colnames(Res) == "Rep1_Max.mean") # N.max1 = 13
N.max2 = which(colnames(Res) == "Rep2_Max.mean") # N.max2 = 31
N.min1 = which(colnames(Res) == "Rep1_Min.mean") # N.min1 = 15
N.min2 = which(colnames(Res) == "Rep2_Min.mean") # N.min2 = 33

##
P.add1 = which(colnames(Res) == "Rep1_p.value_Add") # P.add1 = 17
P.add2 = which(colnames(Res) == "Rep2_p.value_Add") # P.add2 = 35
P.max1 = which(colnames(Res) == "Rep1_p.value_Max") # P.max1 = 18
P.max2 = which(colnames(Res) == "Rep2_p.value_Max") # P.max2 = 36
P.min1 = which(colnames(Res) == "Rep1_p.value_Min") # P.min1 = 19
P.min2 = which(colnames(Res) == "Rep2_p.value_Min") # P.min2 = 37

library(metap)
#Inverse variance weighted average = (mean1/var1 + mean2/var2)/(1/var1 + 1/var2)
Res$gi.Weighted.Avg.Add = apply(Res, 1, function(x) {(x[N.add1]/(x[N.add1+1]^2) + x[N.add2]/(x[N.add2+1]^2))/(1/(x[N.add1+1]^2) + 1/(x[N.add2+1]^2)) })
Res$gi.Weighted.Avg.Max = apply(Res, 1, function(x) {(x[N.max1]/(x[N.max1+1]^2) + x[N.max2]/(x[N.max2+1]^2))/(1/(x[N.max1+1]^2) + 1/(x[N.max2+1]^2)) })
Res$gi.Weighted.Avg.Min = apply(Res, 1, function(x) {(x[N.min1]/(x[N.min1+1]^2) + x[N.min2]/(x[N.min2+1]^2))/(1/(x[N.min1+1]^2) + 1/(x[N.min2+1]^2)) })

#Combine p value by the sum of logs (Fisher's) method. sumz(Stoufer's method) refuse to take p value equal 1
Res$fisher.p.Add = apply(Res, 1, function(x) sumlog(as.numeric(c(x[P.add1], x[P.add2])))[[3]])
Res$fisher.p.Add.adj = p.adjust(Res$fisher.p.Add, method = "BH")
Res$fisher.p.Max = apply(Res, 1, function(x) sumlog(as.numeric(c(x[P.max1], x[P.max2])))[[3]])
Res$fisher.p.Max.adj = p.adjust(Res$fisher.p.Max, method = "BH")
Res$fisher.p.Min = apply(Res, 1, function(x) sumlog(as.numeric(c(x[P.min1], x[P.min2])))[[3]])
Res$fisher.p.Min.adj = p.adjust(Res$fisher.p.Min, method = "BH")

write.csv(Res, "minimal/Final_min_f.avg_GI.csv")

#----------------------------------------------@@@@@@@@@@@@@@@@@@@---------------------------------------------------------
#=================================PART 4: Plots and Figures based on results above================================================
#----------------------------------------------@@@@@@@@@@@@@@@@@@@---------------------------------------------------------
#1. Barplot of 52 SKOs and CTRL
#2. Correlation of gene-level fitness of two replicates
#3. Plotting oncogenic GIs
#4. Correlation of gene-level weighted GI scores of two replicates
#5.GI profiles similarity heatmap
#----------------------------------I. Barplot of 53 single gene fitness--------------------------------
fit = read.csv("minimal/fitness/final_min_combined_fitlm_se_CI_diff_p_1379_gene_pairs_ordered.csv", row.names = 1, header = T)
dat = data.frame(fit[grepl("CTRL", rownames(fit)),])
dat$sname = unlist(lapply(rownames(dat), function(x) strsplit(x, split = "-")[[1]][1]))
write.csv(dat, "minimal/fitness/min_SKO_fitness_gene_level.csv")
positions <- dat$sname[order(dat$f.avg, decreasing = F)]

p <- ggplot(dat,aes(x=sname,y=f.avg))+
  geom_bar(position=position_dodge(), stat="identity",fill = "blue")+
  geom_errorbar(aes(ymin=l, ymax=u), width=.2)+
  scale_x_discrete(limits = positions)+
  theme(#panel.border = element_blank(), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line = element_line(colour = "black"),
    #legend.position = c(0.7,0.7),
    #legend.title = element_text(size = 20, face= "bold"),
    #legend.text = element_text(size=18 ),
    axis.text.y  =element_text(size=12, color = "black"),
    axis.text.x =element_text(size=19, color = "black"),
    axis.title=element_text(size=20, color = "black"),#face="bold")
    plot.title = element_text(size=21,hjust=0.5, color = "black")
  )+
  scale_y_continuous(limits = c(-0.15,0.22), breaks = round(seq(-0.1,0.2,0.1),2))+
  coord_flip()+
  labs(title="Minimal",x ="SKOs", y = "Fitness")
p
ggsave("minimal/fitness/Min_SKO_fitness_53_genes.pdf", width = 3.6, height= 9)

#----------------------------------II. Correlation of gene-level fitness of two replicates--------------------------------
f1 = read.csv("minimal/fitness/min.1.fitlm_se_CI_diff_p_1379_gene_pairs.csv", row.names = 1, header = T)
f2 = read.csv("minimal/fitness/min.2.fitlm_se_CI_diff_p_1379_gene_pairs.csv", row.names = 1, header = T)

df = data.frame(rep1 = f1[,"f.avg"], rep2 = f2[rownames(f1), "f.avg"])
rownames(df) = rownames(f1)
ggscatter(df, x = "rep1" , y = "rep2", size = 3, alpha = 0.5,
          add = "reg.line", conf.int = FALSE, 
          add.params = list(color = "#007AFF", fill = "lightgray"), # Customize reg. line
          cor.coef = TRUE,rug = TRUE,
          cor.coeff.args = list(method = "pearson", label.x = -0.17,label.y = 0.25, color = "blue",label.sep = "\n", size = 7))+
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
    axis.text.y  =element_text(size=24,color = 'black'),
    axis.text.x =element_text(size=24,color = 'black', vjust = 0.5,margin = margin(5,10,0,0)),
    axis.title=element_text(size=26,margin = margin(10,0,0,0)),#face="bold")
    plot.title = element_text(size=26,hjust= 0.5, margin = margin(0,0,10,0))
    #strip.text  = element_text(size = 16, face="bold")
    ##strip.text is for facet
  )+
  scale_x_continuous(limits = c(-0.2,0.3), breaks = round(seq(-0.2,0.3,0.15),2))+
  scale_y_continuous(limits = c(-0.2,0.3), breaks = round(seq(-0.2,0.3,0.1),2))+
  geom_hline(yintercept= 0, linetype="dashed", color = "grey")+
  geom_vline(xintercept= 0, linetype="dashed", color = "grey")+
  labs(title="Minimal",x ="Fitness, Rep 1", y = "Fitness, Rep 2")

ggsave("minimal/fitness/Reproducibility_gene_level_fitnesss_min.pdf", width = 4.4, height = 4)

#---------------------------------- III. Plotting oncogenic GIs--------------------------------
res=  data.frame(read.csv("minimal/Final_min_f.avg_GI.csv", header=T, row.names=1))
res$sname = rownames(res)
res$oncogenic = "Not_oncogenic"
res[res$f.avg > 0  & res$p.adj_I < 0.05, "oncogenic"] = "Yes_oncogenic"

res$Sig = "Not_Sig"
res[res$f.avg > 0 & res$p.adj_I < 0.05 & res$gi.Weighted.Avg.Max >0 & res$fisher.p.Max.adj < 0.05, "Sig"] = "Yes_Sig"
#res[res$f.avg > f.ctrl & res$p.adj_II < 0.05 & res$gi.Weighted.Avg.Max >0 & res$fisher.p.Max.adj < 0.05, "Sig"] = "Yes_Sig"

sel = res[res$f.avg > 0 & res$p.adj_I < 0.05 & res$gi.Weighted.Avg.Max >0 & res$fisher.p.Max.adj < 0.05, ]     
dim(sel)
gene1 <- sapply(rownames(sel), function(x){unlist(strsplit(x,"-")[[1]][1])} )
gene2 <- sapply(rownames(sel), function(x){unlist(strsplit(x,"-")[[1]][2])} )

PPI_net =  cbind.data.frame(gene1, gene2, ggi = sel[,"gi.Weighted.Avg.Max"])
write.csv(PPI_net, "minimal/GI/Oncogenic_GIs_cytoscape_min.csv", quote = FALSE)

p <- ggplot(res, aes(x = gi.Weighted.Avg.Max, y = gi.Weighted.Avg.Min)) +
  geom_point(aes(size=Sig,color=oncogenic,group = Sig, alpha = Sig)) +    # Use hollow circles
  scale_size_manual(values = c(2,4)) +
  scale_color_manual(values= c("grey","darkorchid1"))+
  scale_alpha_manual(values = c(0.2, 1))+
  #geom_smooth(method=lm,linetype = "dashed", color= "black")+ # Add linear regression line #  (by default includes 95% confidence region) 
  theme(panel.border = element_rect(colour = "black",size = 0.6,fill=NA), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.text.y  =element_text(size=24),
        axis.text.x =element_text(size=24,margin = margin(t = 1, r = 1, b = 1, l = 1)),
        axis.title=element_text(size=24, hjust =0.5),#face="bold")
        plot.title = element_text(size=24,hjust=0.5)
        #strip.text.y  = element_text(size = 12, face="bold")
        ##strip.text is for facet
  )+
  labs(title="Growth-promoting interactions\n(Minimal)",x ="Mean: f(DKO) - f(SKO1)", y = "Mean: f(DKO) - f(SKO2)")+
  scale_x_continuous(limits =c(-0.15,0.12), breaks = round(seq(-0.15, 0.1, 0.05),2)) +
  scale_y_continuous(limits= c(-0.1,0.35), breaks = round(seq(-0.1, 0.4, 0.1),2))+
  geom_hline(yintercept= 0, linetype="dashed", color = "grey")+
  geom_vline(xintercept= 0, linetype="dashed", color = "grey")
p
require("ggrepel")
set.seed(42)
data = subset(res, Sig == "Yes_Sig")
d0 = data[with(data, order(-f.avg,-gi.Weighted.Avg.Max)),]
p+ geom_text_repel(
  data = head(d0,5), 
  aes(label = sname),
  size = 6,
  color = "darkblue",
  box.padding = unit(0.4, "lines"),
  point.padding = unit(0.3, "lines"),
  nudge_y =0.005,
  nudge_x = 0.06
)

ggsave("minimal/GI/GI_scatter_plot_min.pdf", width = 6.5, height = 7)

#---------------------------------- IV. Correlation of gene-level weighted GI scores of two replicates--------------------------------
res=  data.frame(read.csv("minimal/Final_min_f.avg_GI.csv", header=T, row.names=1))
ggscatter(res, x = "Rep1_Max.mean" , y = "Rep2_Max.mean", size = 3, 
          add = "reg.line", conf.int = FALSE,  alpha = 0.5,
          add.params = list(color = "#007AFF", fill = "lightgray"), # Customize reg. line
          cor.coef = TRUE,rug = TRUE,
          cor.coeff.args = list(method = "pearson", label.x = -0.14,label.y = 0.075, color = "blue",label.sep = "\n", size = 7))+
  theme(#panel.border = element_blank(), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    #axis.ticks.y = element_line(colour = "black"),
    legend.position = 'none',
    #legend.position = c(0.75,0.85),
    #legend.title = element_text(size = 20, face= "bold"),
    legend.text = element_text(size=20),
    axis.text.y  =element_text(size=24,color = 'black'),
    axis.text.x =element_text(size=24,color = 'black', vjust = 0.5,margin = margin(5,10,0,0)),
    axis.title=element_text(size=26,margin = margin(10,0,0,0)),#face="bold")
    plot.title = element_text(size=26,hjust= 0.5, margin = margin(0,0,10,0))
    #strip.text  = element_text(size = 16, face="bold")
    ##strip.text is for facet
  )+
  scale_x_continuous(limits = c(-0.15,0.1), breaks = round(seq(-0.2,0.1,0.1),2))+
  scale_y_continuous(limits = c(-0.15,0.1), breaks = round(seq(-0.2,0.1,0.05),2))+
  geom_hline(yintercept= 0, linetype="dashed", color = "grey")+
  geom_vline(xintercept= 0, linetype="dashed", color = "grey")+
  labs(title="Minimal",x ="GI score, Rep 1", y = "GI score, Rep 2")

ggsave("minimal/GI/Reproducibility_gene_level_GI_score_min.pdf", width = 4.4, height = 4)

#---------------------------------- V. GI profiles similarity heatmap--------------------------------
#Color brewer
numericToColor = function(x, mycol = c("purple", "grey80", "green"), mybreaks = seq(-0.15, 0.15, 0.01), bkground = "grey80") {
  # x is a list of numeric values. It will assign colors based on given mybreaks points.
  totalN = length(mybreaks) + 1
  allColors = colorRampPalette(mycol)(totalN)
  output = rep(bkground, length(x))
  minP = min(mybreaks)
  maxP = max(mybreaks)
  output[x <= minP] = allColors[1]
  for (index in 2:totalN-1) {
    output[(x > mybreaks[index]) & (x <= mybreaks[index+1])] = allColors[index]
  }
  output[x >= maxP] = allColors[totalN]
  output
}

##Input matrix
input_s = read.csv("minimal/Final_min_f.avg_GI.csv", row.names=1)
name_p = data.frame(do.call(rbind, strsplit(rownames(input_s), split="-")))
colnames(name_p) = c("gene1", "gene2")
name_p2 = name_p[, c("gene2", "gene1")]
colnames(name_p2) = c("gene1", "gene2")
input1 = cbind(name_p, input_s)
input2 = cbind(name_p2, input_s)
rownames(input2) = paste(rownames(input2), "-2", sep="")
input = rbind(input1, input2)

##Cast long-format dataframe into short-format dataframe
totalvalues = c("gi.Weighted.Avg.Max", "gi.Weighted.Avg.Add")
for (svalue in totalvalues) {
  M_p = cast(input, gene1~gene2, value=svalue)
  rownames(M_p) = as.character(M_p[, 1])
  M_p2 = M_p[, 2:dim(M_p)[2]]
  neworder = rownames(M_p)[order(rownames(M_p))]
  M_p_ordered = M_p2[neworder, neworder]
  diag(M_p_ordered) = 0
  GImatrix = data.frame(M_p_ordered)
  corrMatrix = cor(M_p_ordered, method = "pearson")
  colnames(GImatrix) = neworder
  colnames(corrMatrix) = neworder
  rownames(GImatrix) = neworder
  rownames(corrMatrix) = neworder
  write.csv(GImatrix, file=paste("minimal/GI/", svalue, ".GIMatrix.csv",sep=""))
  write.csv(corrMatrix, file=paste("minimal/GI/", svalue, ".correlation.Based.On.GI.csv",sep=""))
  
  ##Determine color of sinlge-gene fitness
  f0 = data.frame(read.csv("minimal/fitness/min_SKO_fitness_gene_level.csv", row.names = 1, header = T))
  f = f0[which(f0$sname != "CTRL"),]
  f = f[order(f$sname, decreasing = F),]
  range(f$f.avg)
  f.color = numericToColor(as.numeric(f$f.avg),mycol = c("purple", "grey80", "green"), mybreaks = seq(-0.1, 0.1, 0.001), bkground = "grey80" )
  
  #GI map
  # png(filename = paste(stype, ".", svalue, ".GIMatrix.png",sep=""),width=1000, height=1000)
  # heatmap.2(as.matrix(GImatrix),
  #            Rowv = TRUE,
  #            Colv = TRUE,
  #            hclustfun = function(x) hclust(x,method = "complete"),
  #            dendrogram="both",
  #            cexRow = 1.4, cexCol = 1.4,
  #            scale = "none",
  #            key=TRUE, keysize=1, symkey=FALSE,
  #            density.info="none", trace="none",
  #            col=colorRampPalette(c("cyan", "black", "yellow")),
  #            breaks = seq(-0.07,0.07, 0.005),
  #            margins=c(10,10)
  #  )
  #  dev.off()
  
  #GI profiling 
  png(filename = paste("minimal/GI/", svalue, ".correlation.Based.On.GI.png",sep=""),width=1000, height=1000)
  heatmap.2(as.matrix(corrMatrix),
            Rowv = TRUE, 
            Colv = TRUE,
            hclustfun = function(x) hclust(x,method = "complete"),
            dendrogram="both",
            cexRow = 1.8, cexCol = 1.8,
            scale = "none", 
            key=TRUE, keysize=0.5, symkey=TRUE,
            density.info="none", trace="none",
            col=colorRampPalette(c("deepskyblue", "black", "yellow")),
            breaks = seq(-0.3,1,0.05),
            margins=c(8,8),
            RowSideColors = f.color,
            colRow = rep("black",52)
            
  )
  dev.off()
}

range(corrMatrix)
range(f$f.avg)

#Scale color for GI correlation coefficient r
##Determine colors for heatmap
breaks = c(seq(-0.5,1,0.01))
my.colors1 <- colorRampPalette(c("deepskyblue", "black"))(length(seq(-0.5, -0.01, 0.01)))
my.colors2 <- colorRampPalette(c("black", "yellow"))(length(seq(0.01,1,0.01)))
my.colors = c(my.colors1, my.colors2)  
df <- reshape2::melt(outer(1:2,1:2), varnames = c("X1", "X2"))
ggplot(df, aes(X1, X2, fill = value)) + 
  geom_raster()+ 
  scale_fill_gradientn(name = "GI profile similarity score",
                       colours=colorRampPalette(c("deepskyblue", "black", "yellow"))(length(seq(-0.3,1,0.05))),na.value = "transparent",
                       breaks=seq(0,1,0.5),labels=seq(0,1,0.5),
                       limits=c(-0.3,1))+
  guides(size = 3, fill = guide_colourbar(ticks = TRUE))+
  theme(legend.text = element_text(size=14),
        legend.title = element_text(size = 14))
##barwidth = 0.5, barheight = 10,label.theme = element_text(colour = "blue", angle = 0))))
ggsave("minimal/GI/min_scale_color_GI_correlation_r.pdf", width = 3, height =3)

#Scale color for fitness 
df <- reshape2::melt(outer(1:2,1:2), varnames = c("X1", "X2"))
ggplot(df, aes(X1, X2, fill = value)) + 
  geom_raster()+ 
  scale_fill_gradientn(name = "Fitness",
                       colours=colorRampPalette(c("purple", "grey80", "green"))(length(seq(-0.1, 0.1, 0.001))), na.value = "transparent",
                       breaks=seq(-0.1,0.1,0.1),labels=seq(-0.1,0.1,0.1),
                       limits=c(-0.1,0.1))+
  guides(size = 3, fill = guide_colourbar(ticks = TRUE))+
  theme(legend.text = element_text(size=14),
        legend.title = element_text(size = 14))
##barwidth = 0.5, barheight = 10,label.theme = element_text(colour = "blue", angle = 0))))
ggsave("minimal/GI/min_scale_color_single_gene_fitness.pdf", width = 3, height =3)






