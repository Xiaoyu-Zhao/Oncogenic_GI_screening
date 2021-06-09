#Xiaoyu Zhao
library(ggplot2)
library(RColorBrewer)
library(Cairo)
setwd("~/Downloads/data/")
ref <- data.frame(read_excel("TableS1_Annotatation_of_asssyed_TSGs_or_TSG_combinations.xlsx", 5))
gene52 = as.character(unique(res[,1]))

#---------------------------------I. Measure fitness = lm(log2RPM ~ time) without any filters------------------------------------------
#Minimal medium
#Input: sheet3 in TableS5 ---> min.1; sheet4 in TableS5 --->min.2
cond = "min.2"
read = data.frame(read_excel("TableS5_Raw_and_processed_data_combinatorial_screening_in_vitro.xlsx",4))
rownames(read) = read[,"gene_names"]
dat00 = as.matrix(read[,3:ncol(read)])
# 1. filter some bad sgRNAs
filtered_list_others = "Random"
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
write.csv(cbind(rpm, fitlm), file=paste("minimal/fitness/",cond, ".RPM_fitlm_no_sgRNA_filter.csv", sep=""))

#TGFbeta medium
#Input: sheet3 in TableS5 ---> min.1; sheet4 in TableS5 --->min.2
cond = "tgf.2"
read = data.frame(read_excel("TableS5_Raw_and_processed_data_combinatorial_screening_in_vitro.xlsx",5))
rownames(read) = read[,"gene_names"]
dat00 = as.matrix(read[,3:ncol(read)])
# 1. filter some bad sgRNAs
filtered_list_others = "Random"
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
x = c(0, 4, 8, 12, 16, 20)
fitgrowth_lm = function(y) {
  outlm = coef(lm(y ~ x))
}
fitlm = t(apply(logrpm, 1, fitgrowth_lm))
#Save results
write.csv(cbind(rpm, fitlm), file=paste("tgf/fitness/",cond, ".RPM_fitlm_no_sgRNA_filter.csv", sep=""))

#-------------------------------------------II. sgRNA density plots------------------------------------------
#Function"sgRNA_d_plot is to get fitness of all the dual-sgRNAs containing the target sgRNA
sgRNA_d_plot = function(medium, fitness = f, totalGenes = gene52, ref = ref, legend.position = c(0.84,0.75)){
  fitness$sgRNAs = NA
  for(i in 1:length(totalGenes)){
    f_sub = fitness[grepl(totalGenes[i], rownames(fitness)),]
    sgs_sub = as.character(ref[grepl(totalGenes[i],ref[,1]),2])
    
    for(j in 1:length(sgs_sub)){
      f_sub[grepl(sgs_sub[j],rownames(f_sub)), "sgRNAs"] = sgs_sub[j]
    }
    p <- ggplot(f_sub, aes(x = f.avg, color = sgRNAs))+
      geom_density(alpha=.4, position="identity", size = 2)+
      scale_color_manual(values= c("red", "orange", "blue", "green", "purple", "yellow", "magenta"))+
      theme(#panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        axis.line = element_line(colour = "black"),
        legend.position = legend.position, #0.855 for sgNF2, 
        legend.title = element_text(size = 14),
        legend.text = element_text(size=14 ),
        axis.text.y  =element_text(size=28,colour = "black"),
        axis.text.x =element_text(size=28,colour = "black"),
        axis.title=element_text(size=30),#face="bold")
        plot.title = element_text(size=28,hjust=0.5)
        #strip.text.y  = element_text(size = 26, face="bold")
        ##strip.text is for facet
      )+
      labs(title= paste0("sg", totalGenes[i]), x ="sgRNA level fitness ", y = "Density")
    ggsave(file = paste0(medium, "/sgRNA_density_plots/", totalGenes[i],".pdf"), width = 5, height = 5 )
  } 
}
#sgRNA density plots for minimal medium
f.min1 = read.csv("minimal/fitness/min.1.RPM_fitlm_no_sgRNA_filter.csv", row.names = 1)
dim(f.min1)
f.min2 = read.csv("minimal/fitness/min.2.RPM_fitlm_no_sgRNA_filter.csv", row.names = 1)
dim(f.min2)
comM = intersect(rownames(f.min1), rownames(f.min2))
f.minimal = data.frame(f1 = f.min1[comM,"x"], f2 = f.min2[comM,"x"])
f.minimal$f.avg = (f.minimal$f1 + f.minimal$f2)/2
rownames(f.minimal) = comM
##legend position for sgNF2 was adjusted to c(0.855 , 0.75)
sgRNA_d_plot(medium = "minimal", fitness = f.minimal, totalGenes = gene52, ref = ref, legend.position = c(0.84,0.75))

#--------------------------------sgRNA density plots for TGFbeta1 medium---------------------------------------
f.tgf1 = read.csv("tgf/fitness/tgf.1.RPM_fitlm_no_sgRNA_filter.csv", row.names = 1)
dim(f.tgf1)
f.tgf2 = read.csv("tgf/fitness/tgf.2.RPM_fitlm_no_sgRNA_filter.csv", row.names = 1)
dim(f.tgf2)
comM = intersect(rownames(f.tgf1), rownames(f.tgf2))
f.tgf = data.frame(f1 = f.tgf1[comM,"x"], f2 = f.tgf2[comM,"x"])
f.tgf$f.avg = (f.tgf$f1 + f.tgf$f2)/2
rownames(f.tgf) = comM
##legend position for sgSMAD4 was adjusted to c(0.2 , 0.75)
sgRNA_d_plot(medium = "tgf", fitness = f.tgf, totalGenes = gene52, ref = ref, legend.position = c(0.8,0.75))

