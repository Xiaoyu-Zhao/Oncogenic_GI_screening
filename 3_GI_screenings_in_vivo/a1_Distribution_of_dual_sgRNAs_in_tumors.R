#Xiaoyu Zhao
library(ggplot2)
library(readxl)
setwd("~/Downloads/data/")
#-----------------------------I. CDF of dual-sgRNAs in plasmid and tumors-----------------------------
ct.plasmid <- data.frame(read_excel("TableS2_Evaluation_of_TSG_plasmid_cell_libraries.xlsx", 1))
rownames(ct.plasmid) = ct.plasmid$gene_names
ct.pten <- data.frame(read_excel("TableS3_Raw_and_processed_data_combinatorial_screening_in_vivo.xlsx", 1))
rownames(ct.pten) = ct.pten$gene_names
ct.pik <- data.frame(read_excel("TableS3_Raw_and_processed_data_combinatorial_screening_in_vivo.xlsx", 2))
rownames(ct.pik) = ct.pik$gene_names
ct.myc <- data.frame(read_excel("TableS3_Raw_and_processed_data_combinatorial_screening_in_vivo.xlsx",3))
rownames(ct.myc) = ct.myc$gene_names

##Combine counts in plasmids and all the tumor samples
cts = as.matrix(cbind(ct.plasmid[rownames(ct.plasmid),3:4], #n_sample = 2
            ct.pten[rownames(ct.plasmid),5:ncol(ct.pten)],  #n_sample = 9
            ct.pik[rownames(ct.plasmid),5:ncol(ct.pik)],    #n_sample = 3
            ct.myc[rownames(ct.plasmid),5:ncol(ct.myc)]) )  #n_sample = 6

#Remove all the dual-sgRNAs containing 'CTRL0003'
bad_candidate = c("CTRL0003")
flag = sapply(bad_candidate, function(x) {grepl(x, rownames(cts))})
counts = cts[!flag, ]
dim(counts)

#Calculate rpm(reads per million)
rpm0 = apply(counts, 2, function(x) {(x/sum(x)*1000000)})
#log2(mean_rpm), "+1" is to avoid log2(0)
logrpm0 = log2(rpm0+1)

#Mean rmp of two replicates
meanP = apply(logrpm0[,1:2], 1, function(x) mean(x))
meanPten = apply(logrpm0[,3:11], 1, function(x) mean(x))
meanPik  = apply(logrpm0[,12:14], 1, function(x) mean(x)) 
meanMyc = apply(logrpm0[,15:20], 1, function(x) mean(x))
log.meanrpm = cbind.data.frame(logrpm0, meanP, meanPten, meanPik, meanMyc)
#----------------------------------------------------------------------------------------------
#Comparison of CDF of dual-sgRNAs in plasmid and tumors
ks.test(log.meanrpm$meanP, log.meanrpm$meanPten) #D = 0.91597, p-value < 2.2e-16
ks.test(log.meanrpm$meanP, log.meanrpm$meanPik)  #D = 0.91597, p-value < 2.2e-16
ks.test(log.meanrpm$meanP, log.meanrpm$meanMyc)  #D = 0.91597, p-value < 2.2e-16
#----------------------------------------------------------------------------------------------
pdf("Cumulative_distribution_plasmid_and_tumors.pdf", width = 7, height = 5.2)
par(mar = c(5, 5, 2, 1))
plot(ecdf(log.meanrpm$meanP),col="red",xlim=c(0,10),verticals=T,do.points=F,lwd=4,col.01line = NULL,
     main = NULL,
     ylab= "Cumulative Distribution",
     xlab = expression(log[2]~RPM~per~dual_sgRNA),
     cex.lab=2, 
     cex.axis=1.8)
lines(ecdf(log.meanrpm$meanPten), col = "magenta",verticals=T,do.points=F,lwd=4,col.01line = NULL)
lines(ecdf(log.meanrpm$meanPik), col = "orange",verticals=T,do.points=F,lwd=4,col.01line = NULL)
lines(ecdf(log.meanrpm$meanMyc), col = "blue",verticals=T,do.points=F,lwd=4,col.01line = NULL)
#abline(v = 1, col = 'black', lty = 3, lwd = 2)
legend(5.3,0.5,c("Plasmid library",
                 "Tumor (PTEN-/-)", 
                 "Tumor (PIK3CA)",
                 "Tumor (MYC)"), 
fill = c("red","magenta","orange","blue"),  
bty = "n", 
pt.cex = 2, 
cex = 1.5, 
text.col = "black")
dev.off()

#----------------------II.Percentage of dual-sgRNAs retained in tumors-----------------------------
#Number of sgRNAs detected in plasmids
keep <- rowSums(counts[,1:2]) > 0
countsM = counts[keep,]
dim(countsM)
#Percentage of dual-sgRNAs retained in tumors
num_retained = apply(countsM[,3:ncol(countsM)], 2, function(x) sum(x > 0)/nrow(countsM)*100)
num_retained
#Results:
#F911.3    F912.2    F912.3    F913.2    F913.4    F912.1    F912.4    F913.1    F913.3    Cell.1    Cell.2    F995.1    F995.2 
#28.926534 21.975783 21.577017 14.573874 15.348120 11.616603 19.853883 12.722669  9.657702 97.770404 96.917569  8.604028 12.222028 
#F996.1    Cell.1    Cell.2    F908.1    F908.3    F908.4    F905.2  F905.4.1  F905.4.2 
#3.635464 97.319246 98.131331 12.373385 67.408895 41.363954 37.708115 21.515892  6.167773 

#Average retaining percentage in 9 of PTEN-/- tumors
mean(num_retained[1:9])  #17.2%
sd(num_retained[1:9])
#Average retaining percentage in 3 of PIK3CA tumors
mean(num_retained[10:12])  #8.093766%
sd(num_retained[10:12])
#Average retaining percentage in 6 of HAMyc tumors
mean(num_retained[13:18])  #30.86158%
sd(num_retained[13:18])

#--------------------------------III.Percentage of log2RPM >= 1 -----------------------------
ct.pten <- data.frame(read_excel("TableS3_Raw_and_processed_data_combinatorial_screening_in_vivo.xlsx", 1))
rownames(ct.pten) = ct.pten$gene_names
ct.pik <- data.frame(read_excel("TableS3_Raw_and_processed_data_combinatorial_screening_in_vivo.xlsx", 2))
rownames(ct.pik) = ct.pik$gene_names
ct.myc <- data.frame(read_excel("TableS3_Raw_and_processed_data_combinatorial_screening_in_vivo.xlsx",3))
rownames(ct.myc) = ct.myc$gene_names

##Combine counts in plasmids and all the tumor samples
cts.CP = as.matrix(cbind(ct.pten[rownames(ct.pten),3:ncol(ct.pten)],  #n_sample = 9
                      ct.pik[rownames(ct.pten),3:ncol(ct.pik)],    #n_sample = 3
                      ct.myc[rownames(ct.pten),3:ncol(ct.myc)]) )  #n_sample = 6

#Remove all the dual-sgRNAs containing 'CTRL0003'
bad_candidate = c("CTRL0003")
flag = sapply(bad_candidate, function(x) {grepl(x, rownames(cts.CP))})
counts.CP = cts.CP[!flag, ]
dim(counts.CP)
#Calculate rpm(reads per million)
rpm0.CP = apply(counts.CP, 2, function(x) {(x/sum(x)*1000000)})
#log2(mean_rpm), "+1" is to avoid log2(0)
logrpm0.CP = log2(rpm0.CP+1)
per = apply(logrpm0.CP, 2, function(x) sum(x >= 1)/nrow(logrpm0.CP)*100)
per

df1 = data.frame(sample = c('cell', 'tumor'),
                 mean = c(mean(c(per[1:2])), mean(c(per[3:11]))),
                 se = c(sd(c(per[1:2]))/sqrt(2), sd(c(per[3:11]))/sqrt(9)),
                 tumor = c('PTEN-/-', 'PTEN-/-'))

df2 = data.frame(sample = c('cell', 'tumor'),
                 mean = c(mean(c(per[12:13])), mean(c(per[14:16]))),
                 se = c(sd(c(per[12:13]))/sqrt(2), sd(c(per[14:16]))/sqrt(3)),
                 tumor = rep('PIK3CA', 2))

df3 = data.frame(sample = c('cell', 'tumor'),
                 mean = c(mean(c(per[17:18])), mean(c(per[19:24]))),
                 se = c(sd(c(per[17:18]))/sqrt(2), sd(c(per[19:24]))/sqrt(6)),
                 tumor = rep('MYC', 2))

df = rbind.data.frame(df1, df2, df3)
write.csv(df, "Percentage_of_dual_sgRNAs_with_log2RPM>=1_in_cells_tumors.csv")

p <- ggplot(df, aes(x=tumor, y=mean, fill = sample)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,position=position_dodge(0.9))+
  #facet_grid(~ tumor)+
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        axis.line = element_line(colour = "black"),
        #legend.position = 'none',
        #legend.position = c(0.75,0.9),
        #legend.title = element_text(size = 20, face= "bold"),
        legend.text = element_text(size=18,color = 'black'),
        axis.text.y  =element_text(size=20,color = 'black'),
        axis.text.x =element_text(size=20,color = 'black', angle = 45, hjust = 1),
        axis.title=element_text(size=25,color = 'black'),#face="bold")
        plot.title = element_text(size=32,hjust=0.5),
        strip.text  = element_text(size = 16, face="bold")
        ##strip.text is for facet
  )+
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,25))+
  #labs(title="",x ="", y = "Percentage(%) of \n log2 RPM > 0 ")+
  labs(title="",x ="", y = expression(paste('% of ',log[2]~RPM, " >= 1")))+
  scale_fill_manual(values=c("#5AC8FA","magenta"), #"#007AFF", #8E8E93
                    name="",
                    #breaks= ,
                    labels=c('Cell', 'Tumor'))
p
ggsave('Percentage_of_dual_sgRNAs_with_log2RPM>=1_in_cells_tumors.pdf', width = 5.8, height = 5.8)

#----------------------IV.Percentage_of_SKOs_and_DKOS in cells and tumors-----------------------------
#"mSD" is a function to calculate log2RPM per SKOs, DKOs and Control dual-sgRNAs
#Input: "read" is log2rpm, "gB" is genetic background
read = rpm0.CP
gB="PTEN-/-"
mSD = function(read, gB){
  double.CTRL = grepl("^CTRL", row.names(read)) & grepl("-CTRL", row.names(read))
  sum(double.CTRL)
  rd0 = read[!double.CTRL,]
  single <- apply(rd0[grepl("CTRL",rownames(rd0)),], 2, function(x) {sum(x)/sum(grepl("CTRL",row.names(rd0)))})
  double <- apply(rd0[!grepl("CTRL",rownames(rd0)),], 2, function(x) {sum(x)/sum(!grepl("CTRL",row.names(rd0)))})
  ctrl <- apply(read[double.CTRL,], 2, function(x) {sum(x)/sum(double.CTRL)})
  
  #data matrix containing log2RPM/dual-sgRNAs for each sample
  #res = rbind(single, double, ctrl)
  res = rbind(single, double)
  totalM = apply(res, 2, function(x) x/sum(x)*100)
#----------------------------------------------------------------------------------
  #data matrix containing mean log2RPM/dual-sgRNAs for cells and tumors
  df1 = data.frame(sample = 'Cells',
                   type = c("SKO", "DKO"), 
                   mean = apply(totalM[1:2,1:2], 1, function(x) mean(x)),
                   se = apply(totalM[1:2,1:2], 1, function(x) sd(x)/sqrt(length(x))),
                   gB = rep(gB, 2))
  df2 = data.frame( sample = 'Tumors',
                   type= c("SKO", "DKO"), 
                   mean = apply(totalM[1:2,3:ncol(totalM)], 1, function(x) mean(x)),
                   se = apply(totalM[1:2,3:ncol(totalM)], 1, function(x) sd(x)/sqrt(length(x))),
                   gB = rep(gB, 2))
  outM = rbind.data.frame(df1, df2)
  return(outM)
}

mSD.Pten = mSD(read = rpm0.CP[,1:11], gB = "PTEN-/-")
mSD.Pik = mSD(read = rpm0.CP[,12:16], gB = "PIK3CA")
mSD.Myc = mSD(read = rpm0.CP[,17:24], gB = "MYC")

data = rbind.data.frame(mSD.Pten, mSD.Pik, mSD.Myc)
write.csv(data, "Relative_percentage_of_SKOs_and_DKOs_in_cells_tumors.csv")

p <- ggplot(data, aes(x=sample, y=mean, fill=type)) +
  geom_bar(stat="identity", position = 'stack')+
  #geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2)+
  facet_grid(~ gB)+
  theme(panel.border = element_blank(), 
         panel.grid = element_blank(),
         panel.background = element_blank(),
         #panel.grid= element_line(linetype = "dashed"),
         axis.line = element_line(colour = "black"),
         #legend.position = c(0.75,0.9),
         legend.title = element_text(size = 20, face= "bold"),
         legend.text = element_text(size=18,color = 'black'),
         axis.text.y  =element_text(size=20,color = 'black'),
         axis.text.x =element_text(size=20,color = 'black',angle = 45,hjust = 1),
         axis.title=element_text(size=25,color = 'black'),#face="bold")
         plot.title = element_text(size=32,hjust=0.5),
         strip.text  = element_text(size = 18,color = 'black')
         ##strip.text is for facet
)+
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,25))+
  labs(title="",x ="", y = paste("Relative percentage(%)"))+
  scale_fill_manual(values=c("magenta","#5AC8FA"),
                    name="",
                    #breaks=c("Cells", "Tumors"),
                    labels=c("DKOs", "SKOs"))
p
ggsave('Relative_percentage_of_SKOs_and_DKOs_in_cells_tumors.pdf', width = 5.8, height = 5.8)





