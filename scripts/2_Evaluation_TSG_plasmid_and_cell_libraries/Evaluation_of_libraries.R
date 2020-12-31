#Xiaoyu Zhao  
library(ggplot2)
library(RColorBrewer)
library(Cairo)
library("ggrepel")
library(readxl)
setwd("~/Downloads/data/")
#------------------------------------I. Dual-sgRNAs in plasmid library-----------------------------
read0 <- data.frame(read_excel("TableS2_Evaluation_of_TSG_plasmid_cell_libraries.xlsx", 1))
counts = read0[,3:ncol(read0)]
rownames(counts) = read0$gene_names
colnames(counts) = c("P1", "P2", "Vector1", "Vector2","PTEN1", "PTEN2", "PI3K1", "PI3K2", "Myc1", "Myc2")
#How many dual-sgRNAs were detected in plasmid library?
sum(counts$P1+counts$P2>0) #34878
#Matrix only containing dual-sgRNAs detected in plasmid library
countsM = counts[(counts$P1+counts$P2)>0, 1:2]
countsM$P12 = countsM$P1 + countsM$P2
countsM$logrpm = log2(countsM$P12/sum(countsM$P12)*1000000)

#Density plot of dual-sgRNAs in plasmid library
theme_set(theme_bw())
p = ggplot(countsM, aes(x=logrpm)) +
  geom_density(alpha=.4,  fill = "blue")+
  scale_x_continuous(limits = c(0,8))+
  theme(#panel.border = element_blank(), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", size = 1.5),
    axis.line = element_line(colour = "black"),
    legend.position = c(0.25,0.8),
    legend.title = element_text(size = 20),
    legend.text = element_text(size=18, color = "black" ),
    axis.text.y  =element_text(size=22, color = "black"),
    axis.text.x =element_text(size=22, color = "black"),
    axis.title=element_text(size=24, color = "black"),#face="bold")
    plot.title = element_text(size=32,hjust=0.5) )+
  labs(title="", x =expression(log[2]~RPM~per~dual_sgRNA), y = "Density")
p
ggsave("Density.plot.plasmid.libarary.pdf", width = 5.2, height =5.8, device = cairo_pdf)


#-------------------II. CDF on log2RPM of plasmid and cell libraries-------------------------
#Calculate rpm(reads per million)
rpm0 = data.frame(apply(counts, 2, function(x) {(x/sum(x)*1000000)}))
#log2(mean_rpm), "+1" is to avoid log2(0)
logrpm0 = log2(rpm0+1)
#Mean rmp of two replicates
logrpm0$meanP = (logrpm0$P1+logrpm0$P2)/2
logrpm0$meanVector = (logrpm0$Vector1+logrpm0$Vector2)/2
logrpm0$meanPTEN  = (logrpm0$PTEN1+logrpm0$PTEN2)/2 
logrpm0$meanPI3K = (logrpm0$PI3K1+logrpm0$PI3K2)/2
logrpm0$meanMyc = (logrpm0$Myc1+logrpm0$Myc2)/2

#No. of dual-sgRNAs detected in each library
sum(counts$P1+counts$P2>0) #34878
sum((counts$Vector1 + counts$Vector2)>0) #34183
sum((counts$PTEN1 + counts$PTEN2)>0)  #34621
sum((counts$PI3K1 +counts$PI3K2)>0)  #34154
sum((counts$Myc1 +counts$Myc2) >0)    #34185

pdf("CDF_of_meanRPM_plasmid_and_cells.pdf", width = 8, height=5.5)
par(mar = c(5, 5, 2, 1))
plot(ecdf(logrpm0$meanP),col="red",xlim=c(0,10), verticals=T,do.points=F,lwd=4,
     main = NULL,
     ylab= "Cumulative Distribution",
     xlab = expression(log[2]~RPM~per~dual_sgRNA),
     cex.lab=2, 
     cex.axis=1.8)
lines(ecdf(logrpm0$meanVector), col = "magenta",verticals=T,do.points=F,lwd=4)
lines(ecdf(logrpm0$meanPTEN), col = "blue",verticals=T,do.points=F,lwd=4)
lines(ecdf(logrpm0$meanPI3K), col = "turquoise1",verticals=T,do.points=F,lwd=4)
lines(ecdf(logrpm0$meanMyc), col = "orange",verticals=T,do.points=F,lwd=4)
legend(4.5,0.4,c("Plasmid library(34878)",
                  "Vector cell library (34183)", 
                  "PTEN-/- cell library(34621)",
                  "PIK3CA cell library (34154)",
                  "HAmyc cell library (34185)"), 
fill = c("red", "blue","magenta","turquoise1","orange"),  ## "royalblue1",
bty = "n", 
pt.cex = 1.5, 
cex = 1.3, 
text.col = "black"
)
dev.off()

#---------------------------III. Scatter plots of cell libraries vs. Plasmid library---------------------------------
# inputs: counts of 
setwd("~/Downloads/data/")
res <- data.frame(read_excel("TableS2_Evaluation_of_TSG_plasmid_cell_libraries.xlsx", 2))
gene52 = as.character(res[,1])

gScatter_cell_vs_plasmid <- function(d1 = rpm0$meanVector, d2 = rpm0$meanP , title = "Vector", fold = "Vector" ){
  setwd("~/Downloads/data/")
  path = paste0(getwd(), "/", "Enrichment_sgTp53_", fold, "/")
  dir.create(path) 
  setwd(path)
  temp1 = NULL
  temp2 = NULL
  coef_1 = NULL
  coef_2 = NULL
  Pvalue = rep(0,52) 
  for(i in 1:length(gene52)){
    dat = cbind.data.frame(d1, d2)
    rownames(dat) = rownames(rpm0)
    freq <- data.frame(dat, group = rep(0, nrow(dat)))
    freq[grepl(gene52[i],rownames(freq)),3] <- paste0("W/ sg", gene52[i])
    freq[!grepl(gene52[i],rownames(freq)),3] <- paste0("W/O sg", gene52[i])
    dat00 = freq
    
    data1 = subset(dat00, group == paste0("W/ sg", gene52[i]))
    data2 = subset(dat00, group == paste0("W/O sg", gene52[i]))
    
    p = anova(lm(d1~d2, data = dat00), lm(d1~d2*group, data = dat00))
    Pvalue[i] = p$`Pr(>F)`[2]
    
    #Save the intercept and slope of two regressions
    temp1 = t(coef(lm(data1$d1 ~ data1$d2)))
    coef_1 = rbind(coef_1, temp1)
    temp2 = t(coef(lm(data2$d1 ~ data2$d2)))
    coef_2 = rbind(coef_2, temp2)
    
    q <- ggplot(dat00, aes(x=d2, y=d1, group = factor(group)))+
      geom_point(aes(shape = group, color= group),size = 2.2, stroke = 2, alpha = 0.5) + 
      geom_smooth(data = subset(dat00, group == paste0("W/ sg", gene52[i])),method=lm,linetype = "dashed", color= "red", size = 2)+
      geom_smooth(data = subset(dat00, group == paste0("W/O sg", gene52[i])), method= lm,linetype = "dashed", color= "blue", size = 2)+
      scale_color_manual(values=c("red", "lightblue"))+
      scale_shape_manual(values=c(0, 21))+
      #scale_size_manual(values=c(2,3))+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.text.y=element_text(size=36,colour = "black"),
            axis.text.x =element_text(size=36, colour = "black"),
            axis.title=element_text(size=40, colour = "black"),
            plot.title = element_text(size=40,hjust=0.5, colour = "black"),
            legend.title = element_blank(),
            legend.text = element_text(color="black", size=28),
            legend.position = c(0.01,0.999),
            legend.justification = c("left", "top"))+
      scale_x_continuous(limits =c(0,160), breaks = round(seq(0,150, 50),0),name="RPM in plasmid library") +
      scale_y_continuous(name="RPM in cell library")+
      ggtitle(title)
    q
    ggsave(paste0("gScatter_sg",gene52[i],"_", fold, ".pdf"), width = 8, height = 6.8)                                    
  }
  padj = p.adjust(Pvalue, method = "BH")
  outM = cbind(Pvalue, padj, coef_1, coef_2, coef_1[,2]-coef_2[,2], log2(coef_1[,2]/coef_2[,2]))
  colnames(outM) = c("p.value", "p.adj", "Intercept_W/", "slope_W/","Intercept_W/O", "slopeW/O","Difference_slopes", "lfc")
  rownames(outM) = gene52
  return(outM)
  }

#1. Scatter plots of Vector_cell library vs. Plasmid library
outM1 = gScatter_cell_vs_plasmid(d1 = rpm0$meanVector,d2 = rpm0$meanP, title = "Vector", fold = "Vector" )
write.csv(outM1, "Regression_differences_Vector.csv")
#2. Scatter plots of PTEN-/-_cell library vs. Plasmid library
outM2 = gScatter_cell_vs_plasmid(d1 = rpm0$meanPTEN, d2 = rpm0$meanP, title = "PTEN-/-", fold = "PTEN" )
write.csv(outM2, "Regression_differences_PTEN.csv")
#3. Scatter plots of PIK3CA_cell library vs. Plasmid library
outM3 = gScatter_cell_vs_plasmid(d1 = rpm0$meanPI3K, d2 = rpm0$meanP, title = "PIK3CA", fold = "PIK3CA" )
write.csv(outM3, "Regression_differences_PIK3CA.csv")
#4. Scatter plots of HAMyc_cell library vs. Plasmid library
outM4 = gScatter_cell_vs_plasmid(d1 = rpm0$meanMyc, d2 = rpm0$meanP, title = "HAMyc", fold = "HAMyc" )
write.csv(outM4, "Regression_differences_HAMyc.csv")

##------------------------------IV. Enrichment of sgTP53 containing dual-sgRNAs-------------------------------------
setwd("~/Downloads/data/")
readV = read.csv("Enrichment_sgTp53_Vector/Regression_differences_Vector.csv", row.names =1 )
readM = read.csv("Enrichment_sgTp53_PTEN/Regression_differences_PTEN.csv", row.names =1  )
readPi3k = read.csv("Enrichment_sgTp53_PIK3CA/Regression_differences_PIK3CA.csv", row.names =1  )
readPten = read.csv("Enrichment_sgTp53_HAMyc/Regression_differences_HAMyc.csv", row.names =1  )

lfc = data.frame(readV[gene52,"lfc"], readM[gene52,"lfc"], readPi3k[gene52,"lfc"], readPten[gene52,"lfc"])
lfc$mean = apply(lfc, 1, function(x) mean(x[1], x[2], x[3], x[4]))

rownames(lfc) = gene52
lfc$sname = gene52
lfc$highlight = "NO"
lfc_order = lfc[order(lfc$mean, decreasing = T),]
head(lfc_order,2) ##TP53, ATM, CSMD1
tail(lfc_order,3) ##BRCA2, USP9X
highlight = c(rownames(head(lfc_order,1)), rownames(tail(lfc_order,2)))
lfc[rownames(lfc) %in% highlight, "highlight"] = "Yes"

Names = factor(gene52, levels= rownames(lfc[order(lfc$mean, decreasing = T),]))
p <- ggplot(lfc, aes(x = Names, y = mean))+
  geom_point(aes(size=highlight, color=highlight, group = highlight)) +
  scale_size_manual(values = c(3,4)) +
  scale_color_manual(values= c("blue","red"))+
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        #legend.title = element_text(size = 14),
        #legend.text = element_text(size=14 ),
        axis.text.y  =element_text(size=18, color = "black"),
        axis.text.x = element_blank(),
        #axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title=element_text(size=20, color = "black"),#face="bold")
        plot.title = element_text(size=32,hjust=0.5)
        #strip.text.y  = element_text(size = 26, face="bold")
        ##strip.text is for facet
  )+geom_hline(yintercept  = c(0, 0.5, -0.5), linetype = "dashed", color = "grey30", size = 1)+
  labs(title= "", x ="sgRNAs at gene level", y = expression(log[2]~(slope1/slope2)))
p
ggsave("Regression_differences_52genes.pdf" , width = 6.5, height=4.3)



 