#Xiaoyu Zhao 
set.seed(100)
library(ggplot2)
library(data.table) ## for read.table("xxx.xls",sep = "\t",header = T) 
library(Cairo) #https://stackoverflow.com/questions/27690729/greek-letters-symbols-and-line-breaks-inside-a-ggplot-legend-label
library(ggpubr)
library(tidyr)
library(dplyr)
library(scales)
library(eulerr)
setwd("~/Downloads/5_Combinatorial_CROP_seq/")
tI_1 =read.csv("CROP_Python/sequencing/outs_S4/results/prediction_and_transcriptomic_differences/diff_NF2-PTEN.csv", row.names = 1)
tI_2 =read.csv("CROP_Python/sequencing/outs_S4/results/prediction_and_transcriptomic_differences/diff_NF2-TP53.csv", row.names = 1)
tI_3 =read.csv("CROP_Python/sequencing/outs_S4/results/prediction_and_transcriptomic_differences/diff_PTEN-TP53.csv", row.names = 1)

#-------------------------------I. Categorization of transcriptomic differences-------------------------------------
#Upregulated DEGs  (LFC > lfc_threshold)are grouped into three categories: Synergistic, Additive and Buffering based on transcriptomic differences
#Downregulated DEGs (LFC < -lfc_threshold) are grouped into three categories: Synergistic, Additive and Buffering based on transcriptomic differences
#Function to count number of each category
countI = function(tI = tI_1, lfc_threshold = 0.03, diff_threshold = 0.05){
  print(paste0("Number of Upregulated gene expressions (LFC(DKO) > ", lfc_threshold, ") : ", length(which(tI$double_lfc > lfc_threshold))))
  
  out = list()
  out[["syn_up"]] = tI[which(tI$double_lfc > lfc_threshold & tI$diff_lm >= diff_threshold),] 
  out[["buff_up"]] = tI[which(tI$double_lfc > lfc_threshold & tI$diff_lm <= -diff_threshold),] 
  out[["add_up"]] = tI[which(tI$double_lfc > lfc_threshold & abs(tI$diff_lm) < diff_threshold),] 
  
  print("**********")
  print(paste0("Number of Downregulated gene expressions (LFC(DKO) < ", lfc_threshold, ") : ", length(which(tI$double_lfc < -lfc_threshold))))
  
  out[["syn_down"]] = tI[which(tI$double_lfc < -lfc_threshold & tI$diff_lm <= -diff_threshold),] 
  out[["buff_down"]] = tI[which(tI$double_lfc < -lfc_threshold & tI$diff_lm >= diff_threshold),] 
  out[["add_down"]] = tI[which(tI$double_lfc < -lfc_threshold & abs(tI$diff_lm) < diff_threshold),] 
  
  tI$Category = paste0("Non_DEGs, abs(LFC) < ", lfc_threshold)
  tI$Category[which(tI$double_lfc > lfc_threshold & tI$diff_lm >= diff_threshold)] = "Synergistic_Up"
  tI$Category[which(tI$double_lfc > lfc_threshold & tI$diff_lm <= -diff_threshold)]  = "Buffering_Up"
  tI$Category[which(tI$double_lfc > lfc_threshold & abs(tI$diff_lm) < diff_threshold)]  = "Additive_Up"
  tI$Category[which(tI$double_lfc < -lfc_threshold & tI$diff_lm <= -diff_threshold)]  = "Synergistic_Dn"
  tI$Category[which(tI$double_lfc < -lfc_threshold & tI$diff_lm >= diff_threshold)]  = "Buffering_Dn"
  tI$Category[which(tI$double_lfc < -lfc_threshold & abs(tI$diff_lm) < diff_threshold)]  = "Additive_Dn"
  return(tI)
}

out1 = countI(tI = tI_1, lfc_threshold = 0.03, diff_threshold = 0.05) #NF2-PTEN
out2 = countI(tI = tI_2, lfc_threshold = 0.03, diff_threshold = 0.05) #NF2-TP53
out3 = countI(tI = tI_3, lfc_threshold = 0.03, diff_threshold = 0.05) #PTEN-TP53

#-------------------------------II. LFCs of Additive, Buffering and Synergistic DEGs in the three DKOs-------------------------------------
#-------------Upregulated DEGs----------------
#LFCs of additive-upregulated DEGs
d1 = data.frame(value = c(mean(out1$pred_lm[which(out1$Category == "Additive_Up")]), mean(out2$pred_lm[which(out2$Category == "Additive_Up")]), mean(out3$pred_lm[which(out3$Category == "Additive_Up")]),
                          mean(out1$double_lfc[which(out1$Category == "Additive_Up")]), mean(out2$double_lfc[which(out2$Category == "Additive_Up")]), mean(out3$double_lfc[which(out3$Category == "Additive_Up")])))
d1$Catergory = rep("Additive", 6)
d1$regression = c(rep("Fit", 3), rep("DKO", 3)) 
d1$name =rep(c("NF2-PTEN", "NF2-TP53", "PTEN-TP53"), 2)
#LFCs of synergistic-upregulated DEGs
d2 = data.frame(value = c(mean(out1$pred_lm[which(out1$Category == "Synergistic_Up")]), mean(out2$pred_lm[which(out2$Category == "Synergistic_Up")]), mean(out3$pred_lm[which(out3$Category == "Synergistic_Up")]),
                          mean(out1$double_lfc[which(out1$Category == "Synergistic_Up")]), mean(out2$double_lfc[which(out2$Category == "Synergistic_Up")]), mean(out3$double_lfc[which(out3$Category == "Synergistic_Up")])))
d2$Catergory = rep("Synergistic", 6)
d2$regression = c(rep("Fit", 3), rep("DKO", 3)) 
d2$name =rep(c("NF2-PTEN", "NF2-TP53", "PTEN-TP53"), 2)
#LFCs of Buffering-upregulated DEGs
d3 = data.frame(value = c(mean(out1$pred_lm[which(out1$Category == "Buffering_Up")]), mean(out2$pred_lm[which(out2$Category == "Buffering_Up")]), mean(out3$pred_lm[which(out3$Category == "Buffering_Up")]),
                          mean(out1$double_lfc[which(out1$Category == "Buffering_Up")]), mean(out2$double_lfc[which(out2$Category == "Buffering_Up")]), mean(out3$double_lfc[which(out3$Category == "Buffering_Up")])))
d3$Catergory = rep("Buffering", 6)
d3$regression = c(rep("Fit", 3), rep("DKO", 3)) 
d3$name =rep(c("NF2-PTEN", "NF2-TP53", "PTEN-TP53"), 2)
#Combine three catergories
lfc_up = rbind.data.frame(d1, d2, d3)
#lfc_up$catergory = factor(lfc_up$Catergory, levels = c( "Additive", "Synergistic", "Buffering"))
lfc_up$Regression = factor(lfc_up$regression, levels = c( "Fit", "DKO"))

ggplot(lfc_up, aes(x = Regression, y = value, fill =  Regression))+
  geom_boxplot( width = 0.6, outlier.size = 0.4, outlier.alpha = 0.2) +
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        axis.line = element_line( color = "black"),
        #legend.position = c(2,0.9),
        legend.position = "None",
        #legend.title = element_blank(),
        legend.text = element_text(size=12 , color = "black"),
        axis.text.y  =element_text(size=14, color = "black"),
        axis.text.x =element_text(size=13,  color = "black"),
        axis.title=element_text(size=14, color = "black"),
        plot.title = element_text(size=14,hjust=0.5,vjust = 2, color = "black"),
        strip.text  = element_text(size = 14, color = "black")
  )+
  scale_fill_brewer(direction = 1)+
  labs(title="Upregulated DEGs",x ="", y = "LFCs")+
  facet_grid(.~Catergory)
ggsave("CROP_R/4_Transcriptional_interactions/LFCs_Additive_Buffering_Synergistic_upregulated_DEGs_NF2-PTEN_TP53.pdf", width = 4.2, height = 3.2)

#-------------Downregulated DEGs----------------
#LFCs of additive-downregulated DEGs
d1 = data.frame(value = c(mean(out1$pred_lm[which(out1$Category == "Additive_Dn")]), mean(out2$pred_lm[which(out2$Category == "Additive_Dn")]), mean(out3$pred_lm[which(out3$Category == "Additive_Dn")]),
                          mean(out1$double_lfc[which(out1$Category == "Additive_Dn")]), mean(out2$double_lfc[which(out2$Category == "Additive_Dn")]), mean(out3$double_lfc[which(out3$Category == "Additive_Dn")])))
d1$Catergory = rep("Additive", 6)
d1$regression = c(rep("Fit", 3), rep("DKO", 3)) 
d1$name =rep(c("NF2-PTEN", "NF2-TP53", "PTEN-TP53"), 2)
#LFCs of Synergistic-downregulated DEGs
d2 = data.frame(value = c(mean(out1$pred_lm[which(out1$Category == "Synergistic_Dn")]), mean(out2$pred_lm[which(out2$Category == "Synergistic_Dn")]), mean(out3$pred_lm[which(out3$Category == "Synergistic_Dn")]),
                          mean(out1$double_lfc[which(out1$Category == "Synergistic_Dn")]), mean(out2$double_lfc[which(out2$Category == "Synergistic_Dn")]), mean(out3$double_lfc[which(out3$Category == "Synergistic_Dn")])))
d2$Catergory = rep("Synergistic", 6)
d2$regression = c(rep("Fit", 3), rep("DKO", 3)) 
d2$name =rep(c("NF2-PTEN", "NF2-TP53", "PTEN-TP53"), 2)
#LFCs of Buffering-downregulated DEGs
d3 = data.frame(value = c(mean(out1$pred_lm[which(out1$Category == "Buffering_Dn")]), mean(out2$pred_lm[which(out2$Category == "Buffering_Dn")]), mean(out3$pred_lm[which(out3$Category == "Buffering_Dn")]),
                          mean(out1$double_lfc[which(out1$Category == "Buffering_Dn")]), mean(out2$double_lfc[which(out2$Category == "Buffering_Dn")]), mean(out3$double_lfc[which(out3$Category == "Buffering_Dn")])))
d3$Catergory = rep("Buffering", 6)
d3$regression = c(rep("Fit", 3), rep("DKO", 3)) 
d3$name =rep(c("NF2-PTEN", "NF2-TP53", "PTEN-TP53"), 2)
#Combine three catergories
lfc_dn = rbind.data.frame(d1, d2, d3)
#lfc_dn$catergory = factor(lfc_dn$Catergory, levels = c( "Additive", "Synergistic", "Buffering"))
lfc_dn$Regression = factor(lfc_dn$regression, levels = c( "Fit", "DKO"))

ggplot(lfc_dn, aes(x = Regression, y = -value, fill =  Regression))+
  geom_boxplot( width = 0.6, outlier.size = 0.4, outlier.alpha = 0.2) +
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        axis.line = element_line( color = "black"),
        #legend.position = c(2,0.9),
        legend.position = "None",
        #legend.title = element_blank(),
        legend.text = element_text(size=12 , color = "black"),
        axis.text.y  =element_text(size=14, color = "black"),
        axis.text.x =element_text(size=13,  color = "black"),
        axis.title=element_text(size=14, color = "black"),
        plot.title = element_text(size=14,hjust=0.5,vjust = 2, color = "black"),
        strip.text  = element_text(size = 14, color = "black")
  )+
  scale_fill_brewer(direction = 1)+
  labs(title="Downregulated DEGs",x ="", y = "- LFCs")+
  facet_grid(.~Catergory)
ggsave("CROP_R/4_Transcriptional_interactions/LFCs_Additive_Buffering_Synergistic_downregulated_DEGs_NF2-PTEN_TP53.pdf", width = 4.2, height = 3.2)


#-----------------------------III. Compare and Analyze the six catergories in NF2-PTEN, NF2-TP53 and PTEN-TP53------------------------------
sname = rownames(out1)
Outs = data.frame(out1$index, out1$first_lfc, out1$second_lfc, out2[sname, "second_lfc"], out1$double_lfc, out2[sname, "double_lfc"], out3[sname, "double_lfc"],
                  out1$pred_lm, out2[sname, "pred_lm"], out3[sname, "pred_lm"], out1$diff_lm, out2[sname, "diff_lm"], out3[sname, "diff_lm"], 
                  out1$Category, out2[sname, "Category"], out3[sname, "Category"])
rownames(Outs) = sname
colnames(Outs) = c("Ensembl_ID", "lfc_NF2", "lfc_PTEN", "lfc_TP53", "double_lfc_NF2-PTEN", "double_lfc_NF2-TP53", "double_lfc_PTEN-TP53","pred_NF2-PTEN", "pred_NF2-TP53",
                   "pred_PTEN-TP53", "diff_NF2-PTEN", "diff_NF2-TP53", "diff_PTEN-TP53", "Category_NF2-PTEN", "Category_NF2-TP53", "Category_PTEN-TP53")
write.csv(Outs, "CROP_R/4_Transcriptional_interactions/TI_categorization_5069_genes_based_on_TDs_3DKOs.csv")

#Function "findCommontI" generate the unioned gene list and venn diagram of each catergory
findCommontI = function(Category = "Synergistic_Up", Outs){
  f1 = Outs[which(Outs$'Category_NF2-PTEN' == Category),]
  f2 = Outs[which(Outs$'Category_NF2-TP53' == Category),]
  f3 = Outs[which(Outs$'Category_PTEN-TP53' == Category),]
  All_names = unique(c(rownames(f1), rownames(f2), rownames(f3)))
  dataM = data.frame(matrix(0, nrow = length(All_names), ncol = 3))
  rownames(dataM) = All_names
  colnames(dataM) = c('NF2-PTEN', 'NF2-TP53', 'PTEN-TP53')
  dataM[match(rownames(f1), All_names), 1] = 1
  dataM[match(rownames(f2), All_names), 2] = 1
  dataM[match(rownames(f3), All_names), 3] = 1
  
  dataN = dataM
  dataN$Group = NA
  dataN$Group[dataN$'NF2-PTEN' == 1 & dataN$'NF2-TP53' == 1 & dataN$'PTEN-TP53' == 1] = paste0(Category,": Common")
  dataN$Group[dataN$'NF2-PTEN' == 1 & dataN$'NF2-TP53' == 1 & dataN$'PTEN-TP53' == 0] = paste0(Category,": NF2_specific")
  dataN$Group[dataN$'NF2-PTEN' == 1 & dataN$'NF2-TP53' == 0 & dataN$'PTEN-TP53' == 1] = paste0(Category,": PTEN_specific")
  dataN$Group[dataN$'NF2-PTEN' == 0 & dataN$'NF2-TP53' == 1 & dataN$'PTEN-TP53' == 1] = paste0(Category,": TP53_specific")
  dataN$Group[dataN$'NF2-PTEN' == 1 & dataN$'NF2-TP53' == 0 & dataN$'PTEN-TP53' == 0] = paste0(Category,": NF2-PTEN_specific")
  dataN$Group[dataN$'NF2-PTEN' == 0 & dataN$'NF2-TP53' == 1 & dataN$'PTEN-TP53' == 0] = paste0(Category,": NF2-TP53_specific")
  dataN$Group[dataN$'NF2-PTEN' == 0 & dataN$'NF2-TP53' == 0 & dataN$'PTEN-TP53' == 1] = paste0(Category,": PTEN-TP53_specific")
  
  write.csv(dataN, file = paste0("CROP_R/4_Transcriptional_interactions/Unioned_gene_list_", Category, ".csv"))
  print(paste0("Total number of unique ", Category, ":", length(All_names)))
  #print(paste0("Number of common ", Category, ":", dim(dataM[which(rowSums(dataM) == 3),])[1]))
  
   #if(overLap(a[1:2],data = dataM) + overLap(a[2:3],data = dataM) + overLap(a[c(1,3)],data = dataM) >= 3 ){
    pdf(paste0("CROP_R/4_Transcriptional_interactions/Venn_diagram_common_",Category,".pdf"), height = 5, width = 5)
    options(warn = -1)
    g = euler(dataM, shape = "ellipse")
    
    p <- plot(g, fills = list(fill = c("mediumorchid1","orange","#5AC8FA"), alpha = 0.6, size = 2), quantities = TRUE, 
         legend = list(labels = c("NF2-PTEN", "NF2-TP53", "PTEN-TP53"), cex = 1.4, alpha = 0.8))
    print(p)
    
    dev.off()
 # } else{
    #print("No venn diagram")
 # }
  return(c(nrow(f1), nrow(f2), nrow(f3), length(All_names), g$original.values))
}

#Function "Count_Common_tIs" finds common genes of all the six categories
Count_Common_tIs = function(tI_1 = tI_1, tI_2 = tI_2, tI_3 = tI_3, lfc_threshold = 0.03, diff_threshold = 0.05){
  print("In DKO1:")
  out1 = countI(tI = tI_1, lfc_threshold = lfc_threshold, diff_threshold = diff_threshold)
  print("In DKO2:")
  out2 = countI(tI = tI_2, lfc_threshold = lfc_threshold, diff_threshold = diff_threshold)
  print("In DKO3:")
  out3 = countI(tI = tI_3, lfc_threshold = lfc_threshold, diff_threshold = diff_threshold)
  
  #Generate three venn_diagrams and saved into pdf
  res1 = findCommontI(Category = "Synergistic_Up", Outs = Outs)
  res2 = findCommontI(Category = "Buffering_Up", Outs = Outs)
  res3 = findCommontI(Category = "Additive_Up", Outs = Outs)
  res4 = findCommontI(Category = "Synergistic_Dn", Outs = Outs)
  res5 = findCommontI(Category = "Buffering_Dn", Outs = Outs)
  res6 = findCommontI(Category = "Additive_Dn", Outs = Outs)
  
  outM = rbind(res1, res2, res3, res4, res5, res6)
  colnames(outM) = c("NF2-PTEN", "NF2-TP53", "PTEN-TP53", "num_union", "num_NF2-PTEN_specific", "num_NF2-TP53_specific", "num_PTEN-TP53_specific",
                     "num_NF2_specific","num_PTEN_specific", "num_TP53_specific", "num_commonn" )
  rownames(outM) = c("Syn_up", "Buff_up", "Add_up", "Syn_down", "Buff_down", "Add_down")
  return(outM)
}

Res = Count_Common_tIs(tI_1 = tI_1, tI_2 = tI_2, tI_3 = tI_3, lfc_threshold = 0.03, diff_threshold = 0.05)
Res
write.csv(Res, "CROP_R/4_Transcriptional_interactions/Num_summary_TI_Categorization_based_on_TDs_3DKOs.csv")

#-------------------------------V. Composition of each category in every DKO-------------------------------------
#Figure 1: Pie chart: composition of each category in every DKO
##Input: Change the row and column to get information of different DKO
#For example:
#Res[1:3,"NF2-PTEN"] Synergistic_upregulation in NF2-PTEN
#Res[4:6,"NF2-PTEN"] Synergistic_upregulation in NF2-PTEN
#NF2-TP53
num = Res[1:3,"NF2-PTEN"] 
num = Res[4:6,"NF2-PTEN"] 
#NF2-PTEN
num = Res[1:3,"NF2-TP53"] 
num = Res[4:6,"NF2-TP53"]
#PTEN-TP53
num = Res[1:3,"PTEN-TP53"] 
num = Res[4:6,"NF2-TP53"] 
sum(num)

dat = data.frame(Category = c( "Synergistic", "Buffering", "Additive"),
                 num,
                 #num = c(7,13,17,8,1),
                 prop = round(num/sum(num),3)
                 #lab.ypos = cumsum(dat$per) - 0.8*dat$per
)

Category0= factor(dat$Category, level = c("Synergistic", "Additive","Buffering"))
ggplot(dat, aes(x = 2, y = prop, fill = Category0))+
  geom_bar(width = 1,stat = "identity", color = "white")+
  coord_polar("y", start = 0)+
  geom_text(aes(label = sprintf("%1.1f%%", 100*prop)), color = "Black", size = 9, position = position_stack(vjust = 0.5))+
  xlim(1.5, 2.5)+
  theme_void()+
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        #axis.line = element_line(colour = "black"),
        legend.position = c(1.2,0.6),
        legend.title = element_text(size = 20),
        legend.text = element_text(size=20 ),
        #axis.text.y  =element_text(size=24,face= "bold"),
        #axis.text.x =element_text(size=24, face = "bold"),#, margin = margin(t =20, r = 10, b = 10, l = 10)),
        #axis.title=element_text(size=30),#face="bold")
        plot.title = element_text(size=22,hjust=0.5)
        #strip.text.y  = element_text(size = 26, face="bold")
        ##strip.text is for facet
  )+
  scale_fill_manual(values = c("#CD534CFF", "#EFC000FF","#0073C2FF"),
                    # values = c("darksalmon", "#EFC000FF","#0073C2FF"),
                    name="Category",
                    #breaks=c("Cells", "Tumors"),
                    labels=c( "Synergistic", "Additive", "Buffering"))
#labs(title="PTEN-/-",x ="", y = "")
ggsave("CROP_R/4_Transcriptional_interactions/Category_TI_downregulated_genes_PTEN-TP53.pdf", width=8, height = 4)
#ggsave("CROP_R/4_Transcriptional_interactions/Category_TI_downregulated_genes_PTEN_TP53.pdf", width=8, height = 4)

#---------------------IV. Group TI_categories in 3SKOs: Commonly shared, gene-specific and DKO-specific--------------------------------------------
#Each catergory of transcriptional differences is divided into "G1: Commonly shared", "G2: NF2 specific", "G3: PTEN specific", 
#"G4: TP53 specific", "G5: NF2-PTEN specific", "G6: NF2-TP53 specific", "G7: PTEN-TP53 specific" based on the venn diagram.
tIs = read.csv("CROP_R/4_Transcriptional_interactions/TI_categorization_5069_genes_based_on_TDs_3DKOs.csv", row.names = 1)
tIs$Group_Additive_Up = NA
tIs$Group_Synergistic_Up = NA
tIs$Group_Buffering_Up = NA
tIs$Group_Additive_Dn = NA
tIs$Group_Synergistic_Dn = NA
tIs$Group_Buffering_Dn= NA

add_up = read.csv("CROP_R/4_Transcriptional_interactions/Unioned_gene_list_Additive_Up.csv", row.names = 1)    #1858
syn_up = read.csv("CROP_R/4_Transcriptional_interactions/Unioned_gene_list_Synergistic_Up.csv", row.names = 1) #1476
buff_up = read.csv("CROP_R/4_Transcriptional_interactions/Unioned_gene_list_Buffering_Up.csv", row.names = 1)  #173
add_down = read.csv("CROP_R/4_Transcriptional_interactions/Unioned_gene_list_Additive_Dn.csv", row.names = 1)  #722
syn_down = read.csv("CROP_R/4_Transcriptional_interactions/Unioned_gene_list_Synergistic_Dn.csv", row.names = 1) #329
buff_down = read.csv("CROP_R/4_Transcriptional_interactions/Unioned_gene_list_Buffering_Dn.csv", row.names = 1) #94

tIs[match(rownames(add_up), rownames(tIs)), "Group_Additive_Up"] = as.character(add_up$Group)
tIs[match(rownames(syn_up), rownames(tIs)), "Group_Synergistic_Up"] = as.character(syn_up$Group)
tIs[match(rownames(buff_up), rownames(tIs)), "Group_Buffering_Up"] = as.character(buff_up$Group)
tIs[match(rownames(add_down), rownames(tIs)), "Group_Additive_Dn"] = as.character(add_down$Group)
tIs[match(rownames(syn_down), rownames(tIs)), "Group_Synergistic_Dn"] = as.character(syn_down$Group)
tIs[match(rownames(buff_down), rownames(tIs)), "Group_Buffering_Dn"] = as.character(buff_down$Group)

write.csv(tIs, "CROP_R/4_Transcriptional_interactions/Groups_of_TI_categorization_5069_genes_based_on_TDs_3DKOs.csv")

