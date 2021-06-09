#Xiaoyu Zhao 
library(ggplot2)
library(Cairo)
library("ggrepel")
library(readxl)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
setwd("~/Downloads/5_Combinatorial_CROP_seq/")

#@@@@@@@@@@@@@---------------I. Synergistic expression of CDK4, DNMT1 and SRPK1 in DKOs--------------@@@@@@@@@@@@@@@
read0 <- data.frame(read_excel("TableS7_Transcriptional_interactions_underlying_GIs.xlsx", 4))
data = read0[,2:ncol(read0)]
rownames(data) = read0[,1]
d = data
#----------------Expression CDK4 in TP53-specific upregulated synergistic genes---------------------------------
d_CDK4 = d["CDK4",c(2:4,6:7)]
colnames(d_CDK4) = c("NF2","PTEN", "TP53",  "NF2-TP53", "PTEN-TP53")
lfc_CDK4 = melt(d_CDK4)
lfc_CDK4$group = factor(lfc_CDK4$variable, levels = c("NF2", "PTEN", "TP53","NF2-TP53", "PTEN-TP53"))
ggplot(lfc_CDK4, aes(x = group, y = value, fill = group))+
  geom_bar(stat="identity", width = 0.7, position=position_dodge(),colour="black")+
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        legend.position = "None",
        axis.line = element_line(colour = "black"),
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size=15 , color = "black"),
        axis.text.y  =element_text(size=15, color = "black"),
        axis.text.x =element_text(size=15, angle = 45, vjust = 1, hjust= 1, color = "black"),
        axis.title.x = element_text(size=16),
        axis.title.y  =element_text(size=17),
        plot.title = element_text(size=16, hjust = 0.5)
  )+
  labs(title="",x ="", y = "LFCs of CDK4")+
  scale_y_continuous(limits = c(0, 0.48), breaks = seq(0, 0.5,0.1))
ggsave("CROP_R/4_Transcriptional_interactions/Expression_CDK4_TP53_specific_upregulated_synergy.pdf", width = 2.4, height =4)

#----------------Expression DNMT1 in TP53-specific upregulated synergistic genes---------------------------------
d_DNMT1 = d["DNMT1",c(2:4,6:7)]
colnames(d_DNMT1) = c("NF2","PTEN", "TP53",  "NF2-TP53", "PTEN-TP53")
lfc_DNMT1 = melt(d_DNMT1)
lfc_DNMT1$group = factor(lfc_DNMT1$variable, levels = c("NF2", "PTEN", "TP53","NF2-TP53", "PTEN-TP53"))
ggplot(lfc_DNMT1, aes(x = group, y = value, fill = group))+
  geom_bar(stat="identity", width = 0.7, position=position_dodge(),colour="black")+
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        legend.position = "None",
        axis.line = element_line(colour = "black"),
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size=15 , color = "black"),
        axis.text.y  =element_text(size=15, color = "black"),
        axis.text.x =element_text(size=15, angle = 45, vjust = 1, hjust= 1, color = "black"),
        axis.title.x = element_text(size=16),
        axis.title.y  =element_text(size=17),
        plot.title = element_text(size=16, hjust = 0.5)
  )+
  labs(title="",x ="", y = "LFCs of DNMT1")+
  scale_y_continuous(limits = c(0, 0.85), breaks = seq(0, 0.8,0.2))
ggsave("CROP_R/4_Transcriptional_interactions/Expression_DNMT1_TP53_specific_upregulated_synergy.pdf", width = 2.4, height =4)

#----------------Expression SRPK1 in TP53-specific upregulated synergistic genes---------------------------------
d_SRPK1 = d["SRPK1",c(2,4,6)]
colnames(d_SRPK1) = c("NF2", "TP53",  "NF2-TP53")
lfc_SRPK1 = melt(d_SRPK1)
lfc_SRPK1$group = factor(lfc_SRPK1$variable, levels = c("NF2", "TP53",  "NF2-TP53"))
ggplot(lfc_SRPK1, aes(x = group, y = value, fill = group))+
  geom_bar(stat="identity", width = 0.7, position=position_dodge(),colour="black")+
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        legend.position = "None",
        axis.line = element_line(colour = "black"),
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size=15 , color = "black"),
        axis.text.y  =element_text(size=15, color = "black"),
        axis.text.x =element_text(size=15, angle = 45, vjust = 1, hjust= 1, color = "black"),
        axis.title.x = element_text(size=16),
        axis.title.y  =element_text(size=17),
        plot.title = element_text(size=16, hjust = 0.5)
  )+
  labs(title="",x ="", y = "LFCs of SRPK1")+
  scale_y_continuous(limits = c(-0.005, 0.08), breaks = seq(0, 0.08, 0.02))
ggsave("CROP_R/4_Transcriptional_interactions/Expression_SRPK1_NF2_TP53_specific_upregulated_synergy.pdf", width = 2, height =4)


#@@@@@@@@@@@@@---------------II. Synergistic expression of CDK4, DNMT1 and SRPK1 in patients-------------@@@@@@@@@@@@@@@
read0 <- data.frame(read_excel("TableS7_Transcriptional_interactions_underlying_GIs.xlsx", 8))
data = read0[,2:ncol(read0)]
rownames(data) = read0[,1]
d = data
#----------------------------------1. Synergy of CDK4 in patients with NF2-TP53 mutations-----------------------------------------------
target = "CDK4"

unique(d$YAP1_CNA)
unique(d$TP53_WT_MUT)
#DKO: CDK4 expression in TP53-Mut YAP1-Gain/Amp or WWTR1-Gain/Amp
exp1 = d[d$TP53_WT_MUT %in% "Mutated" & d$YAP1_CNA %in% c("YAP1: Gain", "YAP1: Amplification"), target]
exp2 = d[d$TP53_WT_MUT %in% "Mutated" & d$WWTR1_CNA %in% c("WWTR1: Gain", "WWTR1: Amplification"), target]
#SKO1: CDK4 expression in TP53-Mut YAP1-Diploid/Amp and WWTR1-Diploid
exp3 = d[(d$TP53_WT_MUT %in% "Mutated"&
            d$YAP1_CNA %in% c("YAP1: Diploid") & 
            d$WWTR1_CNA %in% c("WWTR1: Diploid")), target]
#SKO2: CDK4 expression in TP53-WT YAP1-Gain/Amp or WWTR1-Gain/Amp
exp4 = d[d$TP53_WT_MUT %in% "Wild type" & d$YAP1_CNA %in% c("YAP1: Gain", "YAP1: Amplification"), target]
exp5 = d[d$TP53_WT_MUT %in% "Wild type" & d$WWTR1_CNA %in% c("WWTR1: Gain", "WWTR1: Amplification"), target]

#Control: CDK4 expression in TP53-WT YAP1-Diploid/Amp and WWTR1-Diploid
exp6 = d[(d$TP53_WT_MUT %in% "Wild type"&
            d$YAP1_CNA %in% c("YAP1: Diploid") & 
            d$WWTR1_CNA %in% c("WWTR1: Diploid")) , target]

df = rbind.data.frame(data.frame(expr = c(exp1, exp2), catergory = rep("TP53-Mut YAP1/TAZ-Gain/Amp", length(exp1) + length(exp2))),
                      data.frame( expr = exp3, catergory = rep("TP53-Mut YAP1&TAZ-Diploid", length(exp3))),
                      data.frame( expr = c(exp4, exp5), catergory = rep("TP53-Wt YAP1/TAZ-Gain/Amp", length(exp4) + length(exp5))),
                      data.frame( expr = exp6, catergory = rep("TP53-Wt YAP1&TAZ-Diploid", length(exp6))))

ggplot(df,  aes(x = catergory, y = expr, fill = catergory)) +
  stat_boxplot(geom ='errorbar', outlier.shape = NA, width = 0.3, lwd = 0.2) +
  geom_boxplot(width = 0.6, outlier.shape = NA, lwd = 0.35)+ 
  theme(#panel.border = element_blank(), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line = element_line(colour = "black"),
    # legend.position = c(0.8,0.7),
    legend.position = "None",
    legend.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size=11, color = "black" ),
    axis.text.y = element_text(size=13, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title=element_text(size=14, color = "black"),
    plot.title = element_text(size=14,hjust=0.5, color = "black"),
    #legend.key = element_rect(fill = "lightblue", color = NA),
    # Change legend key size and key width
    legend.key.size = unit(0.25, "cm")
    #strip.text.y  = element_text(size = 26, face="bold")
    ##strip.text is for facet
  )+
  scale_x_discrete(name ="", labels ="") +
  labs(title="",x ="", y = paste0("Expression of ", target))+
  guides(fill=guide_legend(title="Types of Genetic Alterations"))+
  #scale_fill_discrete(#palette =  "Pastel2",
  #name="Types of Genetic Alterations",
  #breaks=c("Cells", "Tumors"),
  #labels=c("TP53- YAP/TAZ+", "20-40%", "40-60%", "60-80%", "80-100%"))+
  geom_hline(yintercept=seq(9.5, 10.5, 1), linetype="dashed", color = "grey20")+
  scale_y_continuous(limits = c(8.7,11), breaks = seq(8.5, 11, 0.5))
ggsave(paste0("CROP_R/4_Transcriptional_interactions/Synergy_", target, "_upregulation_in_NF2_TP53_patients.pdf"), width = 2.1, height = 2.5)

t.test(c(exp1,exp2), exp3, alternative = "greater") #4.408e-06
t.test(c(exp1,exp2), c(exp4, exp5), alternative = "greater") #3.918e-05

#----------------------------------2. Synergy of CDK4 in patients with PTEN-TP53 mutations-----------------------------------------------
unique(d$PTEN_WT_MUT)
unique(d$PTEN_CNA)
#DKO: CDK4 expression in TP53-Mut PTEN-Mut&Homdel
exp1 = d[d$TP53_WT_MUT %in% "Mutated" &
           (d$PTEN_WT_MUT %in% "Mutated" | d$PTEN_CNA %in% "PTEN: Deep Deletion"), target]
#SKO1: CDK4 expression in TP53-Mut PTEN-Wt/Diploid
exp2 = d[d$TP53_WT_MUT %in% "Mutated" &
           d$PTEN_WT_MUT %in% "Wild type" & d$PTEN_CNA %in% "PTEN: Diploid" , target]

#SKO2: CDK4 expression in TP53-Wt PTEN-Mut/Homdel
exp3 = d[d$TP53_WT_MUT %in% "Wild type" & 
           (d$PTEN_WT_MUT %in% "Mutated" | d$PTEN_CNA %in% "PTEN: Deep Deletion"), target]

#Control: CDK4 expression in TP53-Wt PTEN-Wt/Diploid
exp4 = d[d$TP53_WT_MUT %in% "Wild type" &
           d$PTEN_WT_MUT %in% "Wild type" & d$PTEN_CNA %in% "PTEN: Diploid" , target]


df = rbind.data.frame(data.frame(expr = exp1, catergory = rep("TP53-Mut PTEN-Mut/Homdel", length(exp1))),
                      data.frame(expr = exp2, catergory = rep("TP53-Mut PTEN-Wt&Diploid", length(exp2))),
                      data.frame( expr = exp3, catergory = rep("TP53-Wt PTEN-Mut/Homdel", length(exp3))),
                      data.frame( expr = exp4, catergory = rep("TP53-Wt PTEN-Wt&Diploid", length(exp4))))

ggplot(df,  aes(x = catergory, y = expr, fill = catergory)) +
  stat_boxplot(geom ='errorbar', outlier.shape = NA, width = 0.3, lwd = 0.2) +
  geom_boxplot(width = 0.6, outlier.shape = NA, lwd = 0.35)+ 
  theme(#panel.border = element_blank(), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line = element_line(colour = "black"),
    #legend.position = c(0.8,0.7),
    legend.position = "None",
    legend.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size=11, color = "black" ),
    axis.text.y = element_text(size=13, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title=element_text(size=14, color = "black"),
    plot.title = element_text(size=14,hjust=0.1, color = "black"),
    #legend.key = element_rect(fill = "lightblue", color = NA),
    # Change legend key size and key width
    legend.key.size = unit(0.25, "cm")
    #strip.text.y  = element_text(size = 26, face="bold")
    ##strip.text is for facet
  )+
  scale_x_discrete(name ="", labels ="") +
  labs(title="",x ="", y = "Expression of CDK4")+
  guides(fill=guide_legend(title="Types of Genetic Alterations"))+
  geom_hline(yintercept=seq(9.5, 10.5, 1), linetype="dashed", color = "grey20")+
  scale_y_continuous(limits = c(8.7,11.2), breaks = seq(8.5, 11.5, 0.5))
ggsave(paste0("CROP_R/4_Transcriptional_interactions/Synergy_", target, "_upregulation_in_PTEN_TP53_patients.pdf"), width = 2.1, height = 2.5)

t.test(exp1, exp2, alternative = "greater") #0.01752
t.test(exp1, exp3, alternative = "greater") #0.001604


#----------------------------------1. Synergy of DNMT1 in patients with NF2-TP53 mutations-----------------------------------------------
target = "DNMT1"

unique(d$YAP1_CNA)
unique(d$TP53_WT_MUT)
#DKO: CDK4 expression in TP53-Mut YAP1-Gain/Amp or WWTR1-Gain/Amp
exp1 = d[d$TP53_WT_MUT %in% "Mutated" & d$YAP1_CNA %in% c("YAP1: Gain", "YAP1: Amplification"), target]
exp2 = d[d$TP53_WT_MUT %in% "Mutated" & d$WWTR1_CNA %in% c("WWTR1: Gain", "WWTR1: Amplification"), target]
#SKO1: CDK4 expression in TP53-Mut YAP1-Diploid/Amp and WWTR1-Diploid
exp3 = d[(d$TP53_WT_MUT %in% "Mutated"&
            d$YAP1_CNA %in% c("YAP1: Diploid") & 
            d$WWTR1_CNA %in% c("WWTR1: Diploid")), target]
#SKO2: CDK4 expression in TP53-WT YAP1-Gain/Amp or WWTR1-Gain/Amp
exp4 = d[d$TP53_WT_MUT %in% "Wild type" & d$YAP1_CNA %in% c("YAP1: Gain", "YAP1: Amplification"), target]
exp5 = d[d$TP53_WT_MUT %in% "Wild type" & d$WWTR1_CNA %in% c("WWTR1: Gain", "WWTR1: Amplification"), target]

#Control: CDK4 expression in TP53-WT YAP1-Diploid/Amp and WWTR1-Diploid
exp6 = d[(d$TP53_WT_MUT %in% "Wild type"&
            d$YAP1_CNA %in% c("YAP1: Diploid") & 
            d$WWTR1_CNA %in% c("WWTR1: Diploid")) , target]

df = rbind.data.frame(data.frame(expr = c(exp1, exp2), catergory = rep("TP53-Mut YAP1/TAZ-Gain/Amp", length(exp1) + length(exp2))),
                      data.frame( expr = exp3, catergory = rep("TP53-Mut YAP1&TAZ-Diploid", length(exp3))),
                      data.frame( expr = c(exp4, exp5), catergory = rep("TP53-Wt YAP1/TAZ-Gain/Amp", length(exp4) + length(exp5))),
                      data.frame( expr = exp6, catergory = rep("TP53-Wt YAP1&TAZ-Diploid", length(exp6))))

ggplot(df,  aes(x = catergory, y = expr, fill = catergory)) +
  stat_boxplot(geom ='errorbar', outlier.shape = NA, width = 0.3, lwd = 0.2) +
  geom_boxplot(width = 0.6, outlier.shape = NA, lwd = 0.35)+ 
  theme(#panel.border = element_blank(), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line = element_line(colour = "black"),
    # legend.position = c(0.8,0.7),
    legend.position = "None",
    legend.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size=11, color = "black" ),
    axis.text.y = element_text(size=13, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title=element_text(size=14, color = "black"),
    plot.title = element_text(size=14,hjust=0.5, color = "black"),
    #legend.key = element_rect(fill = "lightblue", color = NA),
    # Change legend key size and key width
    legend.key.size = unit(0.25, "cm")
    #strip.text.y  = element_text(size = 26, face="bold")
    ##strip.text is for facet
  )+
  scale_x_discrete(name ="", labels ="") +
  labs(title="",x ="", y = paste0("Expression of ", target))+
  guides(fill=guide_legend(title="Types of Genetic Alterations"))+
  #scale_fill_discrete(#palette =  "Pastel2",
  #name="Types of Genetic Alterations",
  #breaks=c("Cells", "Tumors"),
  #labels=c("TP53- YAP/TAZ+", "20-40%", "40-60%", "60-80%", "80-100%"))+
  geom_hline(yintercept=seq(10, 11, 1), linetype="dashed", color = "grey20")+
  scale_y_continuous(limits = c(8.8,12), breaks = seq(8.5, 12, 0.5))
ggsave(paste0("CROP_R/4_Transcriptional_interactions/Synergy_", target, "_upregulation_in_NF2_TP53_patients.pdf"), width = 2.1, height = 2.5)

t.test(c(exp1,exp2), exp3, alternative = "greater") #0.001394
t.test(c(exp1,exp2), c(exp4, exp5), alternative = "greater") #2.38e-12

#----------------------------------2. Synergy of DNMT1 in patients with PTEN-TP53 mutations-----------------------------------------------
unique(d$PTEN_WT_MUT)
unique(d$PTEN_CNA)
#DKO: CDK4 expression in TP53-Mut PTEN-Mut&Homdel
exp1 = d[d$TP53_WT_MUT %in% "Mutated" &
           (d$PTEN_WT_MUT %in% "Mutated" | d$PTEN_CNA %in% "PTEN: Deep Deletion"), target]
#SKO1: CDK4 expression in TP53-Mut PTEN-Wt/Diploid
exp2 = d[d$TP53_WT_MUT %in% "Mutated" &
           d$PTEN_WT_MUT %in% "Wild type" & d$PTEN_CNA %in% "PTEN: Diploid" , target]

#SKO2: CDK4 expression in TP53-Wt PTEN-Mut/Homdel
exp3 = d[d$TP53_WT_MUT %in% "Wild type" & 
           (d$PTEN_WT_MUT %in% "Mutated" | d$PTEN_CNA %in% "PTEN: Deep Deletion"), target]

#Control: CDK4 expression in TP53-Wt PTEN-Wt/Diploid
exp4 = d[d$TP53_WT_MUT %in% "Wild type" &
           d$PTEN_WT_MUT %in% "Wild type" & d$PTEN_CNA %in% "PTEN: Diploid" , target]


df = rbind.data.frame(data.frame(expr = exp1, catergory = rep("TP53-Mut PTEN-Mut/Homdel", length(exp1))),
                      data.frame(expr = exp2, catergory = rep("TP53-Mut PTEN-Wt&Diploid", length(exp2))),
                      data.frame( expr = exp3, catergory = rep("TP53-Wt PTEN-Mut/Homdel", length(exp3))),
                      data.frame( expr = exp4, catergory = rep("TP53-Wt PTEN-Wt&Diploid", length(exp4))))

ggplot(df,  aes(x = catergory, y = expr, fill = catergory)) +
  stat_boxplot(geom ='errorbar', outlier.shape = NA, width = 0.3, lwd = 0.2) +
  geom_boxplot(width = 0.6, outlier.shape = NA, lwd = 0.35)+ 
  theme(#panel.border = element_blank(), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line = element_line(colour = "black"),
    #legend.position = c(0.8,0.7),
    legend.position = "None",
    legend.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size=11, color = "black" ),
    axis.text.y = element_text(size=13, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title=element_text(size=14, color = "black"),
    plot.title = element_text(size=14,hjust=0.1, color = "black"),
    #legend.key = element_rect(fill = "lightblue", color = NA),
    # Change legend key size and key width
    legend.key.size = unit(0.25, "cm")
    #strip.text.y  = element_text(size = 26, face="bold")
    ##strip.text is for facet
  )+
  scale_x_discrete(name ="", labels ="") +
  labs(title="",x ="", y = "Expression of CDK4")+
  guides(fill=guide_legend(title="Types of Genetic Alterations"))+
  geom_hline(yintercept=seq(10, 11, 1), linetype="dashed", color = "grey20")+
  scale_y_continuous(limits = c(8.8,12), breaks = seq(8.5, 12, 0.5))

ggsave(paste0("CROP_R/4_Transcriptional_interactions/Synergy_", target, "_upregulation_in_PTEN_TP53_patients.pdf"), width = 2.1, height = 2.5)

t.test(exp1, exp2, alternative = "greater") #0.02408
t.test(exp1, exp3, alternative = "greater") #0.0002611


#----------------------------------1. Synergy of SRPK1 in patients with NF2-TP53 mutations-----------------------------------------------
target = "SRPK1"

unique(d$YAP1_CNA)
unique(d$TP53_WT_MUT)
#DKO: CDK4 expression in TP53-Mut YAP1-Gain/Amp or WWTR1-Gain/Amp
exp1 = d[d$TP53_WT_MUT %in% "Mutated" & d$YAP1_CNA %in% c("YAP1: Gain", "YAP1: Amplification"), target]
exp2 = d[d$TP53_WT_MUT %in% "Mutated" & d$WWTR1_CNA %in% c("WWTR1: Gain", "WWTR1: Amplification"), target]
#SKO1: CDK4 expression in TP53-Mut YAP1-Diploid/Amp and WWTR1-Diploid
exp3 = d[(d$TP53_WT_MUT %in% "Mutated"&
            d$YAP1_CNA %in% c("YAP1: Diploid") & 
            d$WWTR1_CNA %in% c("WWTR1: Diploid")), target]
#SKO2: CDK4 expression in TP53-WT YAP1-Gain/Amp or WWTR1-Gain/Amp
exp4 = d[d$TP53_WT_MUT %in% "Wild type" & d$YAP1_CNA %in% c("YAP1: Gain", "YAP1: Amplification"), target]
exp5 = d[d$TP53_WT_MUT %in% "Wild type" & d$WWTR1_CNA %in% c("WWTR1: Gain", "WWTR1: Amplification"), target]

#Control: CDK4 expression in TP53-WT YAP1-Diploid/Amp and WWTR1-Diploid
exp6 = d[(d$TP53_WT_MUT %in% "Wild type"&
            d$YAP1_CNA %in% c("YAP1: Diploid") & 
            d$WWTR1_CNA %in% c("WWTR1: Diploid")) , target]

df = rbind.data.frame(data.frame(expr = c(exp1, exp2), catergory = rep("TP53-Mut YAP1/TAZ-Gain/Amp", length(exp1) + length(exp2))),
                      data.frame( expr = exp3, catergory = rep("TP53-Mut YAP1&TAZ-Diploid", length(exp3))),
                      data.frame( expr = c(exp4, exp5), catergory = rep("TP53-Wt YAP1/TAZ-Gain/Amp", length(exp4) + length(exp5))),
                      data.frame( expr = exp6, catergory = rep("TP53-Wt YAP1&TAZ-Diploid", length(exp6))))

ggplot(df,  aes(x = catergory, y = expr, fill = catergory)) +
  stat_boxplot(geom ='errorbar', outlier.shape = NA, width = 0.3, lwd = 0.2) +
  geom_boxplot(width = 0.6, outlier.shape = NA, lwd = 0.35)+ 
  theme(#panel.border = element_blank(), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line = element_line(colour = "black"),
    # legend.position = c(0.8,0.7),
    legend.position = "None",
    legend.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size=11, color = "black" ),
    axis.text.y = element_text(size=13, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title=element_text(size=14, color = "black"),
    plot.title = element_text(size=14,hjust=0.5, color = "black"),
    #legend.key = element_rect(fill = "lightblue", color = NA),
    # Change legend key size and key width
    legend.key.size = unit(0.25, "cm")
    #strip.text.y  = element_text(size = 26, face="bold")
    ##strip.text is for facet
  )+
  scale_x_discrete(name ="", labels ="") +
  labs(title="",x ="", y = paste0("Expression of ", target))+
  guides(fill=guide_legend(title="Types of Genetic Alterations"))+
  #scale_fill_discrete(#palette =  "Pastel2",
  #name="Types of Genetic Alterations",
  #breaks=c("Cells", "Tumors"),
  #labels=c("TP53- YAP/TAZ+", "20-40%", "40-60%", "60-80%", "80-100%"))+
  geom_hline(yintercept=seq(8.5, 9.5, 1), linetype="dashed", color = "grey20")+
  scale_y_continuous(limits = c(7.5,10.5), breaks = seq(7, 11, 0.5))
ggsave(paste0("CROP_R/4_Transcriptional_interactions/Synergy_", target, "_upregulation_in_NF2_TP53_patients.pdf"), width = 2.1, height = 2.5)

t.test(c(exp1,exp2), exp3, alternative = "greater") #8.063e-08
t.test(c(exp1,exp2), c(exp4, exp5), alternative = "greater") #< 2.2e-16


















