#Xiaoyu Zhao 
library(ggplot2)
library(Cairo)
library("ggrepel")
library(readxl)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
setwd("~/Downloads/5_Combinatorial_CROP_seq/")
read0 <- data.frame(read_excel("TableS7.xlsx", 6))
data = read0[,2:ncol(read0)]
rownames(data) = read0[,1]
d = data[,1:9]
#----------------------------------I. Synergy of CDK4 in NF2-TP53 DKO-----------------------------------------------
unique(d$YAP1_CNA)
unique(d$TP53_WT_MUT)
#DKO: CDK4 expression in TP53-Mut YAP1-Gain/Amp or WWTR1-Gain/Amp
exp1 = d[d$TP53_WT_MUT %in% "Mutated" & d$YAP1_CNA %in% c("YAP1: Gain", "YAP1: Amplification"), "CDK4"]
exp2 = d[d$TP53_WT_MUT %in% "Mutated" & d$WWTR1_CNA %in% c("WWTR1: Gain", "WWTR1: Amplification"), "CDK4"]
#SKO1: CDK4 expression in TP53-Mut YAP1-Diploid/Amp and WWTR1-Diploid
exp3 = d[(d$TP53_WT_MUT %in% "Mutated"&
            d$YAP1_CNA %in% c("YAP1: Diploid") & 
            d$WWTR1_CNA %in% c("WWTR1: Diploid")), "CDK4"]
#SKO2: CDK4 expression in TP53-WT YAP1-Gain/Amp or WWTR1-Gain/Amp
exp4 = d[d$TP53_WT_MUT %in% "Wild type" & d$YAP1_CNA %in% c("YAP1: Gain", "YAP1: Amplification"), "CDK4"]
exp5 = d[d$TP53_WT_MUT %in% "Wild type" & d$WWTR1_CNA %in% c("WWTR1: Gain", "WWTR1: Amplification"), "CDK4"]

#Control: CDK4 expression in TP53-WT YAP1-Diploid/Amp and WWTR1-Diploid
exp6 = d[(d$TP53_WT_MUT %in% "Wild type"&
            d$YAP1_CNA %in% c("YAP1: Diploid") & 
            d$WWTR1_CNA %in% c("WWTR1: Diploid")) , "CDK4"]

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
    #legend.position = "None",
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
  labs(title="",x ="", y = "Expression of CDK4")+
  guides(fill=guide_legend(title="Types of Genetic Alterations"))+
  #scale_fill_discrete(#palette =  "Pastel2",
  #name="Types of Genetic Alterations",
  #breaks=c("Cells", "Tumors"),
  #labels=c("TP53- YAP/TAZ+", "20-40%", "40-60%", "60-80%", "80-100%"))+
  geom_hline(yintercept=seq(7.5, 13.5, 1), linetype="dashed", color = "grey20")+
  scale_y_continuous(limits = c(9,11), breaks = seq(9, 11, 0.5))

ggsave("Synergy_CDK4_upregulation_in_NF2_TP53.pdf", width = 5, height = 3)
#ggsave("Synergy_CDK4_upregulation_in_NF2_TP53_2.pdf", width = 2.1, height = 2.5)

t.test(c(exp1,exp2), exp3, alternative = "greater")
t.test(c(exp1,exp2), exp4, alternative = "greater")
t.test(c(exp1,exp2), c(exp5, exp6), alternative = "greater")

#----------------------------------II. Synergy of CDK4 in PTEN-TP53 DKO-----------------------------------------------
unique(d$PTEN_WT_MUT)
unique(d$PTEN_CNA)
#DKO: CDK4 expression in TP53-Mut PTEN-Mut&Homdel
exp1 = d[d$TP53_WT_MUT %in% "Mutated" &
           (d$PTEN_WT_MUT %in% "Mutated" | d$PTEN_CNA %in% "PTEN: Deep Deletion"), "CDK4"]
#SKO1: CDK4 expression in TP53-Mut PTEN-Wt/Diploid
exp2 = d[d$TP53_WT_MUT %in% "Mutated" &
           d$PTEN_WT_MUT %in% "Wild type" & d$PTEN_CNA %in% "PTEN: Diploid" , "CDK4"]

#SKO2: CDK4 expression in TP53-Wt PTEN-Mut/Homdel
exp3 = d[d$TP53_WT_MUT %in% "Wild type" & 
           (d$PTEN_WT_MUT %in% "Mutated" | d$PTEN_CNA %in% "PTEN: Deep Deletion"), "CDK4"]

#Control: CDK4 expression in TP53-Wt PTEN-Wt/Diploid
exp4 = d[d$TP53_WT_MUT %in% "Wild type" &
           d$PTEN_WT_MUT %in% "Wild type" & d$PTEN_CNA %in% "PTEN: Diploid" , "CDK4"]


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
  geom_hline(yintercept=seq(7.5, 13.5, 1), linetype="dashed", color = "grey20")+
  #scale_y_continuous(limits = c(7.5,13.5), breaks = seq(7.5, 13.5, 1))
  scale_y_continuous(limits = c(9,11), breaks = seq(9, 11, 0.5))
#ggsave("Synergy_CDK4_upregulation_in_PTEN_TP53_with_legend.pdf", width = 4.4, height = 3.2)
ggsave("Synergy_CDK4_upregulation_in_PTEN_TP53_2.pdf", width = 2.1, height = 2.5)

t.test(exp1, exp2, alternative = "greater")
t.test(exp1, exp3, alternative = "greater")


#----------------------------------III. Synergy of CDK4 in NF2-PTEN DKO-----------------------------------------------
unique(d$PTEN_WT_MUT)
unique(d$PTEN_CNA)
#DKO: CDK4 expression in YAP1/WWTR-Gain/Amp PTEN-Mut/Homdel
exp1 = d[d$YAP1_CNA %in% c("YAP1: Gain", "YAP1: Amplification") & (d$PTEN_WT_MUT %in% "Mutated" | d$PTEN_CNA %in% "PTEN: Deep Deletion"), "CDK4"]
exp2 = d[d$WWTR1_CNA %in% c("WWTR1: Gain", "WWTR1: Amplification") & (d$PTEN_WT_MUT %in% "Mutated" | d$PTEN_CNA %in% "PTEN: Deep Deletion"), "CDK4"]

#SKO1: CDK4 expression in YAP1/WWTR-Gain/Amp PTEN-Wt&Diploid
exp3 = d[d$YAP1_CNA %in% c("YAP1: Gain", "YAP1: Amplification") & d$PTEN_WT_MUT %in% "Wild type" & d$PTEN_CNA %in% "PTEN: Diploid", "CDK4"]
exp4 = d[d$WWTR1_CNA %in% c("WWTR1: Gain", "WWTR1: Amplification") & d$PTEN_WT_MUT %in% "Wild type" & d$PTEN_CNA %in% "PTEN: Diploid", "CDK4"]

#SKO2: CDK4 expression in YAP1&WWTR-Diploid/Hetloss PTEN-Mut/Homdel
exp5 = d[(d$YAP1_CNA %in% c("YAP1: Diploid") & 
            d$WWTR1_CNA %in% c("WWTR1: Diploid") &
            (d$PTEN_WT_MUT %in% "Mutated" | d$PTEN_CNA %in% "PTEN: Deep Deletion")), "CDK4"]

#Control: CDK4 expression in YAP1&WWTR-Diploid/Hetloss PTEN-Wt&Diploid
exp6 = d[(d$YAP1_CNA %in% c("YAP1: Diploid") & 
            d$WWTR1_CNA %in% c("WWTR1: Diploid")&
            d$PTEN_WT_MUT %in% "Wild type" & d$PTEN_CNA %in% "PTEN: Diploid") , "CDK4"]


df = rbind.data.frame(data.frame( expr = c(exp1, exp2), catergory = rep("YAP1/TAZ-Gain/Amp PTEN-Mut/Homdel", length(exp1) + length(exp2))),
                      data.frame( expr = c(exp3, exp4), catergory = rep("YAP1/TAZ-Gain/Amp PTEN-Wt&Diploid", length(exp3) + length(exp4))),
                      data.frame( expr = exp5, catergory = rep("YAP1&TAZ-Diploid PTEN-Mut/Homdel", length(exp5) )),
                      data.frame( expr = exp6, catergory = rep("YAP1&TAZ-Diploid PTEN-Wt&Diploid", length(exp6))))

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
  geom_hline(yintercept=seq(7.5, 13.5, 1), linetype="dashed", color = "grey20")+
  #scale_y_continuous(limits = c(8.5,11.5), breaks = seq(8.5, 11.5, 1))
  scale_y_continuous(limits = c(9,11), breaks = seq(9, 11, 0.5))
#scale_y_continuous(limits = quantile(df$CDK4pr, c(0.1, 0.9)))
ggsave("Synergy_CDK4_upregulation_in_NF2_PTEN_with_legned.pdf", width = 5, height = 3)
#ggsave("Synergy_CDK4_upregulation_in_NF2_PTEN_2.pdf", width = 2.1, height = 2.5)

t.test(c(exp1,exp2), c(exp3,exp4), alternative = "greater")
t.test(c(exp1,exp2), exp5, alternative = "greater")






















