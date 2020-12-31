#Xiaoyu Zhao
library(ggplot2)
library(RColorBrewer)
library(Cairo)
library(readxl)
setwd("~/Downloads/data/")

#--------------------------I. Examples of growth curve of a combo of DKO, SKOs and CTRL---------------------------------------
#Example: NF2-TP53 in minimal medium
dat = read.csv("minimal/fitness/min.1.RPM_fitlm.csv", row.names = 1, header = T)
#DKO: RPM of NF2-TP53
d1 = dat[grepl("NF2", rownames(dat)) & grepl("TP53", rownames(dat)),1:6]
colnames(d1) = c( "D0", "D3", "D6", "D9", "D12", "D15")
m1 = apply(d1, 2, function(x) mean(x))
se1 = apply(d1, 2, function(x) sd(x)/sqrt(length(x)))
df_1 = cbind.data.frame(Groups = rep("NF2-TP53", 6),Days = c(0,3,6,9,12,15), mean = m1, se = se1 )
#SKO1: RPM of NF2-CTRL
d2 = dat[grepl("NF2", rownames(dat)) & grepl("CTRL", rownames(dat)),1:6]
colnames(d2) = c( "D0", "D3", "D6", "D9", "D12", "D15")
m2 = apply(d2, 2, function(x) mean(x))
se2 = apply(d2, 2, function(x) sd(x)/sqrt(length(x)))
df_2 = cbind.data.frame(Groups = rep("NF2-CTRL", 6),Days = c(0,3,6,9,12,15), mean = m2, se = se2 )
#SKO2: RPM of TP53-CTRL
d3 = dat[grepl("TP53", rownames(dat)) & grepl("CTRL", rownames(dat)),1:6]
colnames(d3) = c( "D0", "D3", "D6", "D9", "D12", "D15")
m3 = apply(d3, 2, function(x) mean(x))
se3 = apply(d3, 2, function(x) sd(x)/sqrt(length(x)))
df_3 = cbind.data.frame(Groups = rep("TP53-CTRL", 6),Days = c(0,3,6,9,12,15), mean = m3, se = se3 )
#Control: RPM of CTRL-CTRL
d4 = dat[grepl("^CTRL", rownames(dat)) & grepl("-CTRL", rownames(dat)),1:6]
colnames(d4) = c( "D0", "D3", "D6", "D9", "D12", "D15")
m4 = apply(d4, 2, function(x) mean(x))
se4 = apply(d4, 2, function(x) sd(x)/sqrt(length(x)))
df_4 = cbind.data.frame(Groups = rep("CTRL-CTRL", 6),Days = c(0,3,6,9,12,15), mean = m4, se = se4 )

#RPM of the combo
df = rbind.data.frame(df_1, df_2, df_3, df_4)
ggplot(data=df, aes(x=Days, y=mean, color=Groups, shape = Groups, linetype = Groups )) +
  geom_line( size = 1)+
  geom_point( size = 1.5)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=2, size = 1, position=position_dodge(0.05)) +
  theme_classic()+
  theme(#panel.border = element_blank(), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #panel.grid= element_line(linetype = "dashed"),
    axis.line = element_line(colour = "black"),
    #legend.position = "none",
    legend.title = element_text(size = 16,color = 'black'),
    legend.text = element_text(size=14 ),
    axis.text.y  =element_text(size=16, color = 'black'),
    axis.text.x =element_text(size=16, color = 'black'),
    axis.title=element_text(size=18, color = 'black'),#face="bold")
    plot.title = element_text(size=18,hjust=0.5),
    strip.text  = element_text(size = 16, color = 'black')
  )+    
  scale_y_continuous(limits = c(0, 1000), breaks = seq(0,1000, by =250))+
  scale_x_continuous(name = "Days",breaks = seq(0,15, by =3),labels = seq(0,15, by= 3))+
  scale_color_manual(values=c("magenta","#FF9500","#007AFF", "black"))+
  labs(title="Minimal",x ="", y = "RPM")
ggsave('GI_example_NF2_TP53/Growth_curve_NF2_TP53_minimal.pdf', width = 4.8, height = 4)

#statistical analysis
t.test(d1[,2]/d1[,1], d2[,2]/d2[,1], alternative = "greater") #p-value = 0.2822
t.test(d1[,3]/d1[,1], d2[,3]/d2[,1], alternative = "greater") #p-value = 0.04412
t.test(d1[,4]/d1[,1], d2[,4]/d2[,1], alternative = "greater") #p-value = 8.597e-05
t.test(d1[,5]/d1[,1], d2[,5]/d2[,1], alternative = "greater") #p-value = 0.001505
t.test(d1[,6]/d1[,1], d2[,6]/d2[,1], alternative = "greater") #p-value = 0.0008468

#-------------------------------------II. cell cycle phase of NF2-TP53 GI---------------------------------------
#DKO: cell cycle phase percentage of NF2-TP53
d1 = data.frame(Group = rep("NF2-TP53", 3), Phase = c("G1/G0","S", "G2/M"), Percentage = c(77.41, 17.07, 5.52))
#SKO1: cell cycle phase percentage of NF2-CTRL
d2 = data.frame(Group = rep("NF2-CTRL", 3), Phase = c("G1/G0","S", "G2/M"), Percentage = c(84.9, 12.68, 2.41))
#SKO2: cell cycle phase percentage of TP53-CTRL
d3 = data.frame(Group = rep("CTRL-TP53", 3), Phase = c("G1/G0","S", "G2/M"), Percentage = c(94.37, 3.24, 2.39))
#Control: cell cycle phase percentage of CTRL-CTRL
d4 = data.frame(Group = rep("CTRL-CTRL", 3), Phase = c("G1/G0","S", "G2/M"), Percentage = c(91.8, 4.23, 3.97))

data = rbind.data.frame(d1, d2, d3, d4)

p <- ggplot(data, aes(x=Group, y=Percentage, fill= factor(Phase, levels = c("G2/M", "S", "G1/G0")))) +
  geom_bar(stat="identity", position = "stack", alpha = 0.9, width = 0.8)+
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        #panel.grid= element_line(linetype = "dashed"),
        axis.line = element_line(colour = "black"),
        #legend.position = c(1,0.8),
        legend.title = element_text(size = 21,color = 'black'),
        legend.text = element_text(size=21,color = 'black'),
        axis.text.y  =element_text(size=23,color = 'black'),
        axis.text.x =element_text(size=20,color = 'black', angle = 45, vjust = 1, hjust = 1),
        axis.title=element_text(size=25,color = 'black')
  )+
  #scale_y_continuous(limits = c(0,103), breaks = seq(0,100,25))+
  labs(title="",x ="", y = "% of cell cycle phases")+
  scale_fill_manual(values=c("red","yellow","blue"),
                    name="Phase",
                    labels=c("G2/M", "S", "G1/G0"))
p
ggsave('GI_example_NF2_TP53/barplot_of_cell_cycle_NF2_TP53_minimal.pdf', width = 4.4, height = 5.5)
