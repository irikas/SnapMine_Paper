#Set up Environment
library(tidyverse)
library(ggplot2)
setwd("/Users/irika/Library/CloudStorage/OneDrive-JohnsHopkins/WongLing/Scripts/231220_MultiExonSearchMouse")

#Upload data frame
df <- read.csv("AllCEPSI_Mouse_Meta.csv", row.names = 1)

#Reformat data for visualization
df <- df %>% select("sampleID",ends_with("PSI"))
colnames(df)[2:ncol(df)] = sapply(strsplit(colnames(df)[2:ncol(df)],"__"),"[[",1)
df <- df %>% pivot_longer(cols=!contains("sampleID"),names_to = "Junction",values_to = "PSI") 

#Plot data
df$Junction <- factor(df$Junction, levels=sort(unique(df$Junction),decreasing = T))
p <- ggplot(df, aes(x = Junction, y = PSI)) +  
  geom_jitter(alpha = 0.2, na.rm = T,color="#A364AC",width = 0.2, height = 0,size=1) +
  geom_violin(alpha = 1,width=1, lwd=0.7)+ 
  geom_hline(yintercept=0,linetype=1, color = "black", lwd=1) +
  geom_hline(yintercept=5,linetype=2, color = "#C12126", lwd=0.5) +
  geom_hline(yintercept=15,linetype=2, color = "#C12126", lwd=0.5) +
  theme_classic() + theme(axis.text = element_text(family="Helvetica"))+ 
  scale_y_continuous(limits = c(0, 101), breaks = c(0, 5, 15, 50, 100),expand = c(0, 0))+
  xlab("Gene")+ylab("Average Percent Spliced-In (PSI) of Cryptic Junctions") + coord_flip()
ggsave(plot = p, "231221_srav1m.png",device = "png",width = 30, height = 20, units = "cm",limitsize = FALSE)

#Information about data
df$PSI = as.numeric(df$PSI)
df$category <- cut(df$PSI, 
                   breaks=c(-Inf, 5, 15, Inf), 
                   labels=c("Less than 5%","Between 5% and 15%","Greater than 15%"))
table(df$Junction,df$category)
apply(table(df$Junction,df$category),1,sum)
