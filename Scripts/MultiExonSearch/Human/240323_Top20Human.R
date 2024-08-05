#Set up Environment
library(tidyverse)
library(ggplot2)
setwd("/Users/irika/Library/CloudStorage/OneDrive-JohnsHopkins/WongLing/Scripts/240319_GenerateSearch")

#Upload data frame
df <- read.csv("240318_AllCEPSI_Human_Meta.csv", row.names = 1)

#Reformat data for averages only
df <- df %>% select(!contains("PSI_Inc")) %>% arrange(-avgbothJunc_NAas0)
df <- df[1:20,]

write.csv(df,"Top20Human.csv")
df <- read.csv("Top20Human.csv", row.names = 1)

#Reformat to Visualize
top20 <- df$sampleID
colnames(df)[2:27] = sapply(strsplit(colnames(df)[2:27],"__"),"[[",1)
df <- df %>% pivot_longer(cols=c(2:27),names_to = "Gene",values_to = "PSI") %>% 
  filter(Gene %in% c("ACTL6B","AGRN_long","ATG4B_long","ATG4B","G3BP1","MYO18A","PFKP","RANBP1","STMN2","UNC13A","UNC13B"))

df$Gene[which(df$Gene == "AGRN_long")] = "AGRN"
df$Gene[which(df$Gene == "ATG4B_long")] = "ATG4B"

#Visualize
df$sampleID = as.character(df$sampleID)
df$sampleID <- factor(df$sampleID, levels= rev(as.character(top20)))

ggplot(df, aes(Gene, sampleID, fill= PSI)) + geom_tile(colour="white", linewidth=0.25) + 
  scale_fill_gradient2(low="white",mid = "#A364AC", high = "#412234", midpoint=75,na.value = "#d3d3d3",
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
  theme_classic() + 
  theme(axis.text = element_text(family="Helvetica", color = "black"), axis.ticks = element_blank(),
        axis.line = element_blank(),legend.text = element_text(family="Helvetica", color = "black"),
        axis.text.x = element_text(angle = 90)) +
  scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+
  xlab("Gene") + ylab("Sample") 
ggsave("240318_srav3h_top20.png",device = "png",width = 30, height = 20, units = "cm",limitsize = FALSE)
