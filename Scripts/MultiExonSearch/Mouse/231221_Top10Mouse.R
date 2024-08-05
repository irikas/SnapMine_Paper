#Set up Environment
library(tidyverse)
library(ggplot2)
setwd("/Users/irika/Library/CloudStorage/OneDrive-JohnsHopkins/WongLing/Scripts/231220_MultiExonSearchMouse")

#Upload data frame
df <- read.csv("AllCEPSI_Mouse_Meta.csv", row.names = 1)

#Reformat data for averages only
df <- df %>% select(!contains("PSI_Inc")) %>% arrange(-avgbothJunc_NAas0)
df <- df[1:20,]

write.csv(df,"Top20Mouse.csv")

#Reformat to Visualize
top20 <- df$sampleID
colnames(df)[2:15] = sapply(strsplit(colnames(df)[2:15],"__"),"[[",1)
df <- df %>% pivot_longer(cols=c(2:15),names_to = "Gene",values_to = "PSI") 

#Visualize
ggplot()
df$sampleID = as.character(df$sampleID)
df$sampleID <- factor(df$sampleID, levels= rev(as.character(top20)))

ggplot(df, aes(Gene, sampleID, fill= PSI)) + geom_tile(colour="white", linewidth=0.25) + 
  scale_fill_gradient2(low="white",mid = "#A364AC", high = "#412234", midpoint=75,na.value = "#d3d3d3",
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
  theme_classic() + 
  theme(axis.text = element_text(family="Helvetica", color = "black"), axis.ticks = element_blank(),
        axis.line = element_blank(),legend.text = element_text(family="Helvetica", color = "black")) +
  scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+
  xlab("Gene") + ylab("Sample") 
ggsave("231221_srav1m_top20.png",device = "png",width = 30, height = 20, units = "cm",limitsize = FALSE)
