library(tidyverse)
library(ggplot2)

#Set Environment
setwd("/Users/irika/Library/CloudStorage/OneDrive-JohnsHopkins/WongLing/Scripts/231220_MultiExonSearch")
info <- read.csv("info.csv")


#Set up functions
regExSamples <- function(df){
  #create df for samples + count
  allSamples <- df$samples
  allSamples<- unlist(strsplit(allSamples,",")) #separate character into string of characters
  allSamples <- allSamples[allSamples!=""] #remove blank samples
  
  #RegEx information: https://www.journaldev.com/36776/regular-expressions-in-r
  sampleID <- gsub(":[^:]+$","",allSamples) #replace everything after and including colon with nothing
  sampleCount <- gsub("^[^:]+:", "",allSamples) #replace everything before and including colon with nothing
  sampleCount <- as.numeric(sampleCount)
  sampleDF <- data.frame(sampleID, sampleCount)
  colnames(sampleDF)=c("sampleID","junctionCount")
  
  return(sampleDF)
}
calcPSI <- function(chr,incA,incB,incC = NA,incD = NA,excA,excD,genome, totalCountMin = 15){
  #Download snaptron junction info around region bounded by exclusion junction
  url <- paste0("https://snaptron.cs.jhu.edu/",genome,"/snaptron?regions=chr",
                chr,":",excA-100,"-",excD + 100)
  options(timeout=300)
  snaptronQuery <- read.csv(url, sep="\t",quote = "")
  
  #Create Count DFs
  if(!is.na(incC) & !is.na(incD)){
    junctions <- list(inc1 = c(incA,incB), 
                      inc2 = c(incC, incD),
                      exc = c(excA, excD))
    
    for(i in 1:3){
      junctionRow <- which(snaptronQuery$start==junctions[[i]][1] & snaptronQuery$end==junctions[[i]][2]) #determine relevant row in Snaptron file
      junctionDF <- snaptronQuery[junctionRow,]
      junctionDF <- regExSamples(junctionDF)
      colnames(junctionDF)[2] = paste0(names(junctions)[i],"_count")
      assign(paste0(names(junctions)[i],"_count"),junctionDF)
    }
    
    countdf <- merge(inc1_count,inc2_count,sort = T, all = T) #include all samples
    countdf <- merge(countdf,exc_count,sort = T, all = T) #include all samples
    countdf[is.na(countdf)] = 0 #any junctions with no counts of a type will be set to 0
    
    #Add PSI values
    countdf <- countdf %>% mutate(totalExonInc1 = inc1_count + exc_count) %>% mutate(PSI_Inc1 = inc1_count/totalExonInc1)
    countdf <- countdf %>% mutate(totalExonInc2 = inc2_count + exc_count) %>% mutate(PSI_Inc2 = inc2_count/totalExonInc2)
    countdf$PSI_Inc1 <- round(countdf$PSI_Inc1*100,2)
    countdf$PSI_Inc2 <- round(countdf$PSI_Inc2*100,2)
    
    #Set PSI to NA if not enough counts of junction to indicate real measurements
    countdf$PSI_Inc1[which(countdf$totalExonInc1 < totalCountMin)] = 0
    countdf$PSI_Inc2[which(countdf$totalExonInc2 < totalCountMin)] = 0
    
    #Generate average
    countdf <- countdf %>% mutate(avgPSI = rowMeans(select(countdf, c(PSI_Inc1,PSI_Inc2)), na.rm = TRUE))
    
    #rename sampleID columns to numeric
    countdf$sampleID <- as.numeric(countdf$sampleID)
    
  } else{
    junctions <- list(inc1 = c(incA,incB), 
                      exc = c(excA, excD))
    
    for(i in 1:2){
      junctionRow <- which(snaptronQuery$start==junctions[[i]][1] & snaptronQuery$end==junctions[[i]][2]) #determine relevant row in Snaptron file
      junctionDF <- snaptronQuery[junctionRow,]
      junctionDF <- regExSamples(junctionDF)
      colnames(junctionDF)[2] = paste0(names(junctions)[i],"_count")
      assign(paste0(names(junctions)[i],"_count"),junctionDF)
    }
    
    countdf <- merge(inc1_count,exc_count,sort = T, all = T) #include all samples
    countdf[is.na(countdf)] = 0 #any junctions with no counts of a type will be set to 0
    
    #Add PSI values
    countdf <- countdf %>% mutate(totalExonInc1 = inc1_count + exc_count) %>% mutate(PSI_Inc1 = inc1_count/totalExonInc1)
    countdf$PSI_Inc1 <- round(countdf$PSI_Inc1*100,2)
    
    #Set PSI to NA if not enough counts of junction to indicate real measurements
    countdf$PSI_Inc1[which(countdf$totalExonInc1 < totalCountMin)] = NA
    countdf$avgPSI = countdf$PSI_Inc1
    
    #rename sampleID columns to numeric
    countdf$sampleID <- as.numeric(countdf$sampleID)
  }
  
  return(countdf)
}
samplesInfo <- function(genome, sampleIDs){
  #Download samples to local folder: "curl https://snaptron.cs.jhu.edu/data/srav1m/samples.tsv >srav1m_samples.tsv"
  fname <- paste0(genome, "_samples.tsv")
  sampleInfo <- read_tsv(fname, lazy = TRUE,quote = "")
  sampleInfo <- sampleInfo %>% select(c("rail_id","external_id","study","study_title",
                                        "library_layout","sample_description","study_abstract",
                                        "study_description","sample_name","sample_title")) %>% filter(rail_id %in% sampleIDs) 
  return(sampleInfo)
}

#Run across info file to populate PSI df
df <- data.frame("sampleID"=matrix(ncol = 1, nrow = 0))

for(i in 1:nrow(info)){
  dftemp <- calcPSI(chr = info$chromosome[i], incA = info$hg38_inc1_A[i], incB = info$hg38_inc1_B[i],
                    incC = info$hg38_inc2_C[i], incD = info$hg38_inc2_D[i], 
                    excA = info$hg38_exc_A[i], excD = info$hg38_exc_D[i], genome = "srav3h")
  dftemp <- dftemp %>% select(sampleID,starts_with("PSI"),ends_with("PSI")) %>% mutate_all(~ifelse(is.nan(.), NA, .))
  colnames(dftemp)[2:ncol(dftemp)] = paste0(info$gene[i],"__",colnames(dftemp)[2:ncol(dftemp)])
  
  df <- merge(df, dftemp,all = T)
  print(i)
  rm(dftemp,i)
}

df <- df %>% mutate(avgbothJunc_rmNA = rowMeans(select(df, ends_with("PSI")), na.rm = TRUE)) %>%
  mutate(avgJunc_rmNA = rowMeans(select(df, contains("PSI_Inc")), na.rm = TRUE))

#Calculate avg junction inclusion across CEs setting NAs as 0s
df2 <- (df %>% mutate_all(~ifelse(is.na(.), 0, .))) 
df2 <- df2 %>% mutate(avgbothJunc_NAas0 = rowMeans(select(df2, ends_with("PSI")))) %>% 
  mutate(avgJunc_NAas0 = rowMeans(select(df2, contains("PSI_Inc"))))
df$avgbothJunc_NAas0 = df2$avgbothJunc_NAas0
df$avgJunc_NAas0 = df2$avgJunc_NAas0
rm(df2)

#Rough normalization of each CE inclusion values and finding avg inclusion across junctions
normdf <- df[,1:(ncol(df)-4)] %>% mutate_all(~ifelse(is.na(.), 0, .))
means = apply(normdf, 2, mean,na.rm=TRUE)
sds = apply(normdf, 2, sd,na.rm=TRUE)

for(i in 2:ncol(normdf)){
  normdf[,i] <- (normdf[,i] - means[i])/sds[i]
}

normdf <- normdf %>% mutate(avgCE_rmNA = rowMeans(select(normdf, ends_with("PSI")), na.rm = TRUE))
df$normAvg = normdf$avgCE_rmNA

#Incorporate metadata and save CSV file for junctions
sampleMeta <- samplesInfo(genome = "srav3h", sampleIDs = df$sampleID)
df <- merge(df,sampleMeta,by.x = "sampleID", by.y = "rail_id")

write.csv(df, "AllCEPSI_Human_Meta.csv")

visualDF <- df %>% select(!contains("PSI_Inc"))
write.csv(visualDF, "AvgJuncCEPSI_Human_Meta.csv")
