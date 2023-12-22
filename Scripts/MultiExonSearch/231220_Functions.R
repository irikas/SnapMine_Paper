#Load required packages
library(tidyverse)

#Define functions

#Name: regExSamples
#Purpose: Reformat samples column of Snaptron junction table in order to create new dataframe with 1 column of rail_ids and 1 column of counts
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


#Name: calcPSI
#Purpose: Given junction inclusion and exclusion coordinates and genome, create a dataframe which includes the inclusion and exclusion junction counts and inclusion junction PSI(s) for each sample with counts of either. Can calculate either one or two inclusion junctions with the same exclusion junction. Requires regExSamples function.
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

#Name: samplesInfo
#Purpose: Use lazy reading to only import sample metadata for the samples identified with inclusion > 15. Requires sampleIDs from output dataframe of calcPSIs but could also just call metadata for a random string of rail_ids if desired.
samplesInfo <- function(genome, sampleIDs){
  #Download samples to local folder: "curl https://snaptron.cs.jhu.edu/data/srav1m/samples.tsv >srav1m_samples.tsv"
  fname <- paste0(genome, "_samples.tsv")
  sampleInfo <- read_tsv(fname, lazy = TRUE,quote = "")
  sampleInfo <- sampleInfo %>% select(c("rail_id","external_id","study","study_title",
                                        "library_layout","sample_description","study_abstract",
                                        "study_description","sample_name","sample_title")) %>% filter(rail_id %in% sampleIDs) 
  return(sampleInfo)
}