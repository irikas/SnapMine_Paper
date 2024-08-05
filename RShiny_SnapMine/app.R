#Load required packages
library(shiny)
library(DT)
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggthemr)
library(shinyWidgets)
library(shinybusy)


ggthemr('dust')

#Load Functions
regExSamples_splice <- function(df){
  #create df for samples + count
  allSamples <- df$samples
  allSamples<- unlist(strsplit(allSamples,",")) #separate character into striling of characters
  allSamples <- allSamples[allSamples!=""] #remove blank samples
  
  #RegEx information: https://www.journaldev.com/36776/regular-expressions-in-r
  sampleID <- gsub(":[^:]+$","",allSamples) #replace everything after and including colon with nothing
  sampleCount <- gsub("^[^:]+:", "",allSamples) #replace everything before and including colon with nothing
  sampleCount <- as.numeric(sampleCount)
  sampleDF <- data.frame(sampleID, sampleCount)
  colnames(sampleDF)=c("sampleID","junctionCount")
  
  return(sampleDF)
}
regExSamples_gex <- function(df){
  #create df for samples + count
  allSamples <- df$samples
  allSamples<- unlist(strsplit(allSamples,",")) #separate character into string of characters
  allSamples <- allSamples[allSamples!=""] #remove blank samples
  
  #RegEx information: https://www.journaldev.com/36776/regular-expressions-in-r
  sampleID <- gsub(":[^:]+$","",allSamples) #replace everything after and including colon with nothing
  sampleCount <- gsub("^[^:]+:", "",allSamples) #replace everything before and including colon with nothing
  sampleCount <- as.numeric(sampleCount)
  sampleDF <- data.frame(sampleID, sampleCount)
  colnames(sampleDF)=c("sampleID","gexCount")
  
  return(sampleDF)
}
calcPSI <- function(chr,incA,incB,incC = NA,incD = NA,excA,excD,genome, totalCountMin = 15){
  #Download snaptron junction info around region bounded by exclusion junction
  url <- paste0("https://snaptron.cs.jhu.edu/",genome,"/snaptron?regions=",
                chr,":",excA-100,"-",excD + 100)
  options(timeout=300)
  #snaptronQuery <- read.csv(url, sep="\t",quote = "")
  snaptronQuery <- fread(url,sep="\t",quote = "")
  
  #Create Count DFs
  if(!is.na(incC) & !is.na(incD)){
    junctions <- list(inc1 = c(incA,incB), 
                      inc2 = c(incC, incD),
                      exc = c(excA, excD))
    
    for(i in 1:3){
      junctionRow <- which(snaptronQuery$start==junctions[[i]][1] & snaptronQuery$end==junctions[[i]][2]) #determine relevant row in Snaptron file
      junctionDF <- snaptronQuery[junctionRow,]
      junctionDF <- regExSamples_splice(junctionDF)
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
    countdf$avgPSI = round(countdf$avgPSI,2)
    
    #rename sampleID columns to numeric
    countdf$sampleID <- as.numeric(countdf$sampleID)
    
  } else{
    junctions <- list(inc1 = c(incA,incB), 
                      exc = c(excA, excD))
    
    for(i in 1:2){
      junctionRow <- which(snaptronQuery$start==junctions[[i]][1] & snaptronQuery$end==junctions[[i]][2]) #determine relevant row in Snaptron file
      junctionDF <- snaptronQuery[junctionRow,]
      junctionDF <- regExSamples_splice(junctionDF)
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
calcGEX <- function(genome,geneName){
  # Download gene counts
  url = paste0("https://snaptron.cs.jhu.edu/",genome,"/genes?regions=", 
               geneName)
  options(timeout=300)
  snaptronQuery <- fread(url,sep="\t",quote = "")
  snaptronQuery <- snaptronQuery[which(toupper(str_split_i(snaptronQuery$`gene_id:gene_name:gene_type:bp_length`,":",2))==toupper(geneName)),]
  #snaptronQuery <- snaptronQuery[which(grepl(geneName,str_split_i(snaptronQuery$`gene_id:gene_name:gene_type:bp_length`,":",2),ignore.case = T)),] 
  
  #Create sample df
  countdf <- regExSamples_gex(snaptronQuery)
  return(countdf)
}
samplesInfo <- function(genome, sampleIDs){
  #Download samples to local folder: "curl https://snaptron.cs.jhu.edu/data/srav1m/samples.tsv >srav1m_samples.tsv"
  
  if(genome == "mouse (SRA)"){
    sampleInfo <- read_tsv("srav1m_samples.tsv", lazy = TRUE,quote = "")
    sampleInfo <- sampleInfo %>% select(c("rail_id","external_id","study","study_title",
                                          "library_layout","sample_description", 
                                          #"study_abstract","study_description",
                                          "sample_name","sample_title")) %>% filter(rail_id %in% sampleIDs)
  } else if(genome == "human (SRA)"){
    sampleInfo <- read_tsv("srav3h_samples.tsv", lazy = TRUE,quote = "")
    sampleInfo <- sampleInfo %>% select(c("rail_id","external_id","study","study_title",
                                          "library_layout","sample_description",
                                          #"study_abstract","study_description",
                                          "sample_name","sample_title")) %>% filter(rail_id %in% sampleIDs)
  } else if(genome == "human (TCGA)"){
    sampleInfo <- read_tsv("tcgav2_samples.tsv", lazy = TRUE,quote = "")
    sampleInfo <- sampleInfo %>% select(c("rail_id","tcga_barcode", "study","gdc_cases.project.name",
                                          "gdc_cases.project.primary_site", "cgc_sample_sample_type",
                                          "gdc_state","gdc_cases.demographic.race",
                                          "gdc_cases.demographic.ethnicity",
                                          "gdc_cases.diagnoses.tumor_stage",
                                          "gdc_cases.diagnoses.vital_status","gdc_cases.samples.oct_embedded",
                                          "gdc_cases.samples.is_ffpe","gdc_cases.samples.sample_type",
                                          "cgc_sample_country_of_sample_procurement","cgc_case_tumor_status",
                                          "cgc_drug_therapy_pharmaceutical_therapy_type","cgc_follow_up_tumor_status")) %>% filter(rail_id %in% sampleIDs)
    
  } else if(genome == "human (GTEX)"){
    sampleInfo <- read_tsv("gtexv2_samples.tsv", lazy = TRUE,quote = "")
    sampleInfo <- sampleInfo %>% select(c("rail_id","run_acc", "study","SEX","AGE", 
                                          "SAMPID","SMTS","SMTSD")) %>% filter(rail_id %in% sampleIDs)
  }
  
  return(sampleInfo)
}

# Parameter tabs
parameter_tabs <- tabsetPanel(
  id = "params_splice",
  type = "hidden",
  tabPanel("onejunction",
           textInput("inc1", "Inclusion Junction 1", value = "chr19:4491836-4492014"),
           textInput("exc", "Exclusion Junction", value = "chr19:4491836-4493702")
  ),
  tabPanel("twojunctions",
           textInput("inc1_twojunc", "Inclusion Junction 1", value = "chr19:4491836-4492014"),
           textInput("inc2_twojunc", "Inclusion Junction 2", value = "chr19:4492153-4493702"),
           textInput("exc_twojunc", "Exclusion Junction", value = "chr19:4491836-4493702")
  )
)

parameter_GEX_tabs <- tabsetPanel(
  id = "params_gex",
  type = "hidden",
  tabPanel("rawdata",
           textInput("geneName_GEX", "Gene", value = "TARDBP")
  ),
  tabPanel("normalized",
           HTML(r"(<h4>Query Information</h4>)"),
           textInput("geneName_GEX_qNorm", "Gene", value = "TARDBP"),
           br(),
           HTML(r"(<h4>Normalization Information</h4>)"),
           textInput("geneName_GEX_Norm", "Gene", value = "EDF1"),
  )
)

# UI
ui <- navbarPage(
  theme = bslib::bs_theme(bootswatch = "sandstone"),
  
  "SnapMine",
  tabPanel("Splice Query",
           p(HTML(r"(<b>07/19/2024 Update: This website has moved to <a href="https://snapmine.idies.jhu.edu/">snapmine.idies.jhu.edu</a>. Please use that version of the website!  - Irika</b>)"),style="color:#A364AC"),
           p(HTML(r"(<b>Contact info:</b> <br> Irika - isinha1 [@] jhu [.] edu <br> Dr. Jonathan Ling - jling [@] jhu [.] edu)"),style="color:#000000"),
           
           ## Sidebar
           sidebarLayout(
             sidebarPanel(
               selectInput(inputId = "genome_splice", label = "Select the organism of interest", 
                           choices = c("human (SRA)", "human (TCGA)", "human (GTEX)", "mouse (SRA)")),
               hr(), 
               HTML(r"(<h3>Junction Information</h3>)"),
               
               fixedRow(column(5,align="center","one incl. junction"),
                        column(2,materialSwitch(inputId = "junctionNum_splice", 
                                                status = "primary")),
                        column(5,align="center","two incl. junctions")),

               #materialSwitch(inputId = "junctionNum_splice", label = "two inc junctions", 
                              #                status = "primary", right = TRUE),
               hr(),
               parameter_tabs,
               actionButton("calculate_splice", "Calculate PSI")
             ),
             
             ## Main panel
             mainPanel(
               tabsetPanel(
                 tabPanel("Information",imageOutput("photo"),p(HTML(r"(<i>Note: The default values are set for the HDGFL2 cryptic exon in human</i>)"),style="color:#000000; font-size:15px")),
                 tabPanel("Results", fluidRow(splitLayout(cellWidths = c("50%", "50%"),
                                                          plotOutput("plot1_splice"), 
                                                          plotOutput("plot2_splice"))))
               )
             )
           ),
           br(),
           "Wait to click download until the table shows up please",
           downloadButton("download_splice"),
           add_busy_spinner(spin = "fading-circle",position="full-page",timeout = 2,color = "#002d72"),
           dataTableOutput("static_splice"),
           p(HTML(r"(<b>Built using <a href=https://snaptron.cs.jhu.edu/ target="_blank" rel="noopener noreferrer" style="color:#002d72 ">Snaptron</a> and <a href=https://rna.recount.bio/ target="_blank" rel="noopener noreferrer" style="color:#002d72">recount3</a></b>)"))
           ),
  
  ### GEX QUERY
  tabPanel("GEX Query",
           p(HTML(r"(<b>This webpage is still a work in progress and will be updated/replaced as we work on it. Thanks for the support! - Irika</b>)"),style="color:#A364AC"),
           p(HTML(r"(<b>Contact info:</b> <br> Irika - isinha1 [@] jhu [.] edu <br> Dr. Jonathan Ling - jling [@] jhu [.] edu)"),style="color:#A364AC"),
           
           ## Sidebar
           sidebarLayout(
             sidebarPanel(
               selectInput(inputId = "genome_gex", label = "Select the organism of interest", 
                           choices = c("human (SRA)", "human (TCGA)", "human (GTEX)", "mouse (SRA)")),
               hr(), 
               HTML(r"(<h3>Gene Information</h3>)"),
               materialSwitch(inputId = "normalizedOrNo", label = "normalized GEX", 
                              status = "primary", right = TRUE),
               hr(),
               parameter_GEX_tabs,
               actionButton("calculate_gex", "Calculate GEX")
             ),
             
             ## Main panel
             mainPanel(
               tabsetPanel(
                 tabPanel("Information",imageOutput("photo_gex"),
                          p("Gene expression counts can be normalized within the study group to a specified gene.", style = "font-family: 'helvetica'; font-si12pt")),
                 tabPanel("Results", fluidRow(splitLayout(cellWidths = c("50%", "50%"),
                                                          plotOutput("plot1_gex"), 
                                                          plotOutput("plot2_gex"))))
               )
             )
           ),
           br(),
           "Wait to click download until the table shows up please",
           downloadButton("download_gex"),
           add_busy_spinner(spin = "fading-circle",position="full-page",timeout = 2,color = "#002d72"),
           dataTableOutput("static_gex")),
)

#Server
server <- function(input, output, session) {
  
  # SPLICE QUERY
  observeEvent(input$junctionNum_splice, {
    
    if(input$junctionNum_splice == TRUE){
      updateTabsetPanel(inputId = "params_splice",selected = "twojunctions") 
    }else{
      updateTabsetPanel(inputId = "params_splice",selected = "onejunction")
    }})
  
  output$photo <- renderImage({
    list(
      src = file.path("2301018_SnaptronQuery.png"),
      contentType = "image/png",
      height = "100%"
    )
  }, deleteFile = FALSE)
  
  output$static <- DT::renderDataTable(NULL)
  
  observeEvent(input$calculate_splice, {
    
    if(input$junctionNum_splice == TRUE){
      #browser()
      #Input values 
      inc1 <- gsub(" ","",input$inc1_twojunc)
      inc1<- gsub(",","",inc1)
      
      inc2 <- gsub(" ","",input$inc2_twojunc)
      inc2 <- gsub(",","",inc2)
      
      exc <- gsub(" ","",input$exc_twojunc)
      exc <- gsub(",","",exc)
      
      coordinates <- str_split_i(c(inc1,
                                   inc2,
                                   exc),
                                 c(":"),2)
      chr <- str_split_1(gsub(" ","",input$inc1_twojunc),c(":"))[1]
      incA <- as.numeric(str_split_1(string = coordinates[1],pattern = "-")[1])
      incB <- as.numeric(str_split_1(string = coordinates[1],pattern = "-")[2])
      incC <- as.numeric(str_split_1(string = coordinates[2],pattern = "-")[1])
      incD <- as.numeric(str_split_1(string = coordinates[2],pattern = "-")[2])
      excA <- as.numeric(str_split_1(string = coordinates[3],pattern = "-")[1])
      excD <- as.numeric(str_split_1(string = coordinates[3],pattern = "-")[2])
      
      if(input$genome_splice == "mouse (SRA)"){
        finalDF <- calcPSI(chr = chr, incA = incA, incB = incB,
                           incC = incC, incD = incD,
                           excA = excA, excD = excD, genome = "srav1m")
        finalDF <- finalDF %>% select(sampleID,starts_with("PSI"),ends_with("PSI")) %>% mutate_all(~ifelse(is.nan(.), NA, .))
        
      } else if(input$genome_splice == "human (SRA)"){
        finalDF <- calcPSI(chr = chr, incA = incA, incB = incB,
                           incC = incC, incD = incD,
                           excA = excA, excD = excD, genome = "srav3h")
        finalDF <- finalDF %>% select(sampleID,starts_with("PSI"),ends_with("PSI")) %>% mutate_all(~ifelse(is.nan(.), NA, .))
        
      } else if(input$genome_splice == "human (TCGA)"){
        finalDF <- calcPSI(chr = chr, incA = incA, incB = incB,
                           incC = incC, incD = incD,
                           excA = excA, excD = excD, genome = "tcgav2")
        finalDF <- finalDF %>% select(sampleID,starts_with("PSI"),ends_with("PSI")) %>% mutate_all(~ifelse(is.nan(.), NA, .))
        
      } else if(input$genome_splice == "human (GTEX)"){
        finalDF <- calcPSI(chr = chr, incA = incA, incB = incB,
                           incC = incC, incD = incD,
                           excA = excA, excD = excD, genome = "gtexv2")
        finalDF <- finalDF %>% select(sampleID,starts_with("PSI"),ends_with("PSI")) %>% mutate_all(~ifelse(is.nan(.), NA, .))
      }
    }else{
      #Input values 
      inc1 <- gsub(" ","",input$inc1)
      inc1 <- gsub(",","",inc1)
      
      exc <- gsub(" ","",input$exc)
      exc <- gsub(",","",exc)
      
      coordinates <- str_split_i(c(inc1,exc),
                                 c(":"),2)
      chr <- str_split_1(gsub(" ","",input$inc1),c(":"))[1]
      incA <- as.numeric(str_split_1(coordinates[1],"-")[1])
      incB <- as.numeric(str_split_1(coordinates[1],"-")[2])
      excA <- as.numeric(str_split_1(coordinates[2],"-")[1])
      excD <- as.numeric(str_split_1(coordinates[2],"-")[2])
      
      if(input$genome_splice == "mouse (SRA)"){
        finalDF <- calcPSI(chr = chr, incA = incA, incB = incB,
                           excA = excA, excD = excD, genome = "srav1m")
        finalDF <- finalDF %>% select(sampleID,starts_with("PSI"),ends_with("PSI")) %>% mutate_all(~ifelse(is.nan(.), NA, .))
        
      } else if(input$genome_splice == "human (SRA)"){
        finalDF <- calcPSI(chr = chr, incA = incA, incB = incB,
                           excA = excA, excD = excD, genome = "srav3h")
        finalDF <- finalDF %>% select(sampleID,starts_with("PSI"),ends_with("PSI")) %>% mutate_all(~ifelse(is.nan(.), NA, .))
        
      } else if(input$genome_splice == "human (TCGA)"){
        finalDF <- calcPSI(chr = chr, incA = incA, incB = incB,
                           excA = excA, excD = excD, genome = "tcgav2")
        finalDF <- finalDF %>% select(sampleID,starts_with("PSI"),ends_with("PSI")) %>% mutate_all(~ifelse(is.nan(.), NA, .))
        
      } else if(input$genome_splice == "human (GTEX)"){
        finalDF <- calcPSI(chr = chr, incA = incA, incB = incB,
                           excA = excA, excD = excD, genome = "gtexv2")
        finalDF <- finalDF %>% select(sampleID,starts_with("PSI"),ends_with("PSI")) %>% mutate_all(~ifelse(is.nan(.), NA, .))
      }
    }
    
    sampleMeta <- samplesInfo(genome = input$genome_splice, sampleIDs = finalDF$sampleID)
    finalDF <- merge(finalDF,sampleMeta,by.x = "sampleID", by.y = "rail_id")
    
    output$static_splice <- DT::renderDataTable(finalDF)
    
    output$plot1_splice <-  renderPlot(
      ggplot(finalDF %>% dplyr::filter(avgPSI >= 5), aes(x=factor(0),avgPSI))+geom_boxplot() +
        theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) + xlab("")+
        ggtitle("Overall Distribution of Junction Inclusion")+ylab("PSI of Inclusion Junction"),
      res = 96
    )
    
    output$plot2_splice <-  renderPlot(
      ggplot(finalDF %>% dplyr::filter(avgPSI >= 5), aes(x=avgPSI))+geom_histogram(binwidth = 3) +
        ggtitle("Inclusion Junction PSI Distribution",subtitle = "avgPSI >= 5")+ xlab("PSI")+xlim(0,100),
      res = 96
    )
    
    output$download_splice <- downloadHandler(
      filename = function() {
        paste0(format(Sys.time(),"%Y%m%d_%H%M%S"),"JunctionInclusion", ".csv")
      },
      content = function(file) {
        write.csv(finalDF, file)
      }
    )
  })
  
  ### GEX QUERY
  
  output$photo_gex <- renderImage({
    list(
      src = file.path("240318_GEX.png"),
      contentType = "image/png",
      height = "100%"
    )
  }, deleteFile = FALSE)
  
  observeEvent(input$normalizedOrNo, {
    
    if(input$normalizedOrNo == TRUE){
      updateTabsetPanel(inputId = "params_gex",selected = "normalized") 
    }else{
      updateTabsetPanel(inputId = "params_gex",selected = "rawdata")
    }})
  
  output$static_gex <- DT::renderDataTable(NULL)
  
  observeEvent(input$calculate_gex, {
    
    if(input$normalizedOrNo == TRUE){
      #Input values 
      
      if(input$genome_gex == "mouse (SRA)"){
        queryDF <- calcGEX(geneName = input$geneName_GEX_qNorm,genome = "srav1m")
        normDF <- calcGEX(geneName = input$geneName_GEX_Norm,genome = "srav1m")
        colnames(normDF)[2] <- "gexCount_Norm"
        
      } else if(input$genome_gex == "human (SRA)"){
        queryDF <- calcGEX(geneName = input$geneName_GEX_qNorm,genome = "srav3h")
        normDF <- calcGEX(geneName = input$geneName_GEX_Norm,genome = "srav3h")
        colnames(normDF)[2] <- "gexCount_Norm"
        
      } else if(input$genome_gex == "human (TCGA)"){
        queryDF <- calcGEX(geneName = input$geneName_GEX_qNorm,genome = "tcgav2")
        normDF <- calcGEX(geneName = input$geneName_GEX_Norm,genome = "tcgav2")
        colnames(normDF)[2] <- "gexCount_Norm"
        
      } else if(input$genome_gex == "human (GTEX)"){
        queryDF <- calcGEX(geneName = input$geneName_GEX_qNorm,genome = "gtexv2")
        normDF <- calcGEX(geneName = input$geneName_GEX_Norm,genome = "gtexv2")
        colnames(normDF)[2] <- "gexCount_Norm"
        
      }
      
      # Create final DF for samples
      finalDF <-  full_join(queryDF,normDF,by=c("sampleID"="sampleID"))
      sampleMeta <- samplesInfo(genome = input$genome_gex, sampleIDs = finalDF$sampleID)
      
      # Determine the normalization factor for each study group and divide sample values by normalization factor
      normDF <- merge(normDF,sampleMeta[,1:3],by.x = "sampleID", by.y = "rail_id")
      normDFstats <- normDF %>% group_by(study) %>% summarize(max = max(gexCount_Norm))
      normDF <- left_join(normDF,normDFstats,by=c("study"="study"))
      normDF <- normDF %>% mutate(NormFactor = gexCount_Norm/max)
      dfnormNum <- normDF %>% select(c("sampleID","NormFactor"))
      
      finalDF <-  left_join(queryDF,dfnormNum,by=c("sampleID"="sampleID"))
      finalDF <- finalDF %>% mutate(normCount = gexCount/NormFactor)
      
      # Create final df for plotting
      finalDF <- merge(finalDF,sampleMeta,by.x = "sampleID", by.y = "rail_id")
      finalDFnoNorm <- finalDF 
      finalDF <- finalDF %>% filter(!is.na(normCount))
      
      output$static_gex <- DT::renderDataTable(finalDF)
      
      output$plot1_gex <-  renderPlot(
        ggplot(finalDF, aes(x=factor(0),log2(normCount)))+geom_boxplot() +
          theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.title = element_text(hjust = 0.5)) + xlab("")+
          ggtitle("Overall Distribution GEX")+ylab("log2(Normalized GEX count)"),
        res = 96
      )
      
      output$plot2_gex <-  renderPlot(
        ggplot(finalDF, aes(x=factor(0),log2(gexCount)))+geom_boxplot() +
          theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.title = element_text(hjust = 0.5)) + xlab("")+
          ggtitle("Overall Distribution of GEX (Not Normalized)")+ylab("log2(GEX count) Not Normalized"),
        res = 96
      )
      
      output$download_gex <- downloadHandler(
        filename = function() {
          paste0(format(Sys.time(),"%Y%m%d_%H%M%S"),"_GEX", ".csv")
        },
        content = function(file) {
          write.csv(finalDFnoNorm, file)
        }
      )
      
    }else{
      
      if(input$genome_gex == "mouse (SRA)"){
        finalDF <- calcGEX(input$geneName_GEX,genome = "srav1m")
        
      } else if(input$genome_gex == "human (SRA)"){
        finalDF <- calcGEX(input$geneName_GEX,genome = "srav3h")
        
      } else if(input$genome_gex == "human (TCGA)"){
        finalDF <- calcGEX(input$geneName_GEX,genome = "tcgav2")
        
      } else if(input$genome_gex == "human (GTEX)"){
        finalDF <- calcGEX(input$geneName_GEX,genome = "gtexv2")
        
      }
      
      sampleMeta <- samplesInfo(genome = input$genome_gex, sampleIDs = finalDF$sampleID)
      finalDF <- merge(finalDF,sampleMeta,by.x = "sampleID", by.y = "rail_id")
      
      output$static_gex <- DT::renderDataTable(finalDF)
      
      output$plot1_gex <-  renderPlot(
        ggplot(finalDF, aes(x=factor(0),log2(gexCount)))+geom_boxplot() +
          theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.title = element_text(hjust = 0.5)) + xlab("")+
          ggtitle("Overall Distribution of GEX")+ylab("log2(GEX count - Not Normalized)"),
        res = 96
      )
      
      output$download_gex <- downloadHandler(
        filename = function() {
          paste0(format(Sys.time(),"%Y%m%d_%H%M%S"),"_GEX", ".csv")
        },
        content = function(file) {
          write.csv(finalDF, file)
        }
      )
      
    }
    
    
  })
 
}

shinyApp(ui, server)