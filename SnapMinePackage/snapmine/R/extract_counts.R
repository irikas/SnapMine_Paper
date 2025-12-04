#' Download Snaptron junction counts based on given junction, strand, and compilation
#'
#' This function creates the Snaptron url for a given novel_junc_id
#' @param junction Junction coordinates to search Snaptron database. Usually "canon_junc_coord" from info file.
#' @param strand Chromosome strand on which the gene of interest is located (+ or -). Usually "strand" from info file.
#' @param compilation Compilation of interest (sra_human, sra_mouse, tcga, gtex, encode),. Usually "compilation" from info file.
#' @returns A dataframe with five columns ("chromosome", "start", "end", "strand", "samples")
#' @export
#' @references
#' Christopher Wilks, Phani Gaddipati, Abhinav Nellore, Ben Langmead.
#' Snaptron: querying splicing patterns across tens of thousands of RNA-seq samples.
#' Bioinformatics 34, 114-116 (2018). https://doi.org/10.1093/bioinformatics/btx547
#'
#' Sinha, I.R., Sandal, P.S., Spence, H. et al.
#' Large-scale RNA-Seq mining reveals ciclopirox olamine induces TDP-43 cryptic exons.
#' Nat Commun 16, 6878 (2025). https://doi.org/10.1038/s41467-025-62004-5


junctionData <- function(junction, strand, compilation){
  # Identify which compilation is being used and the appropriate folder name in Snaptron
  keyGenome = list("sra_human" = "srav3h",
                   "sra_mouse" = "srav1m",
                   "tcga" = "tcgav2",
                   "gtex" = "gtexv2",
                   "encode" = "encode1159")

  # Identify start and end coordinates of query: 100bp +/- query
  chr = str_split_i(junction,":",1)
  start = as.numeric(str_split_i(str_split_i(junction,":",2),"-",1))-100
  end = as.numeric(str_split_i(str_split_i(junction,":",2),"-",2))+100
  url <- paste0("https://snaptron.cs.jhu.edu/",keyGenome[compilation],"/snaptron?regions=",
                chr,":",start,"-",end,"&rfilter=strand:",strand)

  # Download junction information
  options(timeout=300)
  snaptronQuery <- fread(file = url, sep="\t", quote = "",
                         select = c("chromosome", "start", "end", "strand", "samples")) %>%
    mutate(junction = paste0(chromosome,":", start,"-", end)) %>%
    dplyr::select(c(junction, strand, samples))
  return(snaptronQuery)
}

#' Download junction data from Snaptron for all queried junctions
#'
#' This function extracts junction count information from Snaptron.
#' @param uniqueCanonJunc Path to the dataframe containing four columns: junction (one row per each unique canonical junction), strand, compilation, and allJunc (all junctions queried within that canonical junction)
#' @returns A dataframe with three columns (junction, strand, samples)
#' @export
#' @references
#' Christopher Wilks, Phani Gaddipati, Abhinav Nellore, Ben Langmead.
#' Snaptron: querying splicing patterns across tens of thousands of RNA-seq samples.
#' Bioinformatics 34, 114-116 (2018). https://doi.org/10.1093/bioinformatics/btx547
#'
#' Sinha, I.R., Sandal, P.S., Spence, H. et al.
#' Large-scale RNA-Seq mining reveals ciclopirox olamine induces TDP-43 cryptic exons.
#' Nat Commun 16, 6878 (2025). https://doi.org/10.1038/s41467-025-62004-5


snaptronQuery <- function(uniqueCanonJunc){
  # Download junction data from Snaptron
  # Filter to relevant junctions
  # Bind relevant rows to new dataframe
  snaptronQuery_df <- data.frame()

  for(i in 1:nrow(uniqueCanonJunc)){
    df_pull <- junctionData(junction = uniqueCanonJunc$canon_junc_coord[i],
                            strand = uniqueCanonJunc$strand[i],
                            compilation = uniqueCanonJunc$compilation[i])
    snaptronQuery_df <- rbind(snaptronQuery_df, df_pull[which(df_pull$junction %in% c(na.omit(uniqueCanonJunc$allJunc[[i]]))),])
    rm(i)
  }

  return(snaptronQuery_df)
}


#' Extract Junction Counts for each Sample
#'
#' This function extracts junction count information from Snaptron-derived dataframe. It creates a new dataframe with 1 column of rail_ids and 1 column of counts.
#' @param df Path to the input file which is a dataframe from Snaptron.
#' @returns A dataframe with two columns (rail_ids, counts).
#' @export
#' @references
#' Christopher Wilks, Phani Gaddipati, Abhinav Nellore, Ben Langmead.
#' Snaptron: querying splicing patterns across tens of thousands of RNA-seq samples.
#' Bioinformatics 34, 114-116 (2018). https://doi.org/10.1093/bioinformatics/btx547
#'
#' Sinha, I.R., Sandal, P.S., Spence, H. et al.
#' Large-scale RNA-Seq mining reveals ciclopirox olamine induces TDP-43 cryptic exons.
#' Nat Commun 16, 6878 (2025). https://doi.org/10.1038/s41467-025-62004-5

sampleJunctionCounts <- function(df){
  # Create df for samples + count
  # Separate character into string of characters
  # Remove blanks
  allSamples <- str_split_1(df$samples, ",")
  allSamples <- allSamples[allSamples!=""]

  # Split sample IDs and counts
  sampleDF <- data.frame(sampleID = str_split_i(allSamples,":",1),
                         junctionCount = as.numeric(str_split_i(allSamples,":",2)))

  return(sampleDF)
}
