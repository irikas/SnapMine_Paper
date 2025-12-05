#' Add metadata to file with samples
#'
#' This function downloads metadata for samples of interest.
#' @param df data.frame that includes Snaptron sampleIDs - column should be named "sampleID". Please subset this in some fashion - otherwise it will take forever to download.
#' @param compilation Compilation of interest (sra_human, sra_mouse, tcga, gtex, encode).
#' @returns A dataframe with five columns ("chromosome", "start", "end", "strand", "samples")
#' @export
#' @references
#' Christopher Wilks, Phani Gaddipati, Abhinav Nellore, Ben Langmead.
#' Snaptron: querying splicing patterns across tens of thousands of RNA-seq samples.
#' Bioinformatics 34, 114-116 (2018). https://doi.org/10.1093/bioinformatics/btx547


add_SnaptronMeta <- function(df, compilation) {
  # Identify which compilation is being used and the appropriate folder name in Snaptron
  keyGenome <- list(
    "sra_human" = "srav3h",
    "sra_mouse" = "srav1m",
    "tcga" = "tcgav2",
    "gtex" = "gtexv2",
    "encode" = "encode1159"
  )

  samples <- unique(df$sampleID)

  ## Download metadata 1000 at a time

  cols_needed <- c(
    "rail_id", "run_acc", "study_acc", "sample_acc", "experiment_acc", "submission_acc",
    "study_title", "study_abstract", "study_description", "experiment_title",
    "sample_title", "sample_name", "sample_description", "sample_attributes",
    "library_layout", "library_construction_protocol"
  )

  ## Samples up to 1000
  # Suppress warnings for fread bc there is some kind of colname mistake >152 but not relevant here
  if (length(samples) <= 1000) {
    urlSamples <- str_remove_all(string = str_flatten_comma(samples), pattern = " ")
    url <- paste0("https://snaptron.cs.jhu.edu/", keyGenome[compilation], "/samples?ids=", urlSamples)
    meta <- suppressWarnings(fread(
      file = url, sep = "\t", quote = "",
      select = cols_needed,
      showProgress = F
    ))
  } else { # >1000 samples
    # Split metadata into lists of 1000
    df_sampleID <- data.frame(sampleID = samples)
    n <- 1000
    split_df <- df_sampleID %>%
      group_by(group = ceiling(row_number() / n)) %>%
      group_split()

    # Create metadata file
    meta <- data.frame()

    # Create progress bar
    groups <- length(split_df)
    pb <- progress_bar$new(
      format = "  Downloading metadata from Snaptron [:bar] :percent Complete // Estimated Time Remaining: :eta",
      total = groups, clear = FALSE, width = round(getOption("width") * 0.75, 0),
      complete = "*"
    )

    for (i in 1:length(split_df)) {
      urlSamples <- str_remove_all(string = str_flatten_comma(unlist(split_df[[i]][1])), pattern = " ")
      url <- paste0("https://snaptron.cs.jhu.edu/", keyGenome[compilation], "/samples?ids=", urlSamples)
      sampleInfo <- suppressWarnings(fread(
        file = url, sep = "\t", quote = "",
        select = cols_needed,
        showProgress = F
      ))
      meta <- rbind(meta, sampleInfo)

      # update progress bar
      pb$tick()
      rm(i, sampleInfo)
    }
  }

  return(meta)
}
