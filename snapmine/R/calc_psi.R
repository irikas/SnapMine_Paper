#' Calculate PSI values for a single junction of interest
#'
#' This function calculates PSI values for cryptic exons given junction coordinates
#' @param novel_junc_coord Coordinates of the novel junction of interest. Given as chr(1-22XYM):(start)-(end).
#' @param canon_junc_coord Coordinates of the related canonical junction Given as chr(1-22XYM):(start)-(end).
#' @param snaptron_df Output from snaptronQuery function.
#' @param totalCountMin Sum of reads of canonical and novel junctions. Default: 15.
#' @param filter_noCanon Filter out samples with 0 reads of canonical junction. Default: T.
#' @returns PSI value of given junction of interest within all Snaptron samples
#' @export
#' @references
#' Christopher Wilks, Phani Gaddipati, Abhinav Nellore, Ben Langmead.
#' Snaptron: querying splicing patterns across tens of thousands of RNA-seq samples.
#' Bioinformatics 34, 114-116 (2018). https://doi.org/10.1093/bioinformatics/btx547
#'
#' Sinha, I.R., Sandal, P.S., Spence, H. et al.
#' Large-scale RNA-Seq mining reveals ciclopirox olamine induces TDP-43 cryptic exons.
#' Nat Commun 16, 6878 (2025). https://doi.org/10.1038/s41467-025-62004-5

calculatePSI <- function(novel_junc_coord, canon_junc_coord, snaptron_df, totalCountMin = 15, filter_noCanon = T) {
  if (!is.na(novel_junc_coord)) {
    # Find rows of interest
    canonRow <- which(snaptron_df$junction == canon_junc_coord)
    novelRow <- which(snaptron_df$junction == novel_junc_coord)

    # Subset to rows of interest if they exist.
    # If no relevant row - make empty dataframe
    # Parse the data to generate PSIs of interest
    if (!is_empty(canonRow)) {
      canonRow <- snaptron_df[canonRow, ]
      canonCounts <- sampleJunctionCounts(canonRow)
    } else {
      canonRow <- data.frame()
    }

    if (!is_empty(novelRow)) {
      novelRow <- snaptron_df[novelRow, ]
      novelCounts <- sampleJunctionCounts(novelRow)
    } else {
      novelRow <- data.frame()
    }

    # Bind counts together by sampleID & filter based on minimum set
    counts_df <- canonCounts %>%
      full_join(novelCounts, by = "sampleID", suffix = c("_canon", "_novel")) %>%
      mutate(
        junctionCount_canon = ifelse(is.na(junctionCount_canon), 0, junctionCount_canon),
        junctionCount_novel = ifelse(is.na(junctionCount_novel), 0, junctionCount_novel),
        totalCounts = junctionCount_canon + junctionCount_novel
      ) %>%
      filter(totalCounts >= totalCountMin)

    # Remove samples in which there are no counts of the canonical junction if set filter_noCanon = True
    if (filter_noCanon) {
      counts_df <- counts_df %>% filter(junctionCount_canon != 0)
    }

    # Calculate PSI
    counts_df <- counts_df %>%
      mutate(
        psi = junctionCount_novel / (totalCounts),
        psi = round(psi * 100, 2)
      ) %>%
      dplyr::select(sampleID, psi)

    return(counts_df)
  }
}

#' Calculate PSI values for each junction of interest from an info file
#'
#' This function calculates PSI values for cryptic exons given an info file
#' @param info_df Path to the info file which contains (at least) the following columns: novel_junc_id, compilation (sra_human, sra_mouse, tcga, gtex, encode), strand, novel_junc_left_coord, novel_junc_right_coord, canon_junc_coord.
#' @returns List of daataframes - one per each unique target in info_df.
#' @export
#' @references
#' Christopher Wilks, Phani Gaddipati, Abhinav Nellore, Ben Langmead.
#' Snaptron: querying splicing patterns across tens of thousands of RNA-seq samples.
#' Bioinformatics 34, 114-116 (2018). https://doi.org/10.1093/bioinformatics/btx547
#'
#' Sinha, I.R., Sandal, P.S., Spence, H. et al.
#' Large-scale RNA-Seq mining reveals ciclopirox olamine induces TDP-43 cryptic exons.
#' Nat Commun 16, 6878 (2025). https://doi.org/10.1038/s41467-025-62004-5

calcPSI_bulk <- function(info_df) {
  # Identify unique canonical junctions to minimize duplication of Snaptron queries
  uniqueCanonJunc <- info_df %>% nest(allJunc = c(canon_junc_coord, novel_junc_left_coord, novel_junc_right_coord), .by = c(canon_junc_coord, strand, compilation))
  uniqueCanonJunc$allJunc <- lapply(uniqueCanonJunc$allJunc, function(x) unique(unlist(x)))

  # Run Snaptron query to download all queried junctions
  snaptronQuery_df <- snaptronQuery(uniqueCanonJunc)

  # For each target in info file generate counts file and merge (left junction)
  counts_list_left <- map2(info_df$canon_junc_coord, info_df$novel_junc_left_coord, function(x, y) calculatePSI(y, x, snaptronQuery_df) %>% as.data.frame())
  names(counts_list_left) <- info_df$novel_junc_id

  # For each target in info file generate counts file and merge (right junction)
  counts_list_right <- map2(info_df$canon_junc_coord, info_df$novel_junc_right_coord, function(x, y) calculatePSI(y, x, snaptronQuery_df))
  names(counts_list_right) <- info_df$novel_junc_id

  # For NAs - replace with other side
  counts_list_right[which(is.na(info_df$novel_junc_right_coord))] <- counts_list_left[which(is.na(info_df$novel_junc_right_coord))]
  counts_list_left[which(is.na(info_df$novel_junc_left_coord))] <- counts_list_right[which(is.na(info_df$novel_junc_left_coord))]

  # Merge right and left PSI dataframes
  # Average PSI for right and left junctions
  # Add column identifying ID
  counts_list_merged <- map2(counts_list_right, counts_list_left, function(x, y) full_join(x, y, by = "sampleID", suffix = c("_r", "_l"))) %>%
    map(function(x) x %>% mutate(avgPSI = ifelse(!is.na(psi_r), ifelse(!is.na(psi_l), (psi_r + psi_l) / 2, 0), 0))) %>%
    modify2(info_df$novel_junc_id, function(x, y) x %>% mutate(novel_junc_id = y))

  # Return list of counts dataframes
  return(counts_list_merged)
}
