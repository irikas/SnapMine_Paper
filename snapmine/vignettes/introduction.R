## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, collapse = TRUE, comment = "#>")

## ----packages-----------------------------------------------------------------
remotes::install_github("irikas/SnapMine_Paper", subdir = "snapmine")
library(snapmine)
library(data.table)
library(tidyverse)
library(biomaRt)
library(knitr)

## ----files, eval=T------------------------------------------------------------
## Load df
df <- fread("/Users/irika/Library/CloudStorage/OneDrive-JohnsHopkins/WongLing/Scripts/250325_NMDi/data/251126_CEdf.csv")

## Original df did not have strand info so need to use biomaRt
## Update df
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
info <- getBM(attributes = c("ensembl_gene_id", "strand"), values = df$`Ensembl ID`, mart = mart)
df <- df %>%
  left_join(info, by = join_by(`Ensembl ID` == ensembl_gene_id)) %>%
  mutate(strand = ifelse(strand == 1, "+", "-"))
rm(mart, info)

## Create info_df file
info_df <- data.frame(
  novel_junc_id = paste(df$`Gene Symbol`, df$`ID #`, sep = "_"),
  compilation = "sra_human",
  strand = df$strand,
  novel_junc_left_coord = df$LeftJunction,
  novel_junc_right_coord = df$RightJunction,
  canon_junc_coord = df$`Canonical Splice Junction`
)

## ----query, eval = F----------------------------------------------------------
# ## Run query
# # For the sake of the example, run twice for long and wide dfs.
# # Could also run without flattening - this would lead to a list of dfs and then both flatten functions could be done manually using [flatten_counts]
# result_df_long <- calcPSI_bulk(info_df, flatten_long = T)
# result_df_wide <- calcPSI_bulk(info_df, flatten_wide = T)
# # saveRDS(result_df_long, "/Users/irika/Library/CloudStorage/OneDrive-JohnsHopkins/WongLing/Scripts/250325_NMDi/output_tables/251204_snapmine_results.RDS")
# # saveRDS(result_df_wide, "/Users/irika/Library/CloudStorage/OneDrive-JohnsHopkins/WongLing/Scripts/250325_NMDi/output_tables/251204_snapmine_result_df_wide.RDS")

## ----load_df, include = F-----------------------------------------------------
result_df_long <- readRDS("/Users/irika/Library/CloudStorage/OneDrive-JohnsHopkins/WongLing/Scripts/250325_NMDi/output_tables/251204_snapmine_results.RDS")
result_df_wide <- readRDS("/Users/irika/Library/CloudStorage/OneDrive-JohnsHopkins/WongLing/Scripts/250325_NMDi/output_tables/251204_snapmine_result_df_wide.RDS")

## ----info_df------------------------------------------------------------------
head(info_df)
knitr::kable(head(info_df))

## ----result_df_long-----------------------------------------------------------
head(result_df_long)

## ----result_df_wide-----------------------------------------------------------
head(result_df_wide)

## ----results_explore_all, fig.align = "c", out.width="100%", fig.dpi=300------
## How many different unique samples?
print(paste0("There are ", length(result_df_wide$sampleID), " unique samples in the results dataset"))

ggplot(result_df_long, aes(x = avgPSI)) +
  geom_histogram(binwidth = 1, color = "white", linewidth = 0.1) +
  theme_classic() +
  ggtitle(label = "Novel junctions present at 0% PSI in most samples") +
  theme(legend.position = "none", text = element_text(size = 5))

## Most included CE
result_df_wide_subset <- result_df_wide %>%
  mutate(maxPSI = apply(result_df_wide[, -1],
    MARGIN = 1,
    function(x) max(x, na.rm = T)
  ))

ggplot(result_df_wide_subset, aes(x = maxPSI)) +
  geom_histogram(binwidth = 1, color = "white", linewidth = 0.1) +
  theme_classic() +
  ggtitle(label = "Max PSI for most samples is 0%") +
  theme(legend.position = "none", text = element_text(size = 5))

## Subset to avgPSI > 15 for at least one CE
result_df_wide_subset <- result_df_wide_subset %>%
  filter(maxPSI > 15)

ggplot(result_df_wide_subset, aes(x = maxPSI)) +
  geom_histogram(binwidth = 1, color = "white", linewidth = 0.1) +
  theme_classic() +
  ggtitle(label = "Max PSI for most samples with PSI >15% ranges") +
  theme(legend.position = "none", text = element_text(size = 5))

result_df_long_subset <- result_df_long %>% filter(sampleID %in% result_df_wide_subset$sampleID)
print(paste0("There are ", length(result_df_wide_subset$sampleID), " unique samples in the results dataset with a PSI of at least 15% in one cryptic exon"))

# Subset to avgPSI > 15 for 10% or more of targets
result_df_wide_TDPpos <- result_df_wide_subset %>%
  dplyr::select(-maxPSI) %>%
  mutate(n15 = apply(result_df_wide_subset[, -1],
    MARGIN = 1,
    function(x) length(which(x > 15))
  )) %>%
  filter(n15 > round(0.10 * nrow(info_df), 0))

print(paste0(
  "There are ",
  length(result_df_wide_TDPpos$sampleID),
  " unique samples in the results dataset with a PSI of at least 10% in at least ",
  round(0.1 * nrow(info_df)), " cryptic exons"
))

result_df_long_TDPpos <- result_df_long_subset %>%
  filter(sampleID %in% result_df_wide_TDPpos$sampleID) %>%
  mutate(avgPSI = ifelse(is.na(avgPSI), 0, avgPSI))

ggplot(result_df_long_TDPpos, aes(x = sampleID, y = avgPSI)) +
  geom_jitter(width = 0.1, size = 0.1, alpha = 0.7) +
  geom_boxplot(alpha = 0.5, outliers = F, linewidth = 0.1) +
  theme_classic() +
  ggtitle(
    label = "Samples that are putatively lacking nuclear TDP-43",
    subtitle = "Contains CE in at least 10% of tested genes"
  ) +
  theme(
    legend.position = "none", text = element_text(size = 5),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

