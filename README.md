## Manuscript Information
**Manuscript:** Large-scale RNA-seq mining reveals ciclopirox triggers TDP-43 cryptic exons 

**Available on *Nature Communications*:** [doi.org/10.1038/s41467-025-62004-5](https://doi.org/10.1038/s41467-025-62004-5)

**Contacts:** Irika Sinha at isinha1 [@] jh [.] edu & Dr. Jonathan Ling at jling [@] jhu [.] edu

## Description
This repository includes scripts used to analyze data for the manuscript listed above. It also includes the scripts use to create the original SnapMine interface on R Shiny. Snapmine is now hosted at [https://snapmine.idies.jhu.edu/](https://snapmine.idies.jhu.edu/). 

*Update (12/3/2025)*: This repository also includes an R package that can parse a csv file and identify samples with (pairs of) novel junctions in Snaptron. This functionality is mainly for querying many junctions at once. This package has NOT been robustly tested. Install using:

```{r class.source="bg-info"}
remotes::install_github("irikas/SnapMine_Paper", subdir = "snapmine")
```

## Abstract
Nuclear clearance and cytoplasmic aggregation of TDP-43 in neurons, initially identified in ALS-FTD, are hallmark pathological features observed across a spectrum of neurodegenerative diseases. We previously found that TDP-43 loss-of-function leads to the transcriptome-wide inclusion of deleterious cryptic exons in brains and biofluids post-mortem as well as during the presymptomatic stage of ALS-FTD, but upstream mechanisms that lead to TDP-43 dysregulation remain unclear. Here, we developed a web-based resource (SnapMine) to determine the levels of TDP-43 cryptic exon inclusion across hundreds of thousands of publicly available RNA sequencing datasets. We established cryptic exon inclusion across a variety of human cells and tissues to provide ground truth references for future studies on TDP-43 dysregulation. We then explored studies that were entirely unrelated to TDP-43 or neurodegeneration and found that ciclopirox olamine (CPX), an FDA-approved antifungal, can trigger the inclusion of TDP-43-associated cryptic exons in a variety of mouse and human primary cells. CPX induction of cryptic exon occurs via heavy metal toxicity and oxidative stress, suggesting that similar vulnerabilities could play a role in neurodegeneration. Our work demonstrates how diverse datasets can be linked through common biological features and underscores that public archives of sequencing data represent a vastly underutilized resource with tremendous potential for uncovering novel insights into complex biological mechanisms and diseases.
