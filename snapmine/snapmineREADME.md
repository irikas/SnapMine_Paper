
## snapmine Package Summary

`snapmine` enables querying junction inclusion of cryptic exons or other
junctions of interest in bulk. Its initial use is described in the
manuscript **Large-scale RNA-seq mining reveals ciclopirox triggers
TDP-43 cryptic exons** which was published in [Nature
Communications](https://doi.org/10.1038/s41467-025-62004-5).

## Installation

Please install this package using `remotes`.

``` r
install.packages("remotes")
library(remotes)

remotes::install_github("irikas/SnapMine_Paper", subdir = "snapmine")
library(snapmine)
```

## Usage

An introductory vignette on how to use this package can be found within
the “vignettes” directory. Using `GitHub & BitBucket HTML Preview`, you
can also view it
[here](https://htmlpreview.github.io/?https://github.com/irikas/SnapMine_Paper/blob/390846dad271d9ff50b3fc42bc848adbf88dbb48/snapmine/vignettes/introduction.html).

For querying only a few junctions, it may be easier to use the [online
interface](https://snapmine.idies.jhu.edu/).

## Algorithm

Below is Figure 1 from the publication. 1A is a schematic of the
SnapMine algorithm. ![SnapMine
schematic](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41467-025-62004-5/MediaObjects/41467_2025_62004_Fig1_HTML.png)
