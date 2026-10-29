# Third-party software

This app is Apache-2.0 (see `LICENSE`). The image bundles:

| Component | Version | Licence | Source | Role |
|---|---|---|---|---|
| CytoNorm | 2.0.12, commit 278f352d8edffd7a5d3e6b0b76c402d991a1e7a9 | GPL (>= 2) (DESCRIPTION; the repository has no LICENSE file) | https://github.com/saeyslab/CytoNorm | Training (`CytoNorm.train`, `QuantileNorm.train`), spline plots |
| FlowSOM | 2.20.0 | GPL (>= 2) | Bioconductor | Clustering, mapping new cells (`NewData`) |
| flowCore | 2.24.0 | Artistic-2.0 | Bioconductor | FCS reading |
| cytolib | 2.24.0 | AGPL-3.0-only | Bioconductor | flowCore's C++ FCS backend |
| ggplot2 | 4.0.3 | MIT | CRAN | Plots |
| jsonlite | 2.0.0 | MIT | CRAN | provenance.json |
| R | 4.6.1 | GPL-2 / GPL-3 | rocker/r-ver:4.6.1 | Runtime |
| Ubuntu 24.04 base | | various free licences | rocker/r-ver | OS |

## Licence of the image as a whole

Our R code is Apache-2.0. It loads the packages above into one R process, and several
are GPL (R itself, CytoNorm, FlowSOM, igraph, ggpubr and others in the table) and one is AGPL-3.0-only: **cytolib**, the C++ library behind
flowCore's FCS reader. The image as distributed is therefore a combined work under the
GNU GPL v3 / AGPL v3 terms (Apache-2.0 is compatible with GPLv3 in that direction). The
AGPL's network clause concerns software offered to users over a network; this app is a
batch container without a network interface. Source for every component is available
from the repositories listed (CRAN/Bioconductor sources at the versions shown; GitHub
packages at the pinned commits). R packages are redistributed under their own licences.

## Every R package in the image

R 4.6.1 (rocker/r-ver:4.6.1, Ubuntu 24.04), packages from Posit Public Package Manager
snapshot 2026-10-26 (CRAN) and Bioconductor 3.23. Generated from
`installed.packages()` in the built image (base packages omitted).

| Package | Version | Licence |
|---|---|---|
| abind | 1.4-8 | MIT + file LICENSE |
| ALL | 1.54.0 | Artistic-2.0 |
| backports | 1.5.1 | GPL-2 / GPL-3 |
| base64enc | 0.1-6 | GPL-2 / GPL-3 |
| BH | 1.90.0-1 | BSL-1.0 |
| Biobase | 2.72.0 | Artistic-2.0 |
| BiocGenerics | 0.58.1 | Artistic-2.0 |
| biocmake | 1.4.0 | MIT + file LICENSE |
| BiocManager | 1.30.27 | Artistic-2.0 |
| BiocVersion | 3.23.1 | Artistic-2.0 |
| boot | 1.3-32 | Unlimited |
| broom | 1.0.13 | MIT + file LICENSE |
| car | 3.1-5 | GPL (>= 2) |
| carData | 3.0-6 | GPL (>= 2) |
| class | 7.3-23 | GPL-2 / GPL-3 |
| cli | 3.6.6 | MIT + file LICENSE |
| cluster | 2.1.8.2 | GPL (>= 2) |
| codetools | 0.2-20 | GPL |
| colorRamps | 2.3.4 | GPL |
| colorspace | 2.1-4 | BSD_3_clause + file LICENSE |
| ConsensusClusterPlus | 1.76.0 | GPL version 2 |
| corrplot | 0.95 | MIT + file LICENSE |
| cowplot | 1.2.0 | GPL-2 |
| cpp11 | 0.5.5 | MIT + file LICENSE |
| cytolib | 2.24.0 | AGPL-3.0-only |
| CytoNorm | 2.0.12 | GPL (>= 2) |
| Deriv | 4.3.6 | GPL (>= 3) |
| dir.expiry | 1.20.0 | GPL-3 |
| doBy | 4.7.2 | GPL (>= 2) |
| docopt | 0.7.2 | MIT + file LICENSE |
| dplyr | 1.2.1 | MIT + file LICENSE |
| emdist | 0.3-3 | MIT + file LICENSE |
| farver | 2.1.2 | MIT + file LICENSE |
| filelock | 1.0.3 | MIT + file LICENSE |
| flowCore | 2.24.0 | Artistic-2.0 |
| FlowSOM | 2.20.0 | GPL (>= 2) |
| forecast | 9.0.2 | GPL-3 |
| foreign | 0.8-91 | GPL (>= 2) |
| Formula | 1.2-6 | GPL-2 / GPL-3 |
| fracdiff | 1.5-4 | GPL (>= 2) |
| generics | 0.1.4 | MIT + file LICENSE |
| ggforce | 0.5.0 | MIT + file LICENSE |
| ggnewscale | 0.5.2 | GPL-3 |
| ggplot2 | 4.0.3 | MIT + file LICENSE |
| ggpubr | 1.0.0 | GPL (>= 2) |
| ggrepel | 0.9.8 | GPL-3 / file LICENSE |
| ggridges | 0.5.7 | GPL-2 / file LICENSE |
| ggsci | 5.2.0 | GPL (>= 3) |
| ggsignif | 0.6.4 | GPL-3 / file LICENSE |
| glue | 1.8.1 | MIT + file LICENSE |
| gridExtra | 2.3.1 | GPL (>= 2) |
| gtable | 0.3.6 | MIT + file LICENSE |
| igraph | 2.3.4 | GPL (>= 2) |
| isoband | 0.3.0 | MIT + file LICENSE |
| jsonlite | 2.0.0 | MIT + file LICENSE |
| KernSmooth | 2.23-26 | Unlimited |
| labeling | 0.4.3 | MIT + file LICENSE / Unlimited |
| lattice | 0.22-9 | GPL (>= 2) |
| lifecycle | 1.0.5 | MIT + file LICENSE |
| littler | 0.3.23 | GPL (>= 2) |
| lme4 | 2.1-0 | GPL (>= 2) |
| lmtest | 0.9-40 | GPL-2 / GPL-3 |
| magrittr | 2.0.5 | MIT + file LICENSE |
| MASS | 7.3-65 | GPL-2 / GPL-3 |
| Matrix | 1.7-5 | GPL (>= 2) / file LICENCE |
| MatrixModels | 0.5-4 | GPL (>= 2) |
| matrixStats | 1.5.0 | Artistic-2.0 |
| mgcv | 1.9-4 | GPL (>= 2) |
| minqa | 1.2.8 | GPL-2 |
| modelr | 0.1.11 | GPL-3 |
| nlme | 3.1-169 | GPL (>= 2) |
| nloptr | 2.2.1 | LGPL (>= 3) |
| nnet | 7.3-20 | GPL-2 / GPL-3 |
| numDeriv | 2016.8-1.1 | GPL-2 |
| pbkrtest | 0.5.5 | GPL (>= 2) |
| pheatmap | 1.0.13 | GPL-2 |
| pillar | 1.11.1 | MIT + file LICENSE |
| pkgconfig | 2.0.3 | MIT + file LICENSE |
| polyclip | 1.10-7 | BSL |
| polynom | 1.4-1 | GPL-2 |
| purrr | 1.2.2 | MIT + file LICENSE |
| quantreg | 6.1 | GPL (>= 2) |
| R6 | 2.6.1 | MIT + file LICENSE |
| rbibutils | 2.4.1 | GPL-2 |
| RColorBrewer | 1.1-3 | Apache License 2.0 |
| Rcpp | 1.1.2 | GPL (>= 2) |
| RcppArmadillo | 15.6.2-1 | GPL (>= 2) |
| RcppEigen | 0.3.4.0.2 | GPL (>= 2) / file LICENSE |
| Rdpack | 2.6.6 | GPL (>= 2) |
| reformulas | 0.4.4 | GPL-3 |
| remotes | 2.5.0 | MIT + file LICENSE |
| Rhdf5lib | 2.0.0 | Artistic-2.0 |
| rlang | 1.3.0 | MIT + file LICENSE |
| rpart | 4.1.27 | GPL-2 / GPL-3 |
| RProtoBufLib | 2.24.0 | BSD_3_clause |
| rstatix | 1.1.0 | GPL-2 |
| Rtsne | 0.17 | file LICENSE |
| S4Vectors | 0.50.3 | Artistic-2.0 |
| S7 | 0.2.2 | MIT + file LICENSE |
| scales | 1.4.0 | MIT + file LICENSE |
| SparseM | 1.84-2 | GPL (>= 2) |
| spatial | 7.3-18 | GPL-2 / GPL-3 |
| stringi | 1.8.9 | file LICENSE |
| stringr | 1.6.0 | MIT + file LICENSE |
| survival | 3.8-6 | LGPL (>= 2) |
| systemfonts | 1.3.2 | MIT + file LICENSE |
| tibble | 3.3.1 | MIT + file LICENSE |
| tidyr | 1.3.2 | MIT + file LICENSE |
| tidyselect | 1.2.1 | MIT + file LICENSE |
| timeDate | 4052.112 | GPL (>= 2) |
| tweenr | 2.0.3 | MIT + file LICENSE |
| urca | 1.3-4 | GPL (>= 2) |
| utf8 | 1.2.6 | Apache License (== 2.0) / file LICENSE |
| vctrs | 0.7.3 | MIT + file LICENSE |
| viridisLite | 0.4.3 | MIT + file LICENSE |
| withr | 3.0.3 | MIT + file LICENSE |
| XML | 3.99-0.25 | BSD_3_clause + file LICENSE |
| zoo | 1.9-1 | GPL-2 / GPL-3 |
