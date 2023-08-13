# CMSPlus
An R package for molecular subtyping by integrating the inter-tumor and intra-tumor heterogeneity of colorectal cancer

```
devtools::install_github("yswutan/CMSPlus")
library(CMSPlus)

## exp2symbol: a dataframe with Gene Expression Profiles data values,
##       samples in columns, genes in rows, rownames corresponding to gene symbols


CMSPlusLabels <- CMSPlus(exp2symbol, plot=TRUE, parallel.sz=0)

```