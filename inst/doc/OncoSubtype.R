## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----wilkerson, eval=FALSE----------------------------------------------------
#  library(OncoSubtype)
#  library(tidyverse)
#  set.seed(2121)
#  data <- get_median_centered(example_fpkm)
#  data <- assays(data)$centered
#  rownames(data) <- rowData(example_fpkm)$external_gene_name
#  # use default wilkerson's nearest centroids method
#  output1 <- centroids_subtype(data, disease = 'LUSC')
#  table(output1$subtypes)

## ----rf, eval=FALSE-----------------------------------------------------------
#  output2 <- ml_subtype(data, disease = 'LUSC', method = 'rf')
#  table(output1$subtypes)
#  confusionMatrix(output1, output2)

## ----confusion, eval=FALSE----------------------------------------------------
#  confusionMatrix(output1, output2)

## ----heatmap, eval=FALSE------------------------------------------------------
#  PlotHeat(object = output2, set = 'both', fontsize = 10,
#            show_rownames = FALSE, show_colnames = FALSE)

