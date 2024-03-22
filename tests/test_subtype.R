library(testthat)
library(OncoSubtype)

test_that(desc = "Whether subtype function works properly", {
  data <- get_median_centered(example_fpkm)
  data <- assays(data)$centered
  rownames(data) <- rowData(example_fpkm)$external_gene_name
  res <- ml_subtype(data)
  expect_equal(class(res)[1], "SubtypeClass")
})
