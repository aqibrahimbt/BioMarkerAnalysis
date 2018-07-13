context("test-de_analysis.R")

test_that("get valid comparisions out of conditions", {
  dataset = "GSE15229"
  dataset_data = readRDS(file = "conditions.rds")
  #saveRDS(dataset_data, file = "./gToolbox/tests/testthat/conditions.rds")
  conditions = dataset_data$conditions
  sample_group = dataset_data$sample_group
  comparisons = get_comparisons(conditions, sample_group)
  comparisons_ref =  readRDS(file = "comparisons.rds")
  #saveRDS(comparisons, file = "./gToolbox/tests/testthat/comparisons.rds")
  expect_identical(comparisons$conditions,comparisons_ref$conditions)
  expect_identical(comparisons$samples,comparisons_ref$samples)
})
