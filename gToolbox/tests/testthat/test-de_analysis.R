context("test-de_analysis.R")

test_that("get valid comparisions out of conditions", {
  conditions = readRDS(file = "conditions.rds")
  #saveRDS(dataset_data$conditions, file = "./gToolbox/tests/testthat/conditions.rds")
  comparisons = get_comparisons(conditions)
  comparisons_ref =  readRDS(file = "comparisons.rds")
  #saveRDS(comparisons, file = "./gToolbox/tests/testthat/comparisons.rds")
  expect_identical(comparisons$conditions,comparisons_ref$conditions)
  expect_identical(comparisons$samples,comparisons_ref$samples)
})
