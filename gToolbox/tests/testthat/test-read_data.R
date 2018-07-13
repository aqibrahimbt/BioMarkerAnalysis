context("test-read_data.R")

test_that("create condition out of new.sample.group and add to metadata", {
  dataset = "GSE15229"
  meta_data = readRDS(file = "meta_data.rds")
  #saveRDS(sub_samples, file = "./gToolbox/tests/testthat/meta_data.rds")
  
  dataset_metadata_ref = readRDS(file = "dataset_meta_data.rds")
  dataset_metadata = get.metadata(meta_data)
  #saveRDS(dataset_metadata, file = "./gToolbox/tests/testthat/dataset_meta_data.rds")
  expect_identical(dataset_metadata,dataset_metadata_ref)
})

test_that("check covariates and conditions, update formula", {
  formula = as.formula("~age+gender+condition")
  formula_ref = as.formula("~condition")
  dataset_metadata = readRDS(file = "dataset_meta_data.rds")
  conditions = get.conditons(dataset_metadata, formula)
  conditions_ref = readRDS(file = "conditions.rds")
  #saveRDS(dataset_metadata, file = "./gToolbox/tests/testthat/conditions.rds")
  expect_identical(conditions$conditions,conditions_ref$conditions)
  expect_equal(conditions$formula,formula_ref)
})

test_that("create a formula out of a vector", {
  names = c("hans","wurst","peter")
  formula_ref = as.formula("~hans + wurst + peter")
  formula = create.new.formula(names)
  expect_equal(formula,formula_ref)
})