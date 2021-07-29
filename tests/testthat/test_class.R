context("scAnnotatR class functions")
library(scAnnotatR)

test_that("Set cell type changes cell type", {
  default_models <- load_models('default')
  classifier_B <- default_models[['B cells']]
  
  cell_type(classifier_B) <- "b cells"
  expect_equal(classifier_B@cell_type, "b cells")
})

test_that("Set probability threshold changes probability threshold", {
  default_models <- load_models('default')
  classifier_B <- default_models[['B cells']]
  
  p_thres(classifier_B) <- 0.6
  expect_equal(classifier_B@p_thres, 0.6)
})

test_that("Set classifier changes classifier and marker genes", {
  default_models <- load_models('default')
  classifier_B <- default_models[['B cells']]
  classifier_T <- default_models[['T cells']]
  
  caret_model(classifier_B) <- caret_model(classifier_T)
  expect_true(all.equal(classifier_B@caret_model, classifier_T@caret_model))
  expect_true(all.equal(classifier_B@marker_genes, classifier_T@marker_genes))
})
