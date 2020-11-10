context("SingleCellClassR class functions")
library(SingleCellClassR)

test_that("Set cell type changes cell type", {
  data("default_models")
  clf_B <- default_models[['B cells']]
  
  cell_type(clf_B) <- "b cells"
  expect_equal(clf_B@cell_type, "b cells")
})

test_that("Set probability threshold changes probability threshold", {
  data("default_models")
  clf_B <- default_models[['B cells']]
  
  p_thres(clf_B) <- 0.6
  expect_equal(clf_B@p_thres, 0.6)
})

test_that("Set classifier changes classifier and features", {
  data("default_models")
  clf_B <- default_models[['B cells']]
  clf_T <- default_models[['T cells']]
  
  clf(clf_B) <- clf(clf_T)
  expect_true(all.equal(clf_B@clf, clf(clf_T)))
  expect_true(all.equal(clf_B@features, clf_T@features))
})
