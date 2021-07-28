context("scAnnotatR support functions")
library(scAnnotatR)
library(Seurat)

test_that("Balance dataset balances dataset", {
  data("tirosh_mel80_example")
  example_data <- subset(tirosh_mel80_example, cells = 1:10)
  example_data <- t(as.matrix(GetAssayData(example_data)))
  example_tag <- c(rep('yes', 4), rep('no', 6))
  names(example_tag) <- rownames(example_data)
  
  balanced_ds <- balance_dataset(example_data, example_tag)
  expect_true(length(balanced_ds$tag[balanced_ds$tag == 'yes']) == 
               length(balanced_ds$tag[balanced_ds$tag == 'no']))
})

test_that("Z-score transformation transforms expression to z-score", {
  data("tirosh_mel80_example")
  example_data <- t(as.matrix(GetAssayData(tirosh_mel80_example)))
  
  z_data <- transform_to_zscore(example_data)
  
  sd_col <- unlist(apply(z_data, 2, sd))
  sd_compare <- unlist(lapply(sd_col, function(x) all.equal(x, 1)))
  # all TRUE returns TRUE
  expect_true(all(sd_compare))
  
  mean_col <- unlist(apply(z_data, 2, mean))
  # compare each mean col with 0 -> return list of TRUE/FALSE
  mean_compare <- unlist(lapply(mean_col, function(x) all.equal(x, 0)))
  # all TRUE returns TRUE
  expect_true(all(mean_compare))
})

test_that("marker genes selection selectes correct marker genes", {
  data("tirosh_mel80_example")
  example_data <- GetAssayData(tirosh_mel80_example)
  
  marker_genes <- c('CD19', 'MS4A1', 'CD8A', 'CD8B')
  filtered_mat <- select_marker_genes(example_data, marker_genes)
  expect_equal(sort(rownames(filtered_mat)), sort(marker_genes))
})

test_that("Check parent-child coherence checks parent-child coherence", {
  data("tirosh_mel80_example")
  example_data <- subset(tirosh_mel80_example, cells = 1:10)
  
  parent_pos <- colnames(example_data)[c(1:4, 7:10)]
  parent_cell <- 'ABC'
  
  tag <- c(rep(TRUE, 6), rep(FALSE, 4))
  names(tag) <- colnames(example_data)
  
  mat <- GetAssayData(example_data)
  new_tag <- check_parent_child_coherence(
    mat, tag, parent_pos, 'ABC', 'CDE', 'CDE')
  
  yes <- new_tag[new_tag == 'yes']
  no <- new_tag[new_tag == 'no']
  n_applicable <- new_tag[new_tag == 'not applicable']
  
  expect_equal(length(yes), 4)
  expect_equal(length(no), 4)
  expect_equal(length(n_applicable), 2)
})

test_that("Filter cells filters cells", {
  data("tirosh_mel80_example")
  example_data <- subset(tirosh_mel80_example, cells = 1:10)
  
  tag <- c(rep('ABC', 4), rep('not applicable', 4), rep(' +', 2))
  names(tag) <- colnames(example_data)
  mat <- GetAssayData(example_data)
  
  filtered_data <- filter_cells(mat, tag)
  expect_equal(ncol(filtered_data$mat), 4)
})

test_that("A tag vector construction works", {
  data("tirosh_mel80_example")
  example_data <- subset(tirosh_mel80_example, cells = 1:10)
  
  tag <- c(rep('ABC', 4), rep('CDE', 4), rep('DEF', 2))
  names(tag) <- colnames(example_data)
  
  bin_vect <- c(rep('yes', 4), rep('no', 6))
  names(bin_vect) <- colnames(example_data)
  
  tag_vect <- construct_tag_vect(tag, 'ABC')
  expect_equal(tag_vect, bin_vect)
})

test_that("Process parent classifier works", {
  data("tirosh_mel80_example")
  example_data <- subset(tirosh_mel80_example, cells = 1:10)
  mat <- GetAssayData(example_data)
  
  tag <- c(rep('ABC', 4), rep('CDE', 4), rep('DEF', 2))
  names(tag) <- colnames(example_data)
  
  return_val <- process_parent_classifier(mat, 
                                          parent_tag = tag, 
                                          'ABC', NULL, '.')
  
  expect_equal(sort(return_val$pos_parent), sort(colnames(example_data)[1:4]))
})
