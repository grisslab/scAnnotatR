context("scTypeR support functions")
library(scTypeR)
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

test_that("Features selection selectes correct features", {
  data("tirosh_mel80_example")
  example_data <- GetAssayData(tirosh_mel80_example)
  
  features <- c('CD19', 'MS4A1', 'CD8A', 'CD8B')
  filtered_mat <- select_features(example_data, features)
  expect_equal(sort(rownames(filtered_mat)), sort(features))
})

test_that("Check parent-child coherence checks parent-child coherence", {
  data("tirosh_mel80_example")
  example_data <- subset(tirosh_mel80_example, cells = 1:10)
  
  parent_pos <- colnames(example_data)[c(1:4, 7:10)]
  parent_cell <- 'ABC'
  
  tag <- c(rep(TRUE, 6), rep(FALSE, 4))
  example_data[['CDE']] <- tag
  
  pc_coherence <- check_parent_child_coherence(
    example_data, parent_pos, 'ABC', 'CDE', 'CDE', tag_slot = 'CDE')
  
  new_obj <- pc_coherence$adjusted_object
  new_tag <- new_obj[['new_tag_slot']][, 1]
  
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
  example_data[['example_tag']] <- tag
  
  filtered_data <- filter_cells(example_data, 'example_tag')
  expect_equal(ncol(filtered_data), 4)
})

test_that("A tag vector construction works", {
  data("tirosh_mel80_example")
  example_data <- subset(tirosh_mel80_example, cells = 1:10)
  
  tag <- c(rep('ABC', 4), rep('CDE', 4), rep('DEF', 2))
  example_data[['example_tag']] <- tag
  
  bin_vect <- c(rep('yes', 4), rep('no', 6))
  names(bin_vect) <- colnames(example_data)
  
  tag_vect <- construct_tag_vect(example_data, 'ABC', 'example_tag')
  expect_equal(tag_vect, bin_vect)
})

test_that("Process parent clf works", {
  data("tirosh_mel80_example")
  example_data <- subset(tirosh_mel80_example, cells = 1:10)
  
  tag <- c(rep('ABC', 4), rep('CDE', 4), rep('DEF', 2))
  example_data[['example_tag']] <- tag
  
  return_val <- process_parent_clf(example_data, parent_tag_slot = 'example_tag', 
                                   'ABC', NULL, '.', seurat_assay = 'RNA', 
                                   seurat_slot = 'data')
  expect_equal(sort(return_val$pos_parent), sort(colnames(example_data)[1:4]))
})
