#' Train cell type classifier
#' 
#' @description Train a classifier for a new cell type. 
#' If cell type has a parent, only available for \code{\link{scClassifR}}
#' object as parent cell classifying model.
#' 
#' @param train_obj object that can be used for training the new model. 
#' \code{\link{Seurat}} object or \code{\link{SingleCellExperiment}} object
#' is expected.
#' If the training model has parent, parent_tag_slot may have been indicated. 
#' This field would have been filled out automatically 
#' if user precedently run classify_cells function. 
#' If no (predicted) cell type annotation provided, 
#' the function can be run if 1- parent_cell or 2- parent_clf is provided.
#' @param cell_type string indicating the name of the subtype
#' This must exactly match cell tag/label if cell tag/label is a string.
#' @param features list of features used for the new training model
#' @param parent_cell string indicated the name of the parent cell type, 
#' if parent cell type classifier has already been saved in model database.
#' Adjust path_to_models for exact database.  
#' @param parent_clf classification model for the parent cell type
#' @param path_to_models path to the folder containing the model database. 
#' As default, the pretrained models in the package will be used. 
#' If user has trained new models, indicate the folder containing the 
#' new_models.rda file.
#' @param zscore whether gene expression in train_obj is transformed to zscore
#' @param ... arguments passed to other methods
#' 
#' @return \code{\link{scClassifR}} object
#'
#' @note Only one cell type is expected for each cell in object. 
#' Ambiguous cell type, such as: "T cells/NK cells/ILC", 
#' will be ignored from training.
#' Subtypes used in training model for parent cell types must be indicated
#' as parent cell type. For example, when training for B cells, 
#' plasma cells must be annotated as B cells in order to be used.
#' 
#' @export
setGeneric("train_classifier", 
           function(train_obj, cell_type, features, 
                    parent_cell = NA_character_, 
                    parent_clf = NULL, path_to_models = c("default", "."), 
                    zscore = TRUE, ...) 
             standardGeneric("train_classifier"))

#' @inherit train_classifier
#' 
#' @param seurat_tag_slot string, name of slot in cell meta data 
#' indicating cell tag/label in the training object.
#' Strings indicating cell types are expected in this slot. 
#' For \code{\link{Seurat}} object, default value is "active.ident".  
#' Expected values are string or binary/logical, 
#' 0/"no"/F/FALSE: not being new cell type, 1/"yes"/T/TRUE: being new cell type.
#' @param seurat_parent_tag_slot string, name of a slot in cell meta data 
#' indicating assigned/predicted cell type. Default is "predicted_cell_type". 
#' This slot would have been filled automatically 
#' if user have called classify_cells function.
#' The slot must contain only string values. 
#' @param seurat_assay name of assay to use in training object. 
#' Default to 'RNA' assay.
#' @param seurat_slot type of expression data to use in training object. 
#' For \code{\link{Seurat}} object, available types are: "counts", "data" 
#' and "scale.data". Default to "counts", which contains unnormalized data.
#' 
#' @examples
#' # load small example dataset
#' data("tirosh_mel80_example")
#' 
#' # this dataset already contains pre-defined cell labels
#' table(Seurat::Idents(tirosh_mel80_example))
#' 
#' # define genes to use to classify this cell type (B cells in this example)
#' selected_features_B = c("CD19", "MS4A1", "CD79A")
#' 
#' # train the classifier, the "cell_type" argument must match 
#' # the cell labels in the data, except upper/lower case
#' set.seed(123)
#' clf_b <- train_classifier(train_obj = tirosh_mel80_example, 
#' features = selected_features_B, cell_type = "b cells")
#' 
#' # classify cell types using B cell classifier, 
#' # a test classifier process may be used before applying the classifier 
#' tirosh_mel80_example <- classify_cells(classify_obj = tirosh_mel80_example, 
#' classifiers = c(clf_b))
#' 
#' # tag all cells that are plasma cells (random example here)
#' tirosh_mel80_example[['plasma_cell_tag']] <- c(rep(1, 80), rep(0, 400))
#' 
#' # set new features for the subtype
#' p_features = c("SDC1", "CD19", "CD79A")
#' 
#' # train the classifier, the "B cell" classifier is used as parent. 
#' # This means, only cells already classified as "B cells" will be evaluated.
#' # the "tag_slot" parameter tells the classifier to use this cell meta data
#' # for the training process.
#' set.seed(123)
#' plasma_clf <- train_classifier(train_obj = tirosh_mel80_example, 
#' cell_type = "Plasma cell", features = p_features, parent_clf = clf_b, 
#' seurat_tag_slot = 'plasma_cell_tag')
#' 
#' @importFrom Seurat GetAssayData
#' 
#' @rdname train_classifier
setMethod("train_classifier", c("train_obj" = "Seurat"), 
          function(train_obj, cell_type, features, parent_cell = NA_character_,
                   parent_clf = NULL, path_to_models = c("default", "."), 
                   zscore = TRUE, seurat_tag_slot = "active.ident", 
                   seurat_parent_tag_slot = "predicted_cell_type", 
                   seurat_assay = 'RNA', seurat_slot = 'counts', ...) {
  #--- part of parent cell type
  parent_process <- process_parent_clf(
    train_obj, seurat_parent_tag_slot, parent_cell, parent_clf, 
    path_to_models, zscore, seurat_assay, seurat_slot
  )
  
  # check parent-child coherence
  if (!is.null(parent_process$pos_parent)) {
    check_res <- check_parent_child_coherence(
      train_obj, parent_process$pos_parent, parent_process$parent_cell, 
      cell_type, cell_type, seurat_tag_slot
    )
    train_obj <- check_res$adjusted_object
    seurat_tag_slot <- check_res$adjusted_tag_slot
  }
  
  #--- end of part of parent cell type
  
  #--- part of cell type
  # filter cells
  train_obj <- filter_cells(train_obj, seurat_tag_slot)
  
  # convert Seurat object to matrix
  mat = Seurat::GetAssayData(object = train_obj, 
                             assay = seurat_assay, slot = seurat_slot)
  
  # feature selection
  mat <- select_features(mat, features)
  
  # transpose mat
  mat <- t(as.matrix(mat))
  
  # transform mat to zscore values
  if (zscore == TRUE)
    mat <- transform_to_zscore(mat)
  
  # construct cell tag to yes/no values
  train_tag <- construct_tag_vect(train_obj, cell_type, seurat_tag_slot)
  
  # if cell type not found in tag
  if (all(train_tag != "yes")) {
    stop("Cell type not available in train data. Please verify cell type.", 
         call. = FALSE)
  }

  # transform list to factor
  train_tag <- factor(train_tag, levels = c('yes', 'no'))
  
  # convert hyphen (-) by underscore (_)
  colnames(mat) <- gsub('-', '_', colnames(mat))
  
  # train
  clf <- train_func(mat, train_tag)
  
  # remove this info to reduce memory
  clf$resampledCM <- NULL 
  p_thres <- 0.5
  
  features <- labels(clf$terms)
  features <- gsub('_', '-', features) # convert back underscore to hyphen
  object <- scClassifR(cell_type, clf, features, p_thres, 
                             NA_character_)
  
  # only assign parent if pretrained model for parent cell type is avai
  parent_check <- 
    (
      !is.null(parent_process$parent.clf) && 
        tolower(cell_type(parent_process$parent.clf)) == 
        tolower(parent_process$parent_cell)
    ) || (
      tolower(parent_process$parent_cell) %in% 
        tolower(names(parent_process$model_list))
    )
  if (parent_check) parent(object) <- parent_process$parent_cell
  
  return(object)
})

#' @inherit train_classifier
#' 
#' @param sce_tag_slot string, name of annotation slot indicating 
#' cell tag/label in the training object.
#' For \code{\link{SingleCellExperiment}} object, default value is "ident".  
#' Expected values are string or binary/logical, 
#' 0/"no"/F/FALSE: not being new cell type, 1/"yes"/T/TRUE: being new cell type.
#' @param sce_parent_tag_slot string, name of a slot in cell meta data 
#' indicating pre-assigned/predicted cell type. 
#' Default field is "predicted_cell_type".
#' This field would have been filled automatically 
#' when user called classify_cells function. 
#' The slot must contain only string values. 
#' @param sce_assay name of assay to use in training object. 
#' Default to 'logcounts' assay.
#' 
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment assay
#' 
#' @rdname train_classifier
setMethod("train_classifier", c("train_obj" = "SingleCellExperiment"), 
          function(train_obj, cell_type, features, parent_cell = NA_character_,
                   parent_clf = NULL, path_to_models = c("default", "."), 
                   zscore = TRUE, sce_tag_slot = "ident", 
                   sce_parent_tag_slot = "predicted_cell_type", 
                   sce_assay = 'logcounts', ...) {
  # solve duplication of cell names
  colnames(train_obj) <- make.unique(colnames(train_obj), sep = '_')
            
  #--- part of parent cell type
  parent_process <- process_parent_clf(
    train_obj, sce_parent_tag_slot, parent_cell, parent_clf, 
    path_to_models, zscore, sce_assay
  )
  
  # check parent-child coherence
  if (!is.null(parent_process$pos_parent)) {
    check_res <- check_parent_child_coherence(
      train_obj, parent_process$pos_parent, parent_process$parent_cell,
      cell_type,cell_type,sce_tag_slot
    )
    train_obj <- check_res$adjusted_object
    sce_tag_slot <- check_res$adjusted_tag_slot
  }
  
  #--- end of part of parent cell type
  
  #--- part of cell type
  # filter cells
  train_obj <- filter_cells(train_obj, sce_tag_slot)
  
  # convert Seurat object to matrix
  mat = SummarizedExperiment::assay(train_obj, sce_assay)
  
  # feature selection
  mat <- select_features(mat, features)
  
  # transpose mat
  mat <- t(as.matrix(mat))
  
  # transform mat to zscore values
  if (zscore == TRUE)
    mat <- transform_to_zscore(mat)
  
  # construct cell tag to yes/no values
  train_tag <- construct_tag_vect(train_obj, cell_type, sce_tag_slot)
  
  # if cell type not found in tag
  if (all(train_tag != "yes")) {
    stop("Cell type not available in train data. Please verify cell type.", 
         call. = FALSE)
  }
  
  # transform list to factor
  train_tag <- factor(train_tag, levels = c('yes', 'no'))
  
  # convert hyphen (-) by underscore (_)
  colnames(mat) <- gsub('-', '_', colnames(mat))
  
  # train
  clf <- train_func(mat, train_tag)
  
  # remove this info to reduce memory
  clf$resampledCM <- NULL 
  p_thres <- 0.5
  
  features <- labels(clf$terms)
  features <- gsub('_', '-', features) # convert back underscore to hyphen
  object <- scClassifR(cell_type, clf, features, p_thres, 
                             NA_character_)
  
  # only assign parent if pretrained model for parent cell type is avai
  parent_check <-
    (
      !is.null(parent_process$parent.clf) && 
        tolower(cell_type(parent_process$parent.clf)) == 
        tolower(parent_process$parent_cell)
    ) || (
      tolower(parent_process$parent_cell) %in% 
        tolower(names(parent_process$model_list))
    )
  if (parent_check) parent(object) <- parent_process$parent_cell
  
  return(object)
})

#' Testing process.
#' 
#' @description Testing process. 
#' 
#' @param test_obj xxobject that can be used for testing
#' @param classifier classification model
#' @param target_cell_type vector indicating other cell types than cell labels 
#' that can be considered as the main cell type in classifier, 
#' for example, c("plasma cell", "b cell", "b cells", "activating b cell"). 
#' Default as NULL.
#' @param parent_clf \code{\link{scClassifR}} object
#' corresponding to classification model for the parent cell type
#' @param path_to_models path to the folder containing the list of models. 
#' As default, the pretrained models in the package will be used. 
#' If user has trained new models, indicate the folder containing 
#' the new_models.rda file.
#' @param zscore boolean, whether gene expression is transformed to zscore
#' @param ... arguments passed to other methods
#'
#' @return result of testing process in form of a list, 
#' including predicted values, prediction accuracy at a probability threshold, 
#' and roc curve information.
#' 
#' @note Only one cell type is expected for each cell. 
#' Ambiguous cell type, such as: "T cells/NK cells/ILC", will be ignored.
#' Subtypes used in testing model for parent cell types can be indicated 
#' as parent cell type, or can be indicated in target_cell_type. 
#' For example, when testing for B cells, plasma cells can be annotated as 
#' B cells, or target_cell_type is set c("plasma cells").
#' 
#' @export
setGeneric("test_classifier", function(test_obj, classifier, 
                                       target_cell_type = NULL, 
                                       parent_clf = NULL, 
                                       path_to_models = c("default", "."), 
                                       zscore = TRUE, ...) 
  standardGeneric("test_classifier"))

#' @inherit test_classifier
#' 
#' @param seurat_tag_slot string, name of annotation slot 
#' indicating cell tag/label in the testing object.
#' Strings indicating cell types are expected in this slot. 
#' For \code{\link{Seurat}} object, default value is "active.ident". 
#' Expected values is string or binary/logical, 
#' 0/"no"/F/FALSE: not being new cell type, 
#' 1/"yes"/T/TRUE: being new cell type.
#' @param seurat_parent_tag_slot string, name of tag slot in cell meta data
#' indicating pre-assigned/predicted parent cell type. 
#' Default field is "predicted_cell_type".
#' The slot must contain only string values. 
#' @param seurat_assay name of assay to use in 
#' \code{\link{Seurat}} object, defaults to 'RNA' assay.
#' @param seurat_slot type of expression data to use in 
#' \code{\link{Seurat}} object. 
#' Some available types are: "counts", "data" and "scale.data". 
#' Default to "counts", which contains unnormalized data.
#'  
#' @examples
#' # load small example dataset
#' data("tirosh_mel80_example")
#' 
#' # train the classifier
#' selected_features_B = c("CD19", "MS4A1", "CD79A")
#' set.seed(123)
#' clf_b <- train_classifier(train_obj = tirosh_mel80_example, 
#' features = selected_features_B, cell_type = "B cells")
#' 
#' # test the classifier, target cell type can be in other formats or
#' # alternative cell type that can be considered as the classified cell type 
#' clf_b_test <- test_classifier(test_obj = tirosh_mel80_example, 
#' classifier = clf_b, target_cell_type = c("B cell"))
#' clf_b_test
#' 
#' @importFrom Seurat GetAssayData
#'
#' @rdname test_classifier
setMethod("test_classifier", c("test_obj" = "Seurat", 
                               "classifier" = "scClassifR"), 
          function(test_obj, classifier, target_cell_type = NULL, 
                   parent_clf = NULL, path_to_models = c("default", "."), 
                   zscore = TRUE, seurat_tag_slot = "active.ident", 
                   seurat_parent_tag_slot = "predicted_cell_type", 
                   seurat_assay = 'RNA', seurat_slot = 'counts', ...) {
  . <- fpr <- tpr <- NULL
  
  # target_cell_type check
  if (!tolower(cell_type(classifier)) %in% tolower(target_cell_type)) {
    target_cell_type <- append(target_cell_type, cell_type(classifier))
  }
  
  #--- parent cell type
  # process parent clf
  parent_process <- process_parent_clf(
    test_obj, seurat_parent_tag_slot, parent(classifier), parent_clf,
    path_to_models, zscore, seurat_assay, seurat_slot
  )
  
  # check parent-child coherence
  if (!is.null(parent_process$pos_parent)) {
    check_res <- check_parent_child_coherence(
      test_obj, parent_process$pos_parent, parent(classifier), 
      cell_type(classifier), target_cell_type, seurat_tag_slot
    )
    test_obj <- check_res$adjusted_object
    seurat_tag_slot <- check_res$adjusted_tag_slot
  }
  
  #--- children cell type 
  # filter cells
  test_obj <- filter_cells(test_obj, seurat_tag_slot)
  
  # convert Seurat object to matrix
  test_mat = Seurat::GetAssayData(
    object = test_obj, assay = seurat_assay, slot = seurat_slot)
  
  # perform features selection
  test_mat <- select_features(test_mat, features(classifier))
  
  # transpose mat
  test_mat <- t(as.matrix(test_mat))
  
  # transform mat to zscore values
  if (zscore == TRUE) test_mat <- transform_to_zscore(test_mat)
  
  # construct cell tag to yes/no values
  test_tag <- construct_tag_vect(test_obj, target_cell_type, seurat_tag_slot)
  
  # if cell type not found in tag
  if (all(test_tag != "yes")) {
    stop("Cell type ", cell_type(classifier), 
    " is not available in the test data. 
    Overwrite the target cell type using the target_cell_type parameter 
    or verify that you chose the correct test dataset.", 
         call. = FALSE)
  }
  
  return_val = test_performance(test_mat, classifier, test_tag)
  return(return_val)
})

#' @inherit test_classifier
#' 
#' @param sce_tag_slot string, name of annotation slot 
#' indicating cell tag/label in the testing object.
#' Strings indicating cell types are expected in this slot. 
#' Default value is "ident".  
#' Expected values are string or binary/logical, 
#' 0/"no"/F/FALSE: not being new cell type, 
#' 1/"yes"/T/TRUE: being new cell type.
#' @param sce_parent_tag_slot string, name of tag slot in cell meta data
#' indicating pre-assigned/predicted parent cell type. 
#' Default is "predicted_cell_type".
#' The slot must contain only string values. 
#' @param sce_assay name of assay to use in \code{\link{SingleCellExperiment}}
#' object, defaults to 'logcounts' assay.
#'  
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment assay
#' 
#' @rdname test_classifier
setMethod("test_classifier", c("test_obj" = "SingleCellExperiment", 
                               "classifier" = "scClassifR"), 
          function(test_obj, classifier, target_cell_type = NULL, 
                   parent_clf = NULL, path_to_models = c("default", "."), 
                   zscore = TRUE, sce_tag_slot = "ident", 
                   sce_parent_tag_slot = "predicted_cell_type", 
                   sce_assay = 'logcounts', ...) {
  # solve duplication of cell names
  colnames(test_obj) <- make.unique(colnames(test_obj), sep = '_')
  . <- fpr <- tpr <- NULL
  
  # target_cell_type check
  if (!tolower(cell_type(classifier)) %in% tolower(target_cell_type)) {
    target_cell_type <- append(target_cell_type, cell_type(classifier))
  }
  
  #--- parent cell type
  # process parent clf
  parent_process <- process_parent_clf(
    test_obj, sce_parent_tag_slot, parent(classifier), parent_clf, 
    path_to_models, zscore, sce_assay
  )
  
  # check parent-child coherence
  if (!is.null(parent_process$pos_parent)) {
    check_res <- check_parent_child_coherence(
      test_obj, parent_process$pos_parent, parent(classifier), 
      cell_type(classifier), target_cell_type, sce_tag_slot
    )
    test_obj <- check_res$adjusted_object
    sce_tag_slot <- check_res$adjusted_tag_slot
  }
  
  #--- children cell type 
  # filter cells
  test_obj <- filter_cells(test_obj, sce_tag_slot)
  
  # convert Seurat object to matrix
  test_mat = SummarizedExperiment::assay(test_obj, sce_assay)
  
  # perform features selection
  test_mat <- select_features(test_mat, features(classifier))
  
  # transpose mat
  test_mat <- t(as.matrix(test_mat))
  
  # transform mat to zscore values
  if (zscore == TRUE) test_mat <- transform_to_zscore(test_mat)
  
  # construct cell tag to yes/no values
  test_tag <- construct_tag_vect(test_obj, target_cell_type, sce_tag_slot)
  
  # if cell type not found in tag
  if (all(test_tag != "yes")) {
    stop("Cell type ", cell_type(classifier), " is not available in the test data. 
          Overwrite the target cell type using the target_cell_type parameter or verify that you chose the correct test dataset.", 
         call. = FALSE)
  }
  
  return_val <- test_performance(test_mat, classifier, test_tag)
  
  return(return_val)
})

#' Plot roc curve
#' 
#' @param test_result result of test_classifier function
#' 
#' @return ggplot2 roc curve
#' @examples
#' # load small example dataset
#' data("tirosh_mel80_example")
#' 
#' # train a classifier, for ex: B cell
#' selected_features_B = c("CD19", "MS4A1", "CD79A")
#' set.seed(123)
#' clf_b <- train_classifier(train_obj = tirosh_mel80_example, 
#' features = selected_features_B, cell_type = "B cells")
#' 
#' clf_b_test <- test_classifier(test_obj = tirosh_mel80_example, 
#' classifier = clf_b)
#' 
#' # run plot curve on the test result
#' roc_curve <- plot_roc_curve(test_result = clf_b_test)
#' @import ROCR
#' @import ggplot2
#' @export
plot_roc_curve <- function(test_result) {
  fpr <- tpr <- NULL
  
  data <- rbind(c(0, 0), test_result$overall_roc[, c(2:3)], c(1, 1))
  colnames(data) <- c('fpr', 'tpr')
  
  # plot ROC curve
  q <- ggplot2::ggplot(data=data.frame(data), aes(x=fpr, y=tpr)) 
  q <- q + ggplot2::geom_line() 
  q <- q + ggplot2::xlab("False Positive Rate (1-Specificity)") 
  q <- q + ggplot2::ylab("True Positive Rate (Sensitivity)") 
  q
  return(q)
}

#' Classify cells from multiple models
#' 
#' @param classify_obj the object containing cells to be classified
#' @param classifiers list of classification models. 
#' The model is obtained from train_classifier function or available in current 
#' working space. 
#' Users may test the model using test_classifier before using this function.
#' If classifiers contain classifiers for sub cell types, classifiers for
#' parent cell type must be indicated first in order to be applied before
#' children classifiers.
#' If classifiers is NULL, the method will use all classifiers in database.
#' @param cell_types list of cell types containing models to be used
#' for classification, only applicable if the models have been saved to package.
#' @param path_to_models path to the folder containing the list of models. 
#' As default value, the pretrained models in the package will be used. 
#' If user has trained new models, indicate the folder containing the 
#' new_models.rda file.
#' @param ignore_ambiguous_result return all ambiguous predictions
#' (multiple cell types) to empty
#' When this parameter turns to TRUE, 
#' most probably predicted cell types will be ignored.  
#' @param cluster_slot name of slot in meta data containing cluster 
#' information, in case users want to have additional cluster-level 
#' prediction
#' @param ... arguments passed to other methods
#' 
#' @return the input object with new slots in cells meta data
#' New slots are: predicted_cell_type, most_probable_cell_type,
#' slots in form of [cell_type]_p and [cell_type]_class 
#' 
#' @export
setGeneric("classify_cells", function(classify_obj, classifiers = NULL, 
                                      cell_types = "all", 
                                      path_to_models = c("default", "."),
                                      ignore_ambiguous_result = FALSE, 
                                      cluster_slot = NULL, ...) 
  standardGeneric("classify_cells"))

#' @inherit classify_cells
#' 
#' @param seurat_assay name of assay to use in 
#' \code{\link{Seurat}} object, defaults to 'RNA' assay.
#' @param seurat_slot type of expression data to use in 
#' \code{\link{Seurat}} object. Some available types are: 
#' "counts", "data" and "scale.data". 
#' Default to "counts", which is unnormalized data.
#' 
#' @examples
#' # load small example dataset
#' data("tirosh_mel80_example")
#' 
#' # train one classifier for one cell type, for ex, B cell
#' # define genes to use to classify this cell type
#' selected_features_B = c("CD19", "MS4A1", "CD79A")
#' 
#' # train the classifier
#' set.seed(123)
#' clf_b <- train_classifier(train_obj = tirosh_mel80_example, 
#' features = selected_features_B, cell_type = "B cells")
#' 
#' # do the same thing with other cell types, for example, T cells
#' selected_features_T = c("CD4", "CD8A", "CD8B")
#' set.seed(123)
#' clf_t <- train_classifier(train_obj = tirosh_mel80_example, 
#' features = selected_features_T, cell_type = "T cells")
#' 
#' # create a list of classifiers
#' classifier_ls <- list(clf_b, clf_t)
#' 
#' # classify cells with list of classifiers
#' seurat.obj <- classify_cells(classify_obj = tirosh_mel80_example, 
#' classifiers = classifier_ls)
#' 
#' @importFrom Seurat GetAssayData
#' @import dplyr
#' @importFrom stats predict
#' 
#' @rdname classify_cells
setMethod("classify_cells", c("classify_obj" = "Seurat"), 
          function(classify_obj, classifiers = NULL, cell_types = "all", 
                   path_to_models = c("default", "."), 
                   ignore_ambiguous_result = FALSE, cluster_slot = 'seurat_clusters', 
                   seurat_assay = 'RNA', seurat_slot = 'counts', ...) {
  if (is.null(classifiers)) { 
    model_list <- load_models(path_to_models)
    
    if (cell_types >= 1 && cell_types != "all") 
      classifiers = model_list[cell_types]
    else classifiers <- model_list
  }
  
  union.features <- unique(unname(unlist(lapply(classifiers, 
                                                function(x) features(x)))))
  mat = Seurat::GetAssayData(object = classify_obj, 
                             assay = seurat_assay, slot = seurat_slot)
  # reduce features to reduce computational complexity
  mat <- select_features(mat, union.features)
  mat <- t(transform_to_zscore(t(as.matrix(mat))))
  
  chunk.size = 5000
  nchunks = ceiling(ncol(classify_obj)/chunk.size)
  for (i in seq(1, nchunks)) {
    idx.chunk = seq((i - 1) * chunk.size + 1, 
                    min(ncol(classify_obj), i * chunk.size))
    obj.chunk <- classify_obj[, idx.chunk]
    mat.chunk <- mat[, idx.chunk]
    
    # create an empty cell type for all cells
    pred_cells <- c(rep("", ncol(mat.chunk))) 
    names(pred_cells) <- colnames(mat.chunk)
    
    # run predictors
    for (classifier in classifiers) {
      if (!is.na(parent(classifier))) {
        applicable_mat <- verify_parent(mat.chunk, classifier, obj.chunk[[]])
        if (is.null(applicable_mat)) next # no parent clf provided or no positive to parent clf
      } else applicable_mat <- mat.chunk
  
      filtered_mat <- select_features(applicable_mat, features(classifier))
      filtered_mat <- t(as.matrix(filtered_mat))
      
      prediction <- make_prediction(
        filtered_mat, classifier, pred_cells, ignore_ambiguous_result
      )
      pred <- prediction$pred
      pred_cells <- prediction$pred_cells
      
      # add prediction to meta data: each classifier/celltype has 2 col: p, class
      for (colname in colnames(pred)) {
        obj.chunk[[colname]] <- pred[, colname, drop = FALSE]
      }
    }
    
    if (any(pred_cells != "")) {
      pred_cells <- gsub("/$", "", pred_cells)
      
      # ignore ambiguous results
      if (ignore_ambiguous_result == TRUE) {
        pred_cells <- unlist(
          lapply(pred_cells, function(x) 
            if (length(unlist(strsplit(x, split = '/'))) <= 1) {x} 
            else {'ambiguous'})
        )
      } 
      
      # add cell type to meta data
      obj.chunk[['predicted_cell_type']] <- pred_cells
      
      # simplify result can only happen when not ignore ambiguous results
      if (ignore_ambiguous_result == FALSE) 
        obj.chunk[['most_probable_cell_type']] <- simplify_prediction(
          obj.chunk[[]], as.matrix(mat.chunk), classifiers)
    }
    
    if (i == 1) classified_obj <- obj.chunk
    else classified_obj <- merge(classified_obj, obj.chunk)
  }
  
  if (!is.null(cluster_slot) 
      & cluster_slot %in% colnames(classified_obj[[]]) 
      & !ignore_ambiguous_result) {
    clusts <- as.factor(classified_obj[[]][, cluster_slot])
    classified_obj$clust_pred <- 
      classify_clust(clusts, classified_obj$most_probable_cell_type)
  }
  
  new.cols <- colnames(classified_obj[[]])[!colnames(classified_obj[[]]) 
                                           %in% colnames(classify_obj[[]])]
  for (colname in new.cols)
    classify_obj[[colname]] <- classified_obj[[colname]]
  return(classify_obj)
})

#' @inherit classify_cells
#' 
#' @param sce_assay name of assay to use in 
#' \code{\link{SingleCellExperiment}} object, 
#' defaults to 'logcounts' assay.
#' 
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment assay colData
#' 
#' @rdname classify_cells
setMethod("classify_cells", c("classify_obj" = "SingleCellExperiment"), 
          function(classify_obj, classifiers = NULL, cell_types = "all", 
                   path_to_models = c("default", "."), 
                   ignore_ambiguous_result = FALSE, 
                   sce_assay = 'logcounts', cluster_slot = NULL, ...) {
  # solve duplication of cell names
  colnames(classify_obj) <- make.unique(colnames(classify_obj), sep = '_')
  
  if (is.null(classifiers)) { 
    model_list <- load_models(path_to_models)
    
    if (cell_types>=1 && cell_types!="all") classifiers = model_list[cell_types]
    else classifiers <- model_list
  }
  
  union.features <- unique(unname(unlist(lapply(classifiers, 
                                                function(x) features(x)))))
  mat = SummarizedExperiment::assay(classify_obj, sce_assay)
  # reduce features to reduce computational complexity
  mat <- select_features(mat, union.features)
  mat <- t(transform_to_zscore(t(as.matrix(mat))))
  
  # split dataset into multiple chunks to reduce running time
  chunk.size = 5000
  nchunks = ceiling(ncol(classify_obj)/chunk.size)
  for (i in seq(1, nchunks)) {
    idx.chunk = seq((i - 1) * chunk.size + 1, 
                    min(ncol(classify_obj), i * chunk.size))
    obj.chunk <- classify_obj[, idx.chunk]
    mat.chunk <- mat[, idx.chunk]
    
    # create an empty cell type for all cells
    pred_cells <- c(rep("", ncol(mat.chunk))) 
    names(pred_cells) <- colnames(mat.chunk)
    
    # run predictors
    for (classifier in classifiers) {
      if (!is.na(parent(classifier))) { 
        applicable_mat <- verify_parent(mat.chunk, classifier, colData(obj.chunk))
        if (is.null(applicable_mat)) next # no parent clf provided or no positive to parent clf
      } else applicable_mat <- mat.chunk
      
      filtered_mat <- select_features(applicable_mat, features(classifier))
      filtered_mat <- t(as.matrix(filtered_mat))
      
      prediction <- make_prediction(
        filtered_mat, classifier, pred_cells, ignore_ambiguous_result)
      
      pred <- prediction$pred
      pred_cells <- prediction$pred_cells
      # add prediction to meta data: 2 cols: p, class 
      for (colname in colnames(pred)) {
        obj.chunk[[colname]] <- unlist(
          lapply(colnames(obj.chunk), function(x)
            if (x %in% rownames(pred)) {pred[x, colname]}
            else {NA})
        )
      }
    }
    
    if (any(pred_cells != "")) {
      pred_cells <- gsub("/$", "", pred_cells)
      
      # double check if there is more than one predicted cell type
      if (ignore_ambiguous_result == TRUE) {
        pred_cells <- unlist(
          lapply(pred_cells, function(x) 
            if (length(unlist(strsplit(x, split = '/'))) <= 1) {x} 
            else {'ambiguous'})
          )
      } 
    
      # add cell type to meta data
      obj.chunk$predicted_cell_type <- pred_cells
      # this will be ignored if ignore ambiguous result is on
      if (ignore_ambiguous_result == FALSE) 
        obj.chunk$most_probable_cell_type <- simplify_prediction(
          as.matrix(SummarizedExperiment::colData(obj.chunk)), 
          as.matrix(mat.chunk), classifiers
        )
    }
    
    if (i == 1) classified_obj <- obj.chunk
    else classified_obj <- cbind(classified_obj, obj.chunk)
  }
  
  if (!is.null(cluster_slot) & !ignore_ambiguous_result) {
    clusts <- SummarizedExperiment::colData(classified_obj)[, cluster_slot]
    classified_obj$clust_pred <- 
      classify_clust(clusts, classified_obj$most_probable_cell_type)
  }
  
  return(classified_obj)
})
