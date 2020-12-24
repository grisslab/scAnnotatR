#' Balance training dataset
#'
#' @param mat count matrix of dimension m x n,
#' corresponding to m cells and n features
#' @param tag named list of training tags/labels (yes/no)
#' corresponding to a specific cell type, name and length of
#' list must be coherent with cells in mat
#' 
#' @return a list of balanced count matrix
#' and corresponding tags of balanced count matrix
#' @rdname internal
balance_dataset <- function(mat, tag) {
  message("Imbalanced dataset has: ", toString(nrow(mat)), " cells.")
  
  n_pos = length(tag[tag == 'yes'])
  n_neg = length(tag[tag == 'no'])
  
  if (n_pos >= n_neg) {
    cut_val = 'yes'
    n_cut = n_neg
  } else {
    cut_val = 'no'
    n_cut = n_pos
  }
  
  # get list of index that need to be cut off
  cut_idx = names(tag[tag == cut_val])
  
  # random a list for cut off observations
  random_idx = sample(cut_idx, n_cut)
  
  # subset n_cut observations
  cut_mat = mat[random_idx,, drop = FALSE]
  
  # subset the corresponding tag
  cut_tag = tag[random_idx]
  
  # refabricate the balanced dataset
  balanced_mat = rbind(
    cut_mat, mat[!rownames(mat) %in% cut_idx,, drop = FALSE])
  balanced_tag = append(cut_tag, tag[tag != cut_val])
  
  message("Balanced dataset has: ", toString(nrow(balanced_mat)), " cells.")
  
  return_val = list("mat" = balanced_mat, "tag" = balanced_tag)

  return(return_val)
}

#' Call training method
#' 
#' @param mat count matrix of dimension m x n
#' corresponding to m cells and n features. 
#' The matrix must have been balanced before. 
#' If not, pass it through the balance_dataset function.
#' @param tag named list of training tags/labels (yes/no) 
#' corresponding to a specific cell type, name and length of 
#' list must be coherent with cells in mat
#'  
#' @return the classification model (caret object)
#' 
#' @import caret
#' @import e1071
#' @import ape 
#' @rawNamespace import(kernlab, except = c(alpha, predict))
#' @rdname internal
train_func <- function(mat, tag) {
  
  # calculate sigma
  # calculate var of mat
  mat.vec <- as.vector(mat)
  mat.len <- length(mat.vec)
  mat.var <- var(mat.vec) * (mat.len - 1) / mat.len
  sigma <- 1 / (ncol(mat) * mat.var)

  mat <- as.data.frame(mat)
  mat$tag <- tag
  clf <- caret::train(form = tag ~ ., data = mat, 
                      method = "svmRadial",
                      tuneGrid = data.frame(
                        .C = 1, .sigma = sigma),
                      metrix = "Accuracy",
                      trControl = trainControl(
                        classProbs = TRUE, trim = TRUE,
                        returnData = FALSE, returnResamp = 'none')
                      ) 
  return(clf)
}

#' Transform whole matrix of counts to z-score
#' 
#' @param mat count matrix of dimension m x n 
#' corresponding to m cells and n features
#' 
#' @return row wise center-scaled count matrix
#' 
#' @importFrom stats sd
#' @rdname internal
transform_to_zscore <- function(mat) {
  # z_mat = scale(mat) # this cause NaN when column has zero variance
  z_mat <- apply(mat, 2, function(y) 
    (y - mean(y)) / stats::sd(y) ^ as.logical(stats::sd(y)))
  return(z_mat)
}

#' Load classifiers from databases
#' 
#' @param path_to_models path to databases, or by default
#' 
#' @return list of classifiers
#' 
#' @importFrom utils data
#' @rdname internal
load_models <- function(path_to_models) {
  # prevents R CMD check note
  model_list <- new_models <- default_models <- NULL
  
  if ("default" %in% path_to_models) {
    utils::data("default_models", envir = environment())
    model_list <- default_models
  } else {
    models_path <- paste0(path_to_models, "/new_models.rda")
    if (!file.exists(models_path)) {
      cat("No model found in provided path to models")
    } else {
      load(models_path, envir = environment()) 
      # models are stored in a variable called new_models
      model_list <- new_models
    }
  }
  
  return(model_list)
}

#' Perform features selection and handle missing features
#' 
#' @param mat count matrix of dimension n x m 
#' corresponding to m cells and n features
#' @param features list of selected features
#' 
#' @return filtered matrix
#' @rdname internal
select_features <- function(mat, features) {
  filtered_mat <- mat[rownames(mat) %in% features,, drop = FALSE]
  
  # perform features selection
  if (any(!features %in% rownames(filtered_mat))) { 
    # cannot perform features selection
    # cat("Not enough features for parent classifier. 
    # Cannot rerun pretrained classifier for parent cell type.\n")
    addi_features <- features[!features %in% rownames(filtered_mat)]
    zero_vec <- c(rep(0, ncol(filtered_mat)))
    for (feature in addi_features) {
      filtered_mat <- rbind(filtered_mat, zero_vec)
      rownames(filtered_mat)[nrow(filtered_mat)] <- feature
    }
  }
  
  return(filtered_mat)
}

#' Check label coherence in parent and child cell type
#' 
#' @param obj object
#' @param pos_parent a vector indicating parent clf prediction
#' @param parent_cell name of parent cell type
#' @param cell_type name of child cell type
#' @param target_cell_type alternative cell types (in case of testing clf)
#' @param ... arguments passed to other methods
#' 
#' @return list of adjusted object and adjusted tag slot
#' @rdname internal
setGeneric("check_parent_child_coherence", 
           function(obj, pos_parent, parent_cell, cell_type, 
                    target_cell_type, ...) 
             standardGeneric("check_parent_child_coherence"))

#' @inherit check_parent_child_coherence
#' 
#' @description Check label coherence in parent and 
#' child cell type in a \code{\link{Seurat}} object.
#' 
#' @param tag_slot tag slot in \code{\link{Seurat}} object 
#' indicating cell type
#' 
#' @importFrom Seurat Idents
#' 
#' @rdname check_parent_child_coherence
setMethod("check_parent_child_coherence", c("obj" = "Seurat"), 
          function(obj, pos_parent, parent_cell, cell_type, 
                   target_cell_type, tag_slot) {
  pos.val <- c(1, "yes", TRUE)
  
  # prepare (sub) cell type tag  
  if (tag_slot == "active.ident") {
    subtype <- Seurat::Idents(obj)
    test <- tolower(subtype) %in% tolower(target_cell_type) | 
      subtype %in% pos.val
    pos_subtype <- subtype[test]
  } else {
    subtype <- obj[[tag_slot]]
    test <- tolower(subtype[, 1]) %in% tolower(target_cell_type) | 
      subtype[, 1] %in% pos.val
    pos_subtype <- subtype[test, , drop=FALSE]
  }
  
  #-- compare with cell type with parent cell type, 
  # ie. cell, which is cell_type, must also be cell_parent
  # if not, raise warnings
  if (any(!rownames(pos_subtype) %in% pos_parent)) {
    warning("Some annotated ", cell_type, " are negative to ", 
    parent_cell, " classifier. They are removed from training/testing for ", 
    cell_type, " classifier.\n", call. = FALSE, immediate. = TRUE)
  }
  tag_slot <- "new_tag_slot"
  
  # join parent cell type and child cell type
  test <- 
    (colnames(obj) %in% rownames(pos_subtype)) & 
    (colnames(obj) %in% pos_parent)
  new.tag_slot <- ifelse(test, 'yes', 'no')
  new.tag_slot <- unlist(
    lapply(seq_len(length(new.tag_slot)), function(i) 
      if (!colnames(obj)[i] %in% pos_parent) {"not applicable"} 
      else {new.tag_slot[i]})
    )
  names(new.tag_slot) <- colnames(obj)
  obj[[tag_slot]] <- new.tag_slot
  
  return_val = list('adjusted_object' = obj, 'adjusted_tag_slot' = tag_slot)
  return(return_val)
})

#' @inherit check_parent_child_coherence
#' 
#' @description Check label coherence in parent and child cell type 
#' in a \code{\link{SingleCellExperiment}} object
#' 
#' @param tag_slot tag slot in \code{\link{SingleCellExperiment}} object 
#' indicating cell type
#' 
#' @importFrom SummarizedExperiment colData
#' 
#' @rdname check_parent_child_coherence
setMethod("check_parent_child_coherence", c("obj" = "SingleCellExperiment"), 
          function(obj, pos_parent, parent_cell, cell_type, 
                   target_cell_type, tag_slot) {
  pos.val <- c(1, "yes", TRUE)
  
  # prepare (sub) cell type tag  
  x <- SummarizedExperiment::colData(obj)[, tag_slot]
  subtype.bin <- ( tolower(x) %in% tolower(target_cell_type) | x %in% pos.val)
  pos_subtype.names <- colnames(obj)[subtype.bin]
  
  #-- compare with cell type with parent cell type, 
  # ie. cell, which is cell_type, must also be cell_parent
  # if not, raise warnings
  if (any(!pos_subtype.names %in% pos_parent)) {
    warning("Some annotated ", cell_type, " are negative to ", 
    parent_cell, " classifier. They are removed from training/testing for ", 
    cell_type, " classifier.\n", call. = FALSE, immediate. = TRUE)
  }
  tag_slot <- "new_tag_slot"
  
  # join parent cell type and child cell type
  new.tag_slot <- unlist(lapply(colnames(obj), function(x) 
    if (x %in% pos_subtype.names && x %in% pos_parent) {"yes"} else {"no"}))
  new.tag_slot <- unlist(lapply(seq_len(length(new.tag_slot)), function(i) 
    if (!colnames(obj)[i] %in% pos_parent) {"not applicable"} 
    else {new.tag_slot[i]}))
  SummarizedExperiment::colData(obj)[, tag_slot] <- new.tag_slot
  
  return_val = list('adjusted_object' = obj, 'adjusted_tag_slot' = tag_slot)
  return(return_val)
})

#' Filter cells from ambiguous chars and non applicable cells
#' 
#' @param obj object
#' @param tag_slot slot in cell meta data indicating cell type
#' @rdname internal
setGeneric("filter_cells", function(obj, tag_slot) 
  standardGeneric("filter_cells"))

#' @inherit filter_cells
#' 
#' @description Filter cells from ambiguous chars and 
#' non applicable cells in a \code{\link{Seurat}} object
#' 
#' @return adjusted \code{\link{Seurat}} object
#'
#' @rdname internal
setMethod("filter_cells", c("obj" = "Seurat"), function(obj, tag_slot) {
  # define characters usually included in ambiguous cell types
  # this is to avoid considering ambiguous cell types as negative cell_type
  ambiguous.chars <- c("/", ",", " -", " [+]", "[.]", "and", 
                       "or", "[(]" ,"[)]", "ambiguous")
  
  # only eliminate cell labels containing cell_type and ambiguous.chars
  if (tag_slot == "active.ident") {
    cell.tags <- Seurat::Idents(obj)
    cell.tags <- as.data.frame(cell.tags)
    rownames(cell.tags) <- names(Seurat::Idents(obj))
  } else {
    cell.tags <- obj[[tag_slot]]
  }
  
  # positive.cells <- rownames(cell.tags[cell.tags[, 1] %in% pos.val 
  # | tolower(cell.tags[, 1]) %in% tolower(cell_type),, drop = F])
  # positive cells must not contain ambiguous chars
  ambiguous.cells <- 
    rownames(cell.tags[grepl(paste(ambiguous.chars, collapse="|"), 
                             cell.tags[, 1]),, drop = FALSE])
  n.applicable.cells <- 
    rownames(cell.tags[grepl("not applicable", cell.tags[, 1]) 
                       | is.na(cell.tags[, 1]),, drop = FALSE])
  
  keeping.cells <- 
    colnames(obj)[!((colnames(obj) %in% ambiguous.cells) 
                    | colnames(obj) %in% n.applicable.cells)]
  obj <- subset(obj, cells = keeping.cells)
  
  return(obj)
})

#' @inherit filter_cells
#' 
#' @description Filter cells from ambiguous chars and non applicable cells
#' in a \code{\link{SingleCellExperiment}} object
#'  
#' @return adjusted \code{\link{SingleCellExperiment}} object
#' 
#' @importFrom SingleCellExperiment colData
#' 
#' @rdname internal
setMethod("filter_cells", c("obj" = "SingleCellExperiment"), 
          function(obj, tag_slot) {
  # define characters usually included in ambiguous cell types
  # this is to avoid considering ambiguous cell types as negative cell_type
  ambiguous.chars <- c("/", ",", " -", " [+]", "[.]", "and", 
                       "or", "[(]" ,"[)]", "ambiguous")
  
  # only eliminate cell labels containing cell_type and ambiguous.chars
  cell.tags <- SummarizedExperiment::colData(obj)[, tag_slot]
  
  # positive.cells <- rownames(cell.tags[cell.tags[, 1] %in% pos.val 
  # | tolower(cell.tags[, 1]) %in% tolower(cell_type),, drop = F])
  # positive cells must not contain ambiguous chars
  ambiguous <- grepl(paste(ambiguous.chars, collapse="|"), cell.tags)
  n.applicable <- (grepl("not applicable", cell.tags) | is.na(cell.tags))
  
  obj <- obj[, !(ambiguous | n.applicable)]
  
  return(obj)
})

#' Construct tag vector
#' 
#' @param obj object
#' @param cell_type name of cell type
#' @param ... arguments passed to other methods
#'
#' @importFrom stats setNames
#' 
#' @return a binary vector for cell tag
#' 
#' @rdname internal
setGeneric("construct_tag_vect", 
           function(obj, cell_type, ...) 
  standardGeneric("construct_tag_vect"))

#' @inherit construct_tag_vect
#' 
#' @description Construct a uniform tag vector for all forms of labels 
#' in a \code{\link{Seurat}} object
#'
#' @param tag_slot tag slot in \code{\link{Seurat}} object indicating cell type
#' 
#' @importFrom Seurat Idents
#' 
#' @rdname internal
setMethod("construct_tag_vect", c("obj" = "Seurat"), 
          function(obj, cell_type, tag_slot) {
  pos.val <- c(1, "yes", TRUE)
  
  # construct new tag
  if (tag_slot == "active.ident") 
    x <- Seurat::Idents(obj)
  else 
    x <- obj[[tag_slot]][, 1]
  
  test <- (x %in% pos.val) | (tolower(x) %in% tolower(cell_type))
  tag <- ifelse(test, "yes", "no")
  
  named_tag = setNames(tag, colnames(obj))
  
  return(named_tag)
})

#' @inherit construct_tag_vect
#' 
#' @description Construct a uniform tag vector for all forms of labels
#' in a \code{\link{SingleCellExperiment}} object
#' 
#' @param tag_slot tag slot in \code{\link{SingleCellExperiment}} object 
#' indicating cell type
#' 
#' @importFrom SummarizedExperiment colData
#' 
#' @rdname internal
setMethod("construct_tag_vect", c("obj" = "SingleCellExperiment"), 
          function(obj, cell_type, tag_slot) {
  pos.val <- c(1, "yes", TRUE)
  
  x <- SummarizedExperiment::colData(obj)[, tag_slot] 
  test <- (x %in% pos.val) | (tolower(x) %in% tolower(cell_type))
  tag <- ifelse(test, "yes", "no")
  
  named_tag = setNames(tag, colnames(obj))
  
  return(named_tag)
})

#' Process parent clf
#' 
#' @param obj object
#' @param parent_tag_slot string, name of annotation tag slot in object 
#' indicating pre-assigned/predicted parent cell type
#' @param parent_cell_type name of parent cell type
#' @param parent_clf \code{\link{scTypeR}} object corresponding 
#' to classification model for the parent cell type
#' @param path_to_models path to databases, or by default
#' @param zscore boolean indicating the transformation of gene expression 
#' in object to zscore or not
#' @param ... arguments passed to other methods
#' 
#' @return list of cells which are positive to parent clf
#' 
#' @importFrom stats predict
#' @import dplyr
#' 
#' @rdname internal
setGeneric("process_parent_clf", 
           function(obj, parent_tag_slot, parent_cell_type, parent_clf, 
                    path_to_models, zscore = TRUE, ...) 
             standardGeneric("process_parent_clf"))

#' @inherit process_parent_clf
#' 
#' @description Process parent classifier in a \code{\link{Seurat}} object
#' 
#' @param seurat_assay name of assay to use in \code{\link{Seurat}} object
#' @param seurat_slot type of expression data to use in 
#' \code{\link{Seurat}} object
#' 
#' @importFrom Seurat GetAssayData
#' @importFrom Seurat Idents
#' 
#' @rdname internal
setMethod("process_parent_clf", c("obj" = "Seurat"), 
          function(obj, parent_tag_slot, parent_cell_type, 
                   parent_clf, path_to_models, zscore = TRUE, 
                   seurat_assay, seurat_slot, ...) {
  pos_parent <- parent.clf <- . <- model_list <- NULL
  
  if (is.na(parent_cell_type) && !is.null(parent_clf))
    parent_cell_type <- cell_type(parent_clf)
  
  # if sub cell type is indicated
  if (!is.na(parent_cell_type)) {
    #-- apply parent cell classifier
    if (is.null(parent_clf)) {
      message("Parent classifier not provided. Try finding available model.")
      
      model_list <- load_models(path_to_models)
      
      if (parent_cell_type %in% names(model_list)) {
        parent.clf <- model_list[[parent_cell_type]]
      } else {
        message("No available model for parent cell type")
      }
    }
    else {
      parent.clf <- parent_clf
    }
    
    if (!is.null(parent.clf)) {
      message("Apply pretrained model for parent cell type.")
      
      # convert Seurat object to matrix
      mat = Seurat::GetAssayData(object = obj, 
                                 assay = seurat_assay, slot = seurat_slot)
      
      filtered_mat <- select_features(mat, features(parent.clf))
      
      filtered_mat <- t(as.matrix(filtered_mat))
      
      # transform mat to zscore values
      if (zscore == TRUE) {
        filtered_mat <- transform_to_zscore(filtered_mat)
      }
      
      # predict
      pred = stats::predict(clf(parent.clf), filtered_mat, type = "prob") %>% 
           dplyr::mutate('class' = apply(., 1, 
           function(x) if(x[1] >= p_thres(parent.clf)) {"yes"} else {"no"}))
      rownames(pred) <- rownames(filtered_mat)
      pos_parent <- rownames(pred[pred$class == "yes",])
    } else if (!is.null(parent_tag_slot)) { # try with predicted tag slot
      message("Parent classifier could not be applied. 
              Try with predicted/pre-assigned cell type.")
      if (parent_tag_slot == 'active.ident') {
        cell_type_anno <- Seurat::Idents(obj)
        pos_parent <- names(cell_type_anno[tolower(cell_type_anno) 
                                           == tolower(parent_cell_type)]) 
      } else {
        cell_type_anno <- obj[[parent_tag_slot]][, 1]
        pos_parent <- colnames(obj)[tolower(cell_type_anno) == 
                                      tolower(parent_cell_type)]
      }
    } else { # only parent cell type provided but no parent clf can be used
      stop("Neither parent classifier nor parent tag slot applied. 
           Parent cell type verification failed. 
           Please check parent classifier/parent tag slot
           or remove parent cell type information.", call. = FALSE)
    }
  }
  
  return_val <- list('pos_parent' = pos_parent, 'parent_cell' = parent_cell_type, 
                     'parent.clf' = parent.clf, 'model_list' = model_list)
  return(return_val)
})

#' @inherit process_parent_clf
#' 
#' @description Process parent classifier in 
#' a \code{\link{SingleCellExperiment}} object
#' 
#' @param sce_assay name of assay to use 
#' in \code{\link{SingleCellExperiment}} object
#' 
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment assay
#' 
#' @rdname internal
setMethod("process_parent_clf", c("obj" = "SingleCellExperiment"), 
          function(obj, parent_tag_slot, parent_cell_type, parent_clf, 
                   path_to_models, zscore = TRUE, sce_assay, ...) {
  pos_parent <- parent.clf <- . <- model_list <- NULL
  
  if (is.na(parent_cell_type) && !is.null(parent_clf))
    parent_cell_type <- cell_type(parent_clf)
  
  # if sub cell type is indicated
  if (!is.na(parent_cell_type)) {
    #-- apply parent cell classifier
    # get parent classifier
    if (is.null(parent_clf)) {
      message("Parent classifier not provided. Try finding available model.")
      
      model_list <- load_models(path_to_models)
      
      if (parent_cell_type %in% names(model_list)) {
        parent.clf <- model_list[[parent_cell_type]]
      } else {
        message("No available model for parent cell type")
      }
    }
    else {
      parent.clf <- parent_clf
    }
    
    if (!is.null(parent.clf)) {
      message("Apply pretrained model for parent cell type.\n")
      
      # convert Seurat object to matrix
      mat = SummarizedExperiment::assay(obj, sce_assay)
      
      filtered_mat <- select_features(mat, features(parent.clf))
      filtered_mat <- t(as.matrix(filtered_mat))
      
      # transform mat to zscore values
      if (zscore == TRUE) {
        filtered_mat <- transform_to_zscore(filtered_mat)
      }
      
      # predict
      pred = stats::predict(clf(parent.clf), filtered_mat, type = "prob") %>% 
        dplyr::mutate('class' = apply(., 1, 
        function(x) if(x[1] >= p_thres(parent.clf)) {"yes"} else {"no"}))
      rownames(pred) <- rownames(filtered_mat)
      pos_parent <- rownames(pred[pred$class == "yes",])
    } else if (!is.null(parent_tag_slot)) { # try with predicted tag slot
      message("Parent classifier could not be applied. 
              Try with predicted/pre-assigned cell type.")
      cell_type_anno <- colData(obj)[, parent_tag_slot]
      pos_parent <- colnames(obj)[tolower(cell_type_anno) 
                                  == tolower(parent_cell_type)]
    } else { 
      # only parent cell type provided but no parent clf/tag slot can be used
      stop("Neither parent classifier nor parent tag slot applied. 
           Parent cell type verification failed. 
           Please check parent classifier/parent 
           tag slot or remove parent cell type information.", call. = FALSE)
    }
  }
  
  return_val <- list('pos_parent' = pos_parent, 'parent_cell'= parent_cell_type,
                     'parent.clf' = parent.clf, 'model_list' = model_list)
  return(return_val)
})


#' Make prediction
#'
#' @param mat count matrix used for prediction
#' @param classifier classifier
#' @param pred_cells a whole prediction for all cells
#' @param ignore_ambiguous_result whether ignore ambigouous result
#' 
#' @return prediction
#' 
#' @import dplyr
#' @importFrom stats predict
#' 
#' @rdname internal
make_prediction <- function(mat, classifier, pred_cells, 
                            ignore_ambiguous_result = TRUE) {
  . <- NULL
  cells <- names(pred_cells)
  
  # predict
  pred = stats::predict(clf(classifier), mat, type = "prob") %>%
    dplyr::mutate('class' = apply(., 1, function(x)
      if(x[1] >= p_thres(classifier)) {"yes"} else {"no"}))
  rownames(pred) <- rownames(mat)
  
  # append a summary to whole predicted cell type
  pred_cells <- unlist(lapply(cells,
  function(i)
    if (i %in% rownames(pred) && pred[i, "class"] == "yes") {
      test <- 
        ignore_ambiguous_result == TRUE && 
        !is.na(parent(classifier)) && 
        gsub("/", "", pred_cells[i]) == parent(classifier)
      if (test)
        paste0(cell_type(classifier), "/")
      else
        paste0(pred_cells[i], cell_type(classifier), "/")
    }
    else { pred_cells[i] }))
  names(pred_cells) <- cells

  # remove no column and rename yes column to p
  pred$no <- NULL
  colnames(pred)[1] <- 'p'

  # add cell type to colnames
  colnames(pred) <- unlist(
    lapply(colnames(pred), function(x) 
      paste0(
        c(unlist(strsplit(cell_type(classifier), split = " ")), x), 
        collapse = "_")
      )
    )

  return_val <- list('pred' = pred, 'pred_cells' = pred_cells)
  return(return_val)
}

#' Simplify prediction
#'
#' @param meta.data cell meta data
#' @param mat expression mat
#' @param classifiers classifiers 
#' 
#' @return simplified prediction
#' 
#' @rdname internal
simplify_prediction <- function(meta.data, mat, classifiers) {
  if (is.null(names(classifiers)))
    names(classifiers) <- unlist(lapply(classifiers, function(x) cell_type(x)))
            
  # list of parents named by children
  parents <- unlist(lapply(classifiers, function(x) parent(x)))
  
  # create simplified result for the first level: no parent
  noparent_idx <- which(is.na(parents))
  noparent_cell <- names(noparent_idx)
  noparent_pcol <- paste0(gsub(' ', '_', noparent_cell), '_p')
  noparent_p <- meta.data[, noparent_pcol, drop = FALSE]
  
  # create most probable pred from cell types having no parent
  max_p <- colnames(noparent_p)[unlist(apply(noparent_p, 1, which.max))]
  max_clf <- gsub('_p', '', max_p)
  max_clf <- gsub('_', ' ', max_clf)
  simplified <- unlist(lapply(1:nrow(noparent_p), function(i) 
    if (noparent_p[i, max_p[i]] >= p_thres(classifiers[[max_clf[i]]])) 
    {max_clf[i]} else {'unknown'}))
  names(simplified) <- colnames(mat)
  
  # continue to deeper level: children
  simplified.copy <- NULL
  while (!identical(simplified, simplified.copy)) {
    simplified.copy <- simplified # copy simplified
    for (parent in unique(simplified)) {
      if (parent %in% parents) {
        children <- names(which(parents == parent))
        
        children_pcol <- paste0(gsub(' ', '_', children), '_p')
        # extract prediction probabilities of children
        children_p <- meta.data[simplified == parent, 
                                children_pcol, drop = FALSE]
        max_p <- colnames(children_p)[unlist(apply(children_p, 1, which.max))]
        names(max_p) <- rownames(children_p)
        max_child <- gsub('_p', '', max_p)
        max_child <- gsub('_', ' ', max_child)
        names(max_child) <- rownames(children_p)
        simplified <- unlist(lapply(names(simplified), function(i) 
          if (simplified[i] == parent 
              && children_p[i, max_p[i]] >= p_thres(classifiers[[max_child[i]]])) 
          {max_child[i]} else {simplified[i]})) # change simplified
        names(simplified) <- colnames(mat)
      }
    }
  }
  
  return(simplified)
}

#' Verify parent prediction
#'
#' @param mat expression matrix
#' @param classifier classifier
#' @param meta.data object meta data
#' 
#' @return applicable matrix
#' 
#' @rdname internal
verify_parent <- function(mat, classifier, meta.data) {
  pos_parent <- applicable_mat <- NULL
  
  # parent clf, if avai, always has to be applied before children clf.
  parent_slot <- paste0(
    c(unlist(strsplit(parent(classifier), split = " ")), "class"), 
    collapse = "_")
  if (parent_slot %in% colnames(meta.data)) {
    parent_pred <- meta.data[, parent_slot]
    pos_parent <- colnames(mat)[parent_pred == 'yes'] 
  } else {
    warning('Parent classifier of ', cell_type(classifier), 'cannot be applied.\n 
             Please list/save parent classifier before child(ren) classifier.\n
             Skip applying classification models for ', cell_type(classifier), 
             ' and its parent cell type.\n', call. = FALSE, immediate. = TRUE)
  }
  
  if (!is.null(pos_parent)) {
    applicable_mat <- mat[, colnames(mat) %in% pos_parent, drop = FALSE]
  } # else next
  
  return(applicable_mat)
}

#' Test clf performance
#'
#' @param mat expression matrix
#' @param classifier classifier
#' @param tag tag of data
#' 
#' @return clf performance
#' @import dplyr
#' @import pROC
#' @importFrom stats predict
#' 
#' @rdname internal
test_performance <- function(mat, classifier, tag) {
  overall.roc <- . <- NULL
  
  tag <- unlist(lapply(tag, function(x) if (x == 'yes') {1} else {0}))
  
  iter <- unique(sort(c(p_thres(classifier), seq(0.1, 0.9, by = 0.1))))
  
  # predict
  for (thres in iter) {
    test_pred = stats::predict(clf(classifier), mat, type = "prob") %>% 
      dplyr::mutate('class' = apply(., 1, function(x) 
        if(x[1] >= thres) {1} else {0}))
    rownames(test_pred) <- rownames(mat)
    
    # calculate TPR, FPR 
    pr <- ROCR::prediction(test_pred$class, tag)
    pe <- ROCR::performance(pr, "tpr", "fpr")
    roc.data <- data.frame(fpr=unlist(pe@x.values), tpr=unlist(pe@y.values))
    
    if (thres == p_thres(classifier)) {
      pred <- test_pred
      message('Current probability threshold: ', toString(p_thres(classifier)))
      # accuracy
      message(" ", "\t\tPositive", "\tNegative", "\tTotal")
      message("Actual", "\t\t", toString(length(tag[tag == 1])), 
              "\t\t", toString(length(tag[tag == 0])), 
              "\t\t", toString(length(tag)))
      message("Predicted", "\t", 
              toString(nrow(test_pred[test_pred$class == 1,])),
              "\t\t", toString(nrow(test_pred[test_pred$class == 0,])), 
              "\t\t", toString(nrow(test_pred)), "\n")
      count <- 0
      for (i in seq_len(length(tag))) { # can improve this later
        if (tag[i] == test_pred$class[i]) 
          count <- count + 1 
      }
      acc <- count/length(tag)
      message("Accuracy: ", toString(acc), "\n")
      message("Sensivity (True Positive Rate) for ", 
              cell_type(classifier), ": ", toString(roc.data[2, 2]))
      message("Specificity (1 - False Positive Rate) for ", 
              cell_type(classifier), ": ", toString(1 - roc.data[2, 1]))
    }
    
    # add new result to overall
    overall.roc <- rbind(overall.roc, 
                         c(thres, roc.data[2, 1], roc.data[2, 2]))
  }
  # calculate AUC
  roc_obj <- pROC::roc(tag, test_pred$yes, levels = c(0, 1), direction = "<")
  auc_obj = pROC::auc(roc_obj)
  message("Area under the curve: ", toString(auc_obj))
  
  colnames(overall.roc) <- c('p_thres', 'fpr', 'tpr')
  return_val = list("pred" = pred, "acc" = acc, "test_tag" = tag, 
                    "overall_roc" = overall.roc, 'auc' = auc_obj)
  return(return_val)
}