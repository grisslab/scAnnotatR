#' Balance training dataset
#'
#' @param mat count matrix of dimension m x n,
#' corresponding to m cells and n marker genes
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
#' corresponding to m cells and n marker genes. 
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
  # mat.vec <- as.vector(mat)
  # mat.len <- length(mat.vec)
  # mat.var <- var(mat.vec) * (mat.len - 1) / mat.len
  # sigma <- 1 / (ncol(mat) * mat.var)

  mat <- as.data.frame(mat)
  mat$tag <- tag
  caret_model <- caret::train(form = tag ~ ., data = mat, 
                      method = "svmLinear",
                      tuneGrid = data.frame(.C = 1),
                      metrix = "Accuracy",
                      trControl = trainControl(method = "cv", 
                        classProbs = TRUE, trim = TRUE, sampling = 'down',
                        returnData = FALSE, returnResamp = 'none')
                      ) 
  return(caret_model)
}

#' Transform whole matrix of counts to z-score
#' 
#' @param mat count matrix of dimension m x n 
#' corresponding to m cells and n marker genes
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
  model_list <- NULL
  data_env <- new.env(parent = emptyenv())
  
  if (path_to_models == "default") {
    utils::data("default_models", envir = data_env)
    model_list <- data_env[["default_models"]]
  } else {
    models_path <- file.path(path_to_models, "new_models.rda")
    if (!file.exists(models_path)) {
      cat("No model found in provided path to models")
    } else {
      load(models_path, envir = data_env) 
      # models are stored in a variable called new_models
      model_list <- data_env[["new_models"]]
    }
  }
  
  return(model_list)
}

#' Perform marker genes selection and handle missing marker genes
#' 
#' @param mat count matrix of dimension n x m 
#' corresponding to m cells and n marker genes
#' @param marker_genes list of selected marker genes
#' 
#' @return filtered matrix
#' @rdname internal
select_marker_genes <- function(mat, marker_genes) {
  filtered_mat <- mat[rownames(mat) %in% marker_genes,, drop = FALSE]
  
  # perform marker genes selection
  if (any(!marker_genes %in% rownames(filtered_mat))) { 
    # cannot perform marker genes selection
    addi_marker_genes <- marker_genes[!marker_genes %in% rownames(filtered_mat)]
    zero_vec <- c(rep(0, ncol(filtered_mat)))
    for (marker in addi_marker_genes) {
      filtered_mat <- rbind(filtered_mat, zero_vec)
      rownames(filtered_mat)[nrow(filtered_mat)] <- marker
    }
  }
  
  return(filtered_mat)
}

#' Check label coherence in parent and child cell type
#' 
#' @param mat expression matrix of size n x m, n: genes, m: cells
#' @param tag vector, named list indicating children cell tage
#' @param pos_parent a vector indicating parent classifier prediction
#' @param parent_cell name of parent cell type
#' @param cell_type name of child cell type
#' @param target_cell_type alternative cell types (in case of testing classifier)
#' 
#' @return list of adjusted tag
#' @rdname internal
setGeneric("check_parent_child_coherence", 
           function(mat, tag, pos_parent, parent_cell, cell_type, 
                    target_cell_type) 
             standardGeneric("check_parent_child_coherence"))

#' @inherit check_parent_child_coherence
#' 
#' @rdname internal
setMethod("check_parent_child_coherence", c("mat" = "dgCMatrix", 'tag' = 'vector'), 
          function(mat, tag, pos_parent, parent_cell, cell_type, 
                   target_cell_type) {
  pos.val <- c(1, "yes", TRUE)
  
  # prepare (sub) cell type tag  
  #x <- SummarizedExperiment::colData(obj)[, tag_slot]
  tag.bin <- (tolower(tag) %in% tolower(target_cell_type) | tag %in% pos.val)
  pos_subtype.names <- colnames(mat)[tag.bin]
  
  #-- compare with cell type with parent cell type, 
  # ie. cell, which is cell_type, must also be cell_parent
  # if not, raise warnings
  if (any(!pos_subtype.names %in% pos_parent)) {
    warning("Some annotated ", cell_type, " are negative to ", 
            parent_cell, " classifier. They are removed from training/testing for ", 
            cell_type, " classifier.\n", call. = FALSE, immediate. = TRUE)
  }
  # tag_slot <- "new_tag_slot"
  
  # join parent cell type and child cell type
  new_tag <- unlist(lapply(colnames(mat), function(x) 
    if (x %in% pos_subtype.names && x %in% pos_parent) {"yes"} else {"no"}))
  new_tag <- unlist(lapply(seq_len(length(new_tag)), function(i) 
    if (!colnames(mat)[i] %in% pos_parent) {"not applicable"} 
    else {new_tag[i]}))
  #SummarizedExperiment::colData(obj)[, tag_slot] <- new.tag_slot
  
  return(new_tag)
})

#' Filter cells from ambiguous chars and non applicable cells
#' Ambiguous characters includes: "/", ",", "-", "+", ".", "and", 
#' "or", "(", ")", "ambiguous"
#' 
#' @param mat expression matrix of size n x m, n: genes, m: cells 
#' @param tag named list indicating cell type
#' 
#' @return filtered matrix and corresponding tag
#' @rdname internal
setGeneric("filter_cells", function(mat, tag) 
  standardGeneric("filter_cells"))

#' @inherit filter_cells
#' 
#' @rdname internal
setMethod("filter_cells", c("mat" = "dgCMatrix", "tag" = "vector"), 
          function(mat, tag) {
            # define characters usually included in ambiguous cell types
            # this is to avoid considering ambiguous cell types as negative cell_type
            ambiguous.chars <- c("/", ",", " -", " [+]", "[.]", " and ", 
                                 " or ", "_or_", "-or-", "[(]" ,"[)]", "ambiguous")
            
            # only eliminate cell labels containing cell_type and ambiguous.chars
            ambiguous <- grepl(paste(ambiguous.chars, collapse="|"), tag)
            n.applicable <- (grepl("not applicable", tag) | is.na(tag))
            
            if (any(ambiguous))
              warning('Cell types containing "/", ",", "-", "+", ".", "and", "or", "(", ")", and "ambiguous" are considered as ambiguous. They are removed from training and testing.\n', 
                      call. = FALSE, immediate. = TRUE)
            #obj <- obj[, !(ambiguous | n.applicable)]
            mat <- mat[, !(ambiguous | n.applicable), drop = FALSE]
            tag <- tag[!(ambiguous | n.applicable)]
            
            filtered <- list('mat' = mat, 'tag' = tag)
            return(filtered)
          })

#' Construct tag vector
#' 
#' @param tag named list containing the tag
#' @param cell_type name of cell type
#'
#' @importFrom stats setNames
#' 
#' @return a binary vector for cell tag
#' 
#' @rdname internal
setGeneric("construct_tag_vect", 
           function(tag, cell_type) 
             standardGeneric("construct_tag_vect"))

#' @inherit construct_tag_vect
#' 
#' @rdname internal
setMethod("construct_tag_vect", c("tag" = "vector"), 
          function(tag, cell_type) {
            pos.val <- c(1, "yes", TRUE)
            
            # x <- SummarizedExperiment::colData(obj)[, tag_slot] 
            test <- (tag %in% pos.val) | (tolower(tag) %in% tolower(cell_type))
            new_tag <- ifelse(test, "yes", "no")
            
            named_tag = setNames(new_tag, names(tag))
            
            return(named_tag)
          })

#' Process parent classifier
#' 
#' @param mat expression matrix of size n x m, n: genes, m: cells
#' @param parent_tag vector, named list indicating pre-assigned/predicted 
#' parent cell type
#' @param parent_cell_type name of parent cell type
#' @param parent_classifier \code{\link{scAnnotatR}} object corresponding 
#' to classification model for the parent cell type
#' @param path_to_models path to databases, or by default
#' @param zscore boolean indicating the transformation of gene expression 
#' in object to zscore or not
#' 
#' @return list of cells which are positive to parent classifier
#' 
#' @importFrom stats predict
#' @import dplyr
#' 
#' @rdname internal
setGeneric("process_parent_classifier", 
           function(mat, parent_tag, parent_cell_type, parent_classifier, 
                    path_to_models, zscore = TRUE) 
             standardGeneric("process_parent_classifier"))

#' @inherit process_parent_classifier
#' 
#' @rdname internal
setMethod("process_parent_classifier", c("mat" = "dgCMatrix"), 
          function(mat, parent_tag, parent_cell_type, parent_classifier, 
                   path_to_models, zscore = TRUE) {
    pos_parent <- parent.classifier <- . <- model_list <- NULL
    
    if (is.na(parent_cell_type) && !is.null(parent_classifier))
      parent_cell_type <- cell_type(parent_classifier)
    
    # if sub cell type is indicated
    if (!is.na(parent_cell_type)) {
      #-- apply parent cell classifier
      # get parent classifier
      if (is.null(parent_classifier)) {
        message("Parent classifier not provided. Try finding available model.")
        
        model_list <- load_models(path_to_models)
        
        if (parent_cell_type %in% names(model_list)) {
          parent.classifier <- model_list[[parent_cell_type]]
        } else {
          message("No available model for parent cell type")
        }
      }
      else {
        parent.classifier <- parent_classifier
      }
      
      if (!is.null(parent.classifier)) {
        message("Apply pretrained model for parent cell type.\n")
        
        # convert Seurat object to matrix
        # mat = SummarizedExperiment::assay(obj, sce_assay)
        
        filtered_mat <- select_marker_genes(mat, marker_genes(parent.classifier))
        filtered_mat <- t(as.matrix(filtered_mat))
        
        # transform mat to zscore values
        if (zscore == TRUE) {
          filtered_mat <- transform_to_zscore(filtered_mat)
        }
        
        # to avoid problem triggered by '-' in gene names
        colnames(filtered_mat) <- gsub('-', '_', colnames(filtered_mat))
        
        # predict
        pred = stats::predict(caret_model(parent.classifier), filtered_mat, type = "prob") %>% 
          dplyr::mutate('class' = apply(., 1, 
                                        function(x) if(x[1] >= p_thres(parent.classifier)) {"yes"} else {"no"}))
        rownames(pred) <- rownames(filtered_mat)
        pos_parent <- rownames(pred[pred$class == "yes",])
      } else if (!is.null(parent_tag)) { # try with predicted tag slot
        message("Parent classifier could not be applied. 
      Try with predicted/pre-assigned cell type.")
        #cell_type_anno <- colData(obj)[, parent_tag_slot]
        pos_parent <- names(parent_tag)[tolower(parent_tag) 
                                    == tolower(parent_cell_type)]
      } else { 
        # only parent cell type provided but no parent classifier/tag slot can be used
        stop("Neither parent classifier nor parent tag slot applied. 
   Parent cell type verification failed. 
   Please check parent classifier/parent 
   tag slot or remove parent cell type information.", call. = FALSE)
      }
    }
    
    return_val <- list('pos_parent' = pos_parent, 'parent_cell'= parent_cell_type,
                       'parent.classifier' = parent.classifier, 'model_list' = model_list)
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
  
  # to avoid problem triggered by '-' in gene names
  colnames(mat) <- gsub('-', '_', colnames(mat))
  
  # predict
  pred = stats::predict(caret_model(classifier), mat, type = "prob") %>%
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
        paste0("/", cell_type(classifier))
      else
        paste0(pred_cells[i], "/", cell_type(classifier))
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
#' @param full_pred full prediction
#' @param classifiers classifiers 
#' 
#' @return simplified prediction
#' 
#' @rdname internal
simplify_prediction <- function(meta.data, full_pred, classifiers) {
  if (is.null(names(classifiers)))
    names(classifiers) <- unlist(lapply(classifiers, function(x) cell_type(x)))
  
  # list of parents named by children
  parents <- unlist(lapply(classifiers, function(x) parent(x)))
  simplified <- full_pred
  names(simplified) <- rownames(meta.data)
  
  # parent level
  for (cell in rownames(meta.data)) {
    predicted_types <- unlist(strsplit(full_pred[cell], split = '/'))
    #predicted_parents <- parents[parents %in% predicted_types]
    if (length(predicted_types) >= 2) {
      p.pcol.names <- paste0(gsub(' ', '_', predicted_types), '_p')
      p.prob <- meta.data[cell, p.pcol.names, drop = FALSE]
      simplified[cell] <- colnames(p.prob)[which.max(p.prob)]
    }
  }
  simplified <- gsub('_p$', '', simplified)
  simplified <- gsub('_', ' ', simplified)
  
  # continue to deeper level: children
  simplified.copy <- NULL
  while (!identical(simplified, simplified.copy)) {
    simplified.copy <- simplified # copy simplified
    for (cell in rownames(meta.data)) {
      parent <- simplified[cell]
      if (parent %in% parents){
        children <- names(which(parents == parent))
        predicted_types <- unlist(strsplit(full_pred[cell], split = '/'))
        predicted_children <- children[children %in% predicted_types]
        if (length(predicted_children) >= 2) {
          c.pcol.names <- paste0(gsub(' ', '_', predicted_children), '_p')
          c.prob <- meta.data[cell, c.pcol.names, drop = FALSE]
          simplified[cell] <- colnames(c.prob)[which.max(c.prob)]
        } else if (length(predicted_children) == 1) {
          simplified[cell] <- predicted_children
        } else simplified[cell] <- simplified[cell]
      }
    }
    simplified <- gsub('_p$', '', simplified)
    simplified <- gsub('_', ' ', simplified)
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
  
  # parent classifier, if avai, always has to be applied before children classifier.
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

#' Test classifier performance
#'
#' @param mat expression matrix
#' @param classifier classifier
#' @param tag tag of data
#' 
#' @return classifier performance
#' @import dplyr
#' @import pROC
#' @importFrom stats predict
#' 
#' @rdname internal
test_performance <- function(mat, classifier, tag) {
  overall.roc <- . <- NULL
  
  # to avoid problem triggered by '-' in gene names
  colnames(mat) <- gsub('-', '_', colnames(mat))
    
  tag <- unlist(lapply(tag, function(x) if (x == 'yes') {1} else {0}))
  
  iter <- unique(sort(c(p_thres(classifier), seq(0.1, 0.9, by = 0.1))))
  
  # predict
  for (thres in iter) {
    test_pred = stats::predict(caret_model(classifier), mat, type = "prob") %>% 
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

#' Test classifier performance
#'
#' @param clusts cluster info
#' @param most_probable_cell_type predicted cell type
#' 
#' @rdname internal
classify_clust <- function(clusts, most_probable_cell_type) {
  clust.cell.coor <- table(most_probable_cell_type, clusts) 
  max.val <- apply(clust.cell.coor, 2, function(x) max(x)/sum(x))
  names(max.val) <- 
    unname(apply(clust.cell.coor, 2, 
                 function(x) rownames(clust.cell.coor)[which.max(x)]))
  clust.pred <- paste0(round(max.val * 100, 2), '% ', names(max.val))
  names(clust.pred) <- levels(clusts)
  converted_pred <- unlist(lapply(clusts, function(x) clust.pred[[x]]))
  return(converted_pred)
}