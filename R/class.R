# scAnnotatR class definition ----

setOldClass("train")

#' scAnnotatR class. 
#' 
#' This class is returned by the \code{\link{train_classifier}} function.
#' Generally, scAnnotatR objects are never created directly.
#' 
#' @slot cell_type character. Name of the cell type.
#' @slot clf list. Trained model returned by caret train function.
#' @slot marker_genes vector/character containing marker genes 
#' used for the training.
#' @slot p_thres numeric. 
#' Probability threshold for the cell type to be assigned for a cell.
#' @slot parent character. Parent cell type.
#' @return A scAnnotatR object.
#' @import methods
#' @examples
#' # load small example dataset
#' data("tirosh_mel80_example")
#' 
#' # train a classifier, for ex: B cell
#' selected_marker_genes_B = c("CD19", "MS4A1", "CD79A")
#' set.seed(123)
#' clf_b <- train_classifier(train_obj = tirosh_mel80_example, 
#'                           marker_genes = selected_marker_genes_B, 
#'                           cell_type = "B cells")
#'
#' clf_b
#' @export
scAnnotatR <- setClass("scAnnotatR",
                       slots = list(cell_type = "character", 
                       clf = "train", 
                       marker_genes = "character", 
                       p_thres = "numeric",
                       parent = "character"))

# ---- constructor function

#' @param cell_type character. Name of the cell type.
#' @param clf list. Trained model returned by caret train function.
#' @param marker_genes vector/character containing marker genes used for the training.
#' @param p_thres numeric. 
#' Probability threshold for the cell type to be assigned for a cell.
#' @param parent character. Parent cell type.
#' @export
scAnnotatR <- function(cell_type, clf, marker_genes, p_thres, parent) {
    classifier <- methods::new("scAnnotatR",
                            cell_type = cell_type,
                            clf = clf,
                            marker_genes = marker_genes,
                            p_thres = p_thres,
                            parent = parent)
    return(classifier)
}

#' Internal functions of scAnnotatR package
#'
#' Check if a scAnnotatR object is valid
#'
#' @param object The request classifier to check.
#'
#' @return TRUE if the classifier is valid or the reason why it is not
#' @rdname internal
#' 
checkObjectValidity <- function(object) {
    # check cell_type
    cell_type.val <- checkCellTypeValidity(cell_type(object))
    if (is.character(cell_type.val)) {
    return(cell_type.val)
    }
  
    # check clf
    clf.val <- checkClassifierValidity(clf(object))
    if (is.character(clf.val)) {
      return(clf.val)
    }
    
    # check marker_genes
    marker_genes.val <- checkMarkerGenesValidity(marker_genes(object))
    if (is.character(marker_genes.val)) {
      return(marker_genes.val)
    }
    
    # check p_thres
    p_thres.val <- checkPThresValidity(p_thres(object))
    if (is.character(p_thres.val)) {
      return(p_thres.val)
    }
    
    # check parent
    parent.val <- checkParentValidity(parent(object))
    if (is.character(parent.val)) {
      return(parent.val)
    }
    
    return(TRUE)
}

#' Check validity of classifier cell type.
#'
#' @param cell_type Classifier cell type to check.
#'
#' @return TRUE if the cell type is valid or the reason why it is not.
#' @rdname internal
checkCellTypeValidity <- function(cell_type) {
  # cell_type must be a string
  if (!is(cell_type, "character")) {
    return("'cell_type' must be of class 'character'")
  }
  
  # cell_type must only contain one element
  if (length(cell_type) != 1) {
    return("'cell_type' can contain only one string.")
  }
  
  # make sure that method is not empty
  if (nchar(cell_type) < 1) {
    return("'cell_type' must be set")
  }
  
  return(TRUE)
}

#' Check validity of classifier marker_genes.
#'
#' @param marker_genes Classifier marker genes to check.
#'
#' @return TRUE if the marker_genes is valid or the reason why it is not.
#' @rdname internal
checkMarkerGenesValidity <- function(marker_genes) {
  # marker_genes must be a vector
  if (!is(marker_genes, "character")) {
    return("'marker_genes' must be of class 'character'")
  }
  
  # marker_genes must only contain at least one element
  if (length(marker_genes) < 1) {
    return("'marker_genes' must be a vector containing at least one feature")
  }
  
  # marker_genes must not be empty
  if (any(nchar(marker_genes) < 1)) {
    return("'marker_genes' must be set")
  }
  
  return(TRUE)
}

#' Check validity of classifier parent.
#'
#' @param parent Classifier parent to check.
#'
#' @return TRUE if the parent is valid or the reason why it is not.
#' @rdname internal
checkParentValidity <- function(parent) {
  # parent must be a string/vector
  if (!is(parent, "character")) {
    return("'parent' must be of class 'character'")
  }
  
  # parent must not be empty
  if (!is.na(parent)) {
    if (length(parent) == 1 && nchar(parent) < 1) {
      return("'parent' can be NA but not empty string")
    }
  }
  
  return(TRUE)
}

#' Check validity of classifier probability threshold.
#'
#' @param p_thres Classifier probability threshold to check.
#'
#' @return TRUE if the p_thres is valid or the reason why it is not.
#' @rdname internal
checkPThresValidity <- function(p_thres) {
  # p_thres must be a numeric
  if (!is(p_thres, "numeric")) {
    return("'p_thres' must be of class 'numeric'")
  }
  
  # p_thres must > 0
  if (p_thres <= 0) {
    return("'p_thres' must be positive")
  }
  
  return(TRUE)
}

#' Check validity of classifier
#'
#' @param clf Classifier to check.
#'
#' @return TRUE if the classifier is valid or the reason why it is not.
#' @rdname internal
checkClassifierValidity <- function(clf) {
  # clf must be a list
  if (!is.list(clf)) {
    return("'clf' must be of class 'list'")
  }
  
  # make sure that the object is returned from caret train func
  
  return(TRUE)
}

setValidity("scAnnotatR", checkObjectValidity)

#' Show object
#' 
#' @param object scAnnotatR object
#' 
#' @return print to console information about the object
#' 
#' @examples
#' data("tirosh_mel80_example")
#' selected_marker_genes_B = c("CD19", "MS4A1", "CD79A")
#' set.seed(123)
#' clf_b <- train_classifier(train_obj = tirosh_mel80_example, 
#' marker_genes = selected_marker_genes_B, cell_type = "B cells")
#' clf_b
#' 
#' @export
#' @rdname show
setMethod("show", c("object" = "scAnnotatR"), function(object) {
  cat("An object of class scAnnotatR for", cell_type(object), "\n")
  cat("*", toString(length(marker_genes(object))), "marker genes applied:", 
                     paste(marker_genes(object), collapse = ', '), "\n")
  cat("* Predicting probability threshold:", toString(p_thres(object)), "\n")
  if (!is.na(parent(object)) && length(parent(object)) == 1) {
    cat("* A child model of:", parent(object), "\n")
  } else {
    cat("* No parent model\n")
  }
})

#--- getters

#' cell_type
#' 
#' Returns the cell type for the given classifier.
#' 
#' @param classifier \code{\link{scAnnotatR}} object
#' 
#' @return cell type of object
#' 
#' @examples
#' data("tirosh_mel80_example")
#' selected_marker_genes_B = c("CD19", "MS4A1", "CD79A")
#' set.seed(123)
#' clf_b <- train_classifier(train_obj = tirosh_mel80_example, 
#' cell_type = "B cells", marker_genes = selected_marker_genes_B)
#' cell_type(clf_b)
#' 
#' @export
cell_type <- function(classifier) {
  return(classifier@cell_type) 
}

#' clf
#' 
#' Returns the classifier of the \code{\link{scAnnotatR}} object
#' 
#' @param classifier \code{\link{scAnnotatR}} object
#' 
#' @return Classifier is the object returned by caret SVM learning process.
#' More information about the caret package: https://topepo.github.io/caret/
#' 
#' @examples
#' data("tirosh_mel80_example")
#' selected_marker_genes_B = c("CD19", "MS4A1", "CD79A")
#' set.seed(123)
#' clf_b <- train_classifier(train_obj = tirosh_mel80_example, 
#' marker_genes = selected_marker_genes_B, cell_type = "B cells")
#' clf(clf_b)
#'  
#' @export
clf <- function(classifier) {
  return(classifier@clf)
}

#' marker_genes
#' 
#' Returns the set of marker genes for the given classifier.
#' 
#' @param classifier scAnnotatR object
#' 
#' @return Applied marker genes of object
#' 
#' @examples
#' data("tirosh_mel80_example")
#' selected_marker_genes_B = c("CD19", "MS4A1", "CD79A")
#' set.seed(123)
#' clf_b <- train_classifier(train_obj = tirosh_mel80_example, 
#' marker_genes = selected_marker_genes_B, cell_type = "B cells")
#' marker_genes(clf_b)
#' 
#' @export
marker_genes <- function(classifier) {
  return(classifier@marker_genes)
}

#' p_thres
#' 
#' Returns the probability threshold for the given classifier.
#' 
#' @param classifier scAnnotatR object
#' 
#' @return Predicting probability threshold of object
#' 
#' @examples
#' data("tirosh_mel80_example")
#' selected_marker_genes_B = c("CD19", "MS4A1", "CD79A")
#' set.seed(123)
#' clf_b <- train_classifier(train_obj = tirosh_mel80_example, 
#' marker_genes = selected_marker_genes_B, cell_type = "B cells")
#' p_thres(clf_b)
#' 
#' @export
#' 
p_thres <- function(classifier) {
  return (classifier@p_thres)
}

#' parent
#' 
#' Returns the parent of the cell type corresponding to the given classifier.
#' 
#' @param classifier scAnnotatR object
#' 
#' @return Parent model of object
#' 
#' @examples
#' data("tirosh_mel80_example")
#' selected_marker_genes_B = c("CD19", "MS4A1", "CD79A")
#' set.seed(123)
#' clf_b <- train_classifier(train_obj = tirosh_mel80_example, 
#' marker_genes = selected_marker_genes_B, cell_type = "B cells")
#' parent(clf_b)
#' 
#' @export
parent <- function(classifier) {
  return(classifier@parent)
}

#--- setters

#' Setter for cell_type.
#' Change cell type of a classifier
#' 
#' @param classifier the classifier whose cell type will be changed
#' @param value the new cell type
#' 
#' @return the classifier with the new cell type
#' @export
setGeneric('cell_type<-', function(classifier, value) 
  standardGeneric("cell_type<-"))

#' @inherit cell_type<-
#' @param classifier scAnnotatR object. 
#' The object is returned from the train_classifier function.
#' @return scAnnotatR object with the new cell type
#' @examples
#' data("tirosh_mel80_example")
#' selected_marker_genes_B = c("CD19", "MS4A1", "CD79A")
#' set.seed(123)
#' clf_b <- train_classifier(train_obj = tirosh_mel80_example, 
#' marker_genes = selected_marker_genes_B, cell_type = "B cells")
#' cell_type(clf_b) <- "B cell"
#' @rdname cell_type
setReplaceMethod('cell_type', c("classifier" = "scAnnotatR"), 
                 function(classifier, value) {
  # check if new thres is a string
  if (is.character(value) && nchar(value) > 0 && length(value) == 1)
    classifier@cell_type <- value
  else 
    stop("New cell type must be a non empty string.", call. = FALSE)
  
  # return or not?
  classifier
})

#' Setter for predicting probability threshold
#' 
#' @param classifier the classifier whose predicting probability threshold 
#' will be changed
#' @param value the new threshold
#' 
#' @return classifier with the new threshold.
#' @export
setGeneric('p_thres<-', function(classifier, value) 
  standardGeneric("p_thres<-"))

#' @inherit p_thres<-
#' @param classifier scAnnotatR object. 
#' The object is returned from the train_classifier function.
#' @return scAnnotatR object with the new threshold.
#' @examples
#' data("tirosh_mel80_example")
#' selected_marker_genes_B = c("CD19", "MS4A1", "CD79A")
#' set.seed(123)
#' clf_b <- train_classifier(train_obj = tirosh_mel80_example, 
#' marker_genes = selected_marker_genes_B, cell_type = "B cells")
#' clf_b_test <- test_classifier(test_obj = tirosh_mel80_example, 
#' classifier = clf_b)
#' # assign a new threhold probability for prediction
#' p_thres(clf_b) <- 0.4
#' @rdname p_thres
setReplaceMethod('p_thres', c("classifier" = "scAnnotatR"), 
                 function(classifier, value) {
  # check if new thres > 0
  if (is.numeric(value) && value > 0)
    classifier@p_thres <- value
  else 
    stop("New threshold is not a positive number.", call. = FALSE)
  
  # return or not?
  classifier
})

#' Setter for parent
#' 
#' @param classifier the classifier whose parent will be changed
#' @param value the new parent
#' 
#' @return the classifier with the new parent.
#' @rdname internal
setGeneric('parent<-', function(classifier, value) 
  standardGeneric("parent<-"))

#' @inherit parent<-
#' @param classifier scAnnotatR object. 
#' The object is returned from the train_classifier function.
#' @return scAnnotatR object with the new parent
#' 
#' @rdname internal
#' 
setReplaceMethod('parent', c("classifier" = "scAnnotatR"), 
                 function(classifier, value) {
  # check if new thres > 0
  if (!is.character(value) || nchar(value) == 0 || length(value) != 1)
    stop("New parent must be a non empty string.", call. = FALSE)
  else
    classifier@parent <- value
    
  # return or not?
  classifier
})

#' Setter for clf.
#' Change of clf will also lead to change of marker_genes.
#' @param classifier the classifier whose classifying model will be changed
#' @param value the new core model
#' 
#' @return the classifier with the new core model.
#' @rdname internal
setGeneric('clf<-', function(classifier, value) 
  standardGeneric("clf<-"))

#' @inherit clf<-
#' @param classifier scAnnotatR object. 
#' The object is returned from the train_classifier function.
#' @param value the new classifier
#' 
#' @return scAnnotatR object with the new trained classifier.
#' @rdname internal
setReplaceMethod('clf', c("classifier" = "scAnnotatR"), 
                 function(classifier, value) {
  # set new classifier
  if (is.na(parent(classifier))) {
    classifier@clf <- value
    
    # set new marker_genes
    new_marker_genes <- labels(value$terms)
    # convert underscore to hyphen if exists
    new_marker_genes <- gsub('_', '-', new_marker_genes) 
    marker_genes(classifier) <- new_marker_genes
  } else {
    stop("Can only assign new classifier for a cell type that has no parent.
    For a sub cell type: train a new classifier based on parent classifier.", 
         call. = FALSE)
  }
  
  # return or not?
  classifier
})

#' Setter for marker_genes. Users are not allowed to change marker_genes. 
#' 
#' @param classifier the classifier whose marker genes will be changed
#' @param value the new list of marker genes
#' 
#' @return the classifier with the new marker genes
#' @rdname internal
setGeneric('marker_genes<-', function(classifier, value) 
  standardGeneric("marker_genes<-"))

#' @inherit marker_genes<-
#' @param classifier scAnnotatR object. 
#' The object is returned from the train_classifier function.
#' @param value the new classifier
#' 
#' @return scAnnotatR object with the new marker genes.
#' @rdname internal
#' 
setReplaceMethod('marker_genes', c("classifier" = "scAnnotatR"), 
                 function(classifier, value) {
  # set new marker_genes
  if (is.character(value) && any(nchar(value)) > 0)
    classifier@marker_genes <- value
  
  # return or not?
  classifier
})
