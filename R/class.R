# scTypeR class definition ----

setOldClass("train")

#' scTypeR class. 
#' 
#' This class is returned by the \code{\link{train_classifier}} function.
#' Generally, scTypeR objects are never created directly.
#' 
#' @slot cell_type character. Name of the cell type.
#' @slot clf list. Trained model returned by caret train function.
#' @slot features vector/character containing features 
#' used for the training.
#' @slot p_thres numeric. 
#' Probability threshold for the cell type to be assigned for a cell.
#' @slot parent character. Parent cell type.
#' @return A scTypeR object.
#' @import methods
#' @examples
#' # load small example dataset
#' data("tirosh_mel80_example")
#' 
#' # train a classifier, for ex: B cell
#' selected_features_B = c("CD19", "MS4A1", "CD79A")
#' set.seed(123)
#' clf_b <- train_classifier(train_obj = tirosh_mel80_example, 
#'                           features = selected_features_B, cell_type = "B cells")
#'
#' clf_b
scTypeR <- setClass("scTypeR",
                            slots = list(cell_type = "character", 
                            clf = "train", 
                            features = "character", 
                            p_thres = "numeric",
                            parent = "character"))

# ---- constructor function

#' @param cell_type character. Name of the cell type.
#' @param clf list. Trained model returned by caret train function.
#' @param features vector/character containing features used for the training.
#' @param p_thres numeric. 
#' Probability threshold for the cell type to be assigned for a cell.
#' @param parent character. Parent cell type.
scTypeR <- function(cell_type, clf, features, p_thres, parent) {
    classifier <- methods::new("scTypeR",
                            cell_type = cell_type,
                            clf = clf,
                            features = features,
                            p_thres = p_thres,
                            parent = parent)
    return(classifier)
}

#' Internal functions of scTypeR package
#'
#' Check if a scTypeR object is valid
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
  
  # check features
  features.val <- checkFeaturesValidity(features(object))
  if (is.character(features.val)) {
    return(features.val)
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

#' Check validity of classifier features.
#'
#' @param features Classifier features to check.
#'
#' @return TRUE if the features is valid or the reason why it is not.
#' @rdname internal
checkFeaturesValidity <- function(features) {
  # features must be a vector
  if (!is(features, "character")) {
    return("'features' must be of class 'character'")
  }
  
  # features must only contain at least one element
  if (length(features) < 1) {
    return("'features' must be a vector containing at least one feature")
  }
  
  # features must not be empty
  if (any(nchar(features) < 1)) {
    return("'features' must be set")
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

setValidity("scTypeR", checkObjectValidity)

#' Show object
#' 
#' @param object scTypeR object
#' 
#' @return print to console information about the object
#' 
#' @examples
#' data("tirosh_mel80_example")
#' selected_features_B = c("CD19", "MS4A1", "CD79A")
#' set.seed(123)
#' clf_b <- train_classifier(train_obj = tirosh_mel80_example, 
#' features = selected_features_B, cell_type = "B cells")
#' clf_b
#' 
#' @export
#' @rdname show
setMethod("show", c("object" = "scTypeR"), function(object) {
  cat("An object of class scTypeR for ", cell_type(object), "\n")
  cat("* ", toString(length(features(object))), " features applied: ", 
                     paste(features(object), collapse = ', '), "\n")
  cat("* Predicting probability threshold: ", toString(p_thres(object)), "\n")
  if (!is.na(parent(object)) && length(parent(object)) == 1) {
    cat("* A child model of: ", parent(object), "\n")
  } else {
    cat("* No parent model\n")
  }
})

#--- getters

#' cell_type
#' 
#' Returns the cell type for the given classifier.
#' 
#' @param classifier \code{\link{scTypeR}} object
#' 
#' @return cell type of object
#' 
#' @examples
#' data("tirosh_mel80_example")
#' selected_features_B = c("CD19", "MS4A1", "CD79A")
#' set.seed(123)
#' clf_b <- train_classifier(train_obj = tirosh_mel80_example, 
#' cell_type = "B cells", features = selected_features_B)
#' cell_type(clf_b)
#' 
#' @export
setGeneric("cell_type", function(classifier) standardGeneric("cell_type"))

#' @inherit cell_type
#' @rdname cell_type
setMethod("cell_type", c("classifier" = "scTypeR"), 
          function(classifier) classifier@cell_type)

#' clf
#' 
#' Returns the classifier of the \code{\link{scTypeR}} object
#' 
#' @param classifier \code{\link{scTypeR}} object
#' 
#' @return Classifier is the object returned by caret SVM learning process.
#' More information about the caret package: https://topepo.github.io/caret/
#' 
#' @examples
#' data("tirosh_mel80_example")
#' selected_features_B = c("CD19", "MS4A1", "CD79A")
#' set.seed(123)
#' clf_b <- train_classifier(train_obj = tirosh_mel80_example, 
#' features = selected_features_B, cell_type = "B cells")
#' clf(clf_b)
#'  
#' @export
setGeneric("clf", function(classifier) standardGeneric("clf"))

#' @inherit clf
#' @rdname clf
setMethod("clf", c("classifier" = "scTypeR"), 
          function(classifier) classifier@clf)

#' features
#' 
#' Returns the set of features for the given classifier.
#' 
#' @param classifier scTypeR object
#' 
#' @return Applied features of object
#' 
#' @examples
#' data("tirosh_mel80_example")
#' selected_features_B = c("CD19", "MS4A1", "CD79A")
#' set.seed(123)
#' clf_b <- train_classifier(train_obj = tirosh_mel80_example, 
#' features = selected_features_B, cell_type = "B cells")
#' features(clf_b)
#' 
#' @export
#' 
setGeneric("features", function(classifier) standardGeneric("features"))

#' @inherit features
#' 
#' @rdname features
setMethod("features", c("classifier" = "scTypeR"), 
          function(classifier) classifier@features)

#' p_thres
#' 
#' Returns the probability threshold for the given classifier.
#' 
#' @param classifier scTypeR object
#' 
#' @return Predicting probability threshold of object
#' 
#' @examples
#' data("tirosh_mel80_example")
#' selected_features_B = c("CD19", "MS4A1", "CD79A")
#' set.seed(123)
#' clf_b <- train_classifier(train_obj = tirosh_mel80_example, 
#' features = selected_features_B, cell_type = "B cells")
#' p_thres(clf_b)
#' 
#' @export
#' 
setGeneric("p_thres", function(classifier) standardGeneric("p_thres"))

#' @inherit p_thres
#' 
#' @rdname p_thres
setMethod("p_thres", c("classifier" = "scTypeR"), 
          function(classifier) classifier@p_thres)

#' parent
#' 
#' Returns the parent of the cell type corresponding to the given classifier.
#' 
#' @param classifier scTypeR object
#' 
#' @return Parent model of object
#' 
#' @examples
#' data("tirosh_mel80_example")
#' selected_features_B = c("CD19", "MS4A1", "CD79A")
#' set.seed(123)
#' clf_b <- train_classifier(train_obj = tirosh_mel80_example, 
#' features = selected_features_B, cell_type = "B cells")
#' parent(clf_b)
#' 
#' @export
setGeneric("parent", function(classifier) standardGeneric("parent"))

#' @inherit parent
#' @rdname parent
setMethod("parent", c("classifier" = "scTypeR"), 
          function(classifier) classifier@parent)

#--- setters

#' Setter for cell_type
#' Change cell type for a classifier
#' 
#' @param classifier scTypeR object. 
#' The object is returned from the train_classifier function.
#' @param value the new cell type
#' 
#' @return scTypeR object with the new cell type.
#' @examples
#' data("tirosh_mel80_example")
#' selected_features_B = c("CD19", "MS4A1", "CD79A")
#' set.seed(123)
#' clf_b <- train_classifier(train_obj = tirosh_mel80_example, 
#' features = selected_features_B, cell_type = "B cells")
#' cell_type(clf_b) <- "B cell"
#' @export
setGeneric("cell_type<-", function(classifier, value) 
  standardGeneric("cell_type<-"))

#' @inherit cell_type<-
#' @rdname cell_type-set
setMethod("cell_type<-", c("classifier" = "scTypeR"), 
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
#' @param classifier scTypeR object. 
#' The object is returned from the train_classifier function.
#' @param value the new threshold
#' 
#' @return scTypeR object with the new threshold.
#' @examples
#' data("tirosh_mel80_example")
#' selected_features_B = c("CD19", "MS4A1", "CD79A")
#' set.seed(123)
#' clf_b <- train_classifier(train_obj = tirosh_mel80_example, 
#' features = selected_features_B, cell_type = "B cells")
#' clf_b_test <- test_classifier(test_obj = tirosh_mel80_example, 
#' classifier = clf_b)
#' # assign a new threhold probability for prediction
#' p_thres(clf_b) <- 0.4
#' @export
setGeneric("p_thres<-", function(classifier, value) 
  standardGeneric("p_thres<-"))

#' @inherit p_thres<-
#' @rdname p_thres-set
setMethod("p_thres<-", c("classifier" = "scTypeR"), 
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
#' @param classifier scTypeR object. 
#' The object is returned from the train_classifier function.
#' @param value the new parent
#' 
#' @return scTypeR object with the new parent.
#' @rdname internal
#' 
setGeneric("parent<-", function(classifier, value) standardGeneric("parent<-"))

#' @inherit parent<-
#' 
#' @rdname internal
setMethod("parent<-", c("classifier" = "scTypeR"), 
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
#' Change of clf will also lead to change of features.
#' 
#' @param classifier scTypeR object. 
#' The object is returned from the train_classifier function.
#' @param value the new classifier
#' 
#' @return scTypeR object with the new trained classifier.
#' @rdname internal
setGeneric("clf<-", function(classifier, value) standardGeneric("clf<-"))

#' @inherit clf<-
#' 
#' @rdname internal
#' 
setMethod("clf<-", c("classifier" = "scTypeR"), 
          function(classifier, value) {
  # set new classifier
  if (is.na(parent(classifier))) {
    classifier@clf <- value
    
    # set new features
    new_features <- labels(value$terms)
    features(classifier) <- new_features
  } else {
    stop("Can only assign new classifier for a cell type that has no parent.
    For a sub cell type: train a new classifier based on parent classifier.", 
         call. = FALSE)
  }
  
  # return or not?
  classifier
})

#' Setter for features. Users are not allowed to change features. 
#' 
#' @param classifier scTypeR object. 
#' The object is returned from the train_classifier function.
#' @param value the new classifier
#' 
#' @return scTypeR object with the new features.
#' @rdname internal
#' 
setGeneric("features<-", function(classifier, value) 
  standardGeneric("features<-"))

#' @inherit features<-
#' @rdname internal
#' 
setMethod("features<-", c("classifier" = "scTypeR"), 
          function(classifier, value) {
  # set new features
  if (is.character(value) && any(nchar(value)) > 0)
    classifier@features <- value
  
  # return or not?
  classifier
})
