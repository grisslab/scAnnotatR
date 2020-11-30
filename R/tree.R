#' Save a model to the package
#' 
#' @param new_model new model to be added into the classification tree
#' @param include.default whether include the default models of the package 
#' in the list of new trained models or not. 
#' If users further want to classify cells, they can only use default 
#' pretrained model list or their new model list. 
#' Including default models in new trained models helps users using 
#' both of them once. In addition, default pretrained models
#' of the package cannot be changed or removed. 
#' This can be done with the new trained model list.
#' @param path.to.models path to the folder containing the list of new models.
#' 
#' @return no return value, but the model is now saved to database
#' 
#' @examples
#' # load small example dataset
#' data("tirosh_mel80_example")
#' 
#' # train classifier
#' selected_features_T = c("CD4", "CD8A", "CD8B")
#' set.seed(123)
#' clf_t <- train_classifier(train_obj = tirosh_mel80_example, 
#' features = selected_features_T, cell_type = "t cells")
#' 
#' # save the trained classifier to system 
#' # (test classifier can be used before this step)
#' save_new_model(new_model = clf_t, path.to.models = ".")
#' 
#' # verify if new model has been saved
#' print(names(load("./new_models.rda")))
#' delete_model("t cells")
#' 
#' @export
setGeneric("save_new_model", 
           function(new_model, 
                    include.default = TRUE, 
                    path.to.models = ".") 
  standardGeneric("save_new_model"))

#' @inherit save_new_model
#' 
#' @importFrom utils data
#' @rdname save_new_model
setMethod("save_new_model", c("new_model" = "SingleCellClassR"), 
          function(new_model, 
                   include.default = TRUE, 
                   path.to.models = ".") {
  default_models <- NULL
  
  utils::data("default_models")
  new_models.file.path = paste0(path.to.models, "/new_models.rda")
  
  if (file.exists(new_models.file.path)) {
    load(new_models.file.path)
  } else {
    new_models = NULL
  }
  
  if (include.default == TRUE) {
    # default models not in new_models will be added to new_models
    to.be.added <- default_models[!names(default_models)%in%names(new_models)]
    new_models <- append(to.be.added, new_models)
  }
  
  # check if new cell type already existed
  if (new_model@cell_type %in% names(new_models)) {
    stop("A model already existed for the cell type. 
         Please delete the old model to add a new one.", call. = FALSE)
  } 
  
  if (!is.na(new_model@parent) && !new_model@parent %in% names(new_models)) {
    stop("Model for parent cell type has not been added. 
         Please save the model classifying parent cell type into tree first.", 
         call. = FALSE)
  }
  
  # add new model to list
  names <- append(names(new_models), new_model@cell_type)
  new_models <- append(new_models, new_model)
  names(new_models) <- names
  
  # save to rda file
  save(new_models, file = new_models.file.path, compress = 'xz')
  cat("Finished saving new model\n")
})

#' Plant tree from list of models
#' 
#' @param models.file.path list of models. If not provided, 
#' list of default pretrained models in the package will be used.
#'  
#' @return tree structure and plot of tree 
#'
#' @examples
#' 
#' # to create the tree of classifiers 
#' # (in this example, based on default classifiers)
#' plant_tree()
#' 
#' @export
setGeneric("plant_tree", function(models.file.path = c("default", ".")) 
  standardGeneric("plant_tree"))

#' @inherit plant_tree
#' 
#' @importFrom utils data
setMethod("plant_tree", , function(models.file.path = c("default", ".")) {
  new_models <- default_models <- NULL
  
  root.name <- "cell types"
  if ("default" %in% models.file.path) {
    utils::data("default_models", envir = environment())
    model_list <- default_models
  } else {
    models_file_path <- paste0(models.file.path, "/new_models.rda")
    if (!file.exists(models_file_path)) {
      stop("No file exists in the indicated models file path", 
           call. = FALSE)
    } else {
      load(models_file_path, envir = environment())
      model_list <- new_models
    }
  }
  
  tree <- NULL
  
  if (!is.null(model_list)) {
    for (model in model_list) {
      if (is.na(model@parent))
        parent.pathString = root.name
      else 
        parent.pathString <- tree[tree$cell_type == model@parent,]$pathString
      
      cell.info <- c(model@cell_type, model@parent, 
                     paste0(parent.pathString, "/", model@cell_type))
      cell.info <- as.data.frame(matrix(cell.info, nrow = 1))
      colnames(cell.info) <- c("cell_type", "parent_cell_type", "pathString")
      
      if (is.null(tree))
        tree <- cell.info
      else
        tree <- rbind(tree, cell.info)
    }
  }
  
  if (!is.null(tree)) {
    tree <- data.tree::as.Node(tree)
    print(tree)
  } else stop('Tree not available.')
  
  return(tree)
})

#' Delete model/branch from package
#' 
#' @param cell_type string indicating the cell type of which 
#' the model will be removed from package
#' Attention: deletion of a parent model will also delete all of child model.
#' @param path.to.models path to the folder containing 
#' the list of models in which the to-be-deleted model is.
#' 
#' @return no return value, but the model is deleted from database
#' 
#' @examples 
#' # load small example dataset
#' data("tirosh_mel80_example")
#' 
#' # train a classifier
#' set.seed(123)
#' selected_features_T = c("CD4", "CD8A", "CD8B")
#' clf_t <- train_classifier(train_obj = tirosh_mel80_example, 
#' features = selected_features_T, cell_type = "t cells")
#' 
#' # save a classifier to system
#' save_new_model(new_model = clf_t)
#' 
#' # delete classifier from system
#' delete_model("t cells")
#' @export
setGeneric("delete_model", function(cell_type, path.to.models = ".") 
  standardGeneric("delete_model"))

#' @inherit delete_model
setMethod("delete_model", , function(cell_type, path.to.models = ".") {
  new_models <- NULL
  
  new_models.file.path <- paste0(path.to.models, "/new_models.rda")
  if (!file.exists(new_models.file.path)) {
    stop("No list of models available", call. = FALSE)
  } else {
    load(new_models.file.path, envir = environment())
  }
  
  if (is.null(new_models))
    stop("No model to be deleted", call. = FALSE)
  
  to.be.removed <- c(cell_type)
  
  # verify cell type avaibility in tree
  if (!cell_type %in% names(new_models))
    stop("Cannot delete unavailable model/cell type.", call. = FALSE)
  else {
    # get a list of models to delete
    for (model in new_models) {
      if (model@parent %in% to.be.removed) {
        to.be.removed <- append(to.be.removed, model@cell_type)
      }
    }
    
    new_models[to.be.removed] <- NULL
  }
  
  # save models after remove
  if (!is.null(new_models))
    save(new_models, file = new_models.file.path)
})
