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
#' @param path_to_models path to the folder containing the list of new models.
#' 
#' @return no return value, but the model is now saved to database
#' 
#' @importFrom utils data
#' 
#' @examples
#' # load small example dataset
#' data("tirosh_mel80_example")
#' 
#' # train classifier
#' selected_marker_genes_T = c("CD4", "CD8A", "CD8B")
#' set.seed(123)
#' classifier_t <- train_classifier(train_obj = tirosh_mel80_example,
#' assay = 'RNA', slot = 'counts', marker_genes = selected_marker_genes_T, 
#' cell_type = "t cells", tag_slot = 'active.ident')
#' 
#' # save the trained classifier to system 
#' # test classifier can be used before this step
#' # Note: We do not include the default models here to runtime of the example
#' save_new_model(new_model = classifier_t, path_to_models = tempdir(), 
#'                include.default = FALSE)
#' 
#' # verify if new model has been saved
#' print(names(load(file.path(tempdir(), "new_models.rda"))))
#' delete_model("t cells")
#' 
#' @export
save_new_model <- function(new_model, include.default = TRUE, 
                    path_to_models = tempdir()) {
  # create the model directory if it doesn't exist
  if (!dir.exists(path_to_models)) {
    dir.create(path_to_models)
  }
  
  default_models <- NULL
  data_env <- new.env(parent = emptyenv())
  
  new_models.file.path = file.path(path_to_models, "new_models.rda")
  
  if (file.exists(new_models.file.path)) {
    load(new_models.file.path, data_env)
    new_models <- data_env[["new_models"]]
  } else {
    new_models = NULL
  }
  
  if (include.default == TRUE) {
    # default models not in new_models will be added to new_models
    default_models <- download_data_file()
    to.be.added <- default_models[!names(default_models)%in%names(new_models)]
    new_models <- append(to.be.added, new_models)
  }
  
  # check if new cell type already existed
  if (cell_type(new_model) %in% names(new_models)) {
    stop("A model already existed for the cell type. 
         Please delete the old model to add a new one.", call. = FALSE)
  } 
  
  if (!is.na(parent(new_model)) && !parent(new_model) %in% names(new_models)) {
    stop("Model for parent cell type has not been added. 
         Please save the model classifying parent cell type into tree first.", 
         call. = FALSE)
  }
  
  # add new model to list
  names <- append(names(new_models), cell_type(new_model))
  new_models <- append(new_models, new_model)
  names(new_models) <- names
  
  # save to rda file
  message("Saving new models to ", new_models.file.path, "...")
  save(new_models, file = new_models.file.path, compress = 'xz')
  message("Finished saving new model")
} 

#' Plant tree from list of models
#' 
#' @param path_to_models list of models. If not provided, 
#' list of default pretrained models in the package will be used.
#'  
#' @return tree structure and plot of tree 
#'
#' @importFrom utils data
#' @import data.tree
#'
#' @examples
#' 
#' # to create the tree of classifiers 
#' # (in this example, based on default classifiers)
#' plant_tree()
#' 
#' @export
plant_tree <- function(path_to_models = "default") { 
  data_env <- new.env(parent = emptyenv())
  
  root.name <- "cell types"
  
  model_list <- load_models(path_to_models) 
  tree <- NULL
  
  if (!is.null(model_list)) {
    for (model in model_list) {
      if (is.na(parent(model)))
        parent.pathString = root.name
      else 
        parent.pathString <- tree[tree$cell_type == parent(model),]$pathString
      
      cell.info <- c(
        cell_type(model), parent(model), 
        paste0(parent.pathString, "/", cell_type(model))
      )
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
}

#' Delete model/branch from package
#' 
#' @param cell_type string indicating the cell type of which 
#' the model will be removed from package
#' Attention: deletion of a parent model will also delete all of child model.
#' @param path_to_models path to the folder containing 
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
#' selected_marker_genes_T = c("CD4", "CD8A", "CD8B")
#' classifier_t <- train_classifier(train_obj = tirosh_mel80_example,
#' assay = 'RNA', slot = 'counts', marker_genes = selected_marker_genes_T, 
#' cell_type = "t cells", tag_slot = 'active.ident')
#' 
#' # save a classifier to system
#' save_new_model(new_model = classifier_t, path_to_models = tempdir(), 
#'                include.default = FALSE)
#' 
#' # delete classifier from system
#' delete_model("t cells", path_to_models = tempdir())
#' @export
delete_model <- function(cell_type, path_to_models = tempdir()) {
  new_models <- NULL
  data_env <- new.env(parent = emptyenv())
  
  if (path_to_models == 'default')
    stop("Cannot delete default models.", call. = FALSE)
  
  new_models.file.path <- file.path(path_to_models, "new_models.rda")
  if (!file.exists(new_models.file.path)) {
    stop("No list of models available", call. = FALSE)
  } else {
    load(new_models.file.path, envir = data_env)
    new_models <- data_env[['new_models']]
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
      if (parent(model) %in% to.be.removed) {
        to.be.removed <- append(to.be.removed, cell_type(model))
      }
    }
    
    new_models[to.be.removed] <- NULL
  }
  
  # save models after remove
  if (!is.null(new_models))
    save(new_models, file = new_models.file.path)
}
