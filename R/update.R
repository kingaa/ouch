#' Update and refit a model.
#'
#' \code{update} will update a model and re-fit.
#' This allows one to change the data and/or parameters.
#' 
#' @name update
#' @docType methods
#' @rdname update
#' @family methods
#' @importFrom stats update
#'
#' @return
#' @return \code{update} returns a new fitted-model object of the same class as  \code{object}.
#' @param object fitted model object.
#' @param data data that replace those used in the original fit.
#' @param ... Additional arguments replace the corresponding arguments in the original call.
NULL

