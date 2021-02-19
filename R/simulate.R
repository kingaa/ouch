#' Simulations of a phylogenetic trait model.
#'
#' \code{simulate} generates random deviates from a fitted model.
#' 
#' @name simulate
#' @docType methods
#' @rdname simulate
#' @family methods
#' @importFrom stats simulate
#'
#' @return
#' \code{simulate} returns a list of data-frames, each comparable to the original data.
#'
#' @param object fitted model object to simulate.
#' @param nsim integer; number of independent simulations.
#' @param seed integer; if non-\code{NULL}, the RNG will be initialized with this seed for the simulations.
#' The RNG will be reset to its pre-existing state when \code{simulate} returns.
#' @param ... additional arguments, ignored.
NULL

