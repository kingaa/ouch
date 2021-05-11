#' Ornstein-Uhlenbeck methods for comparative phylogenetic hypotheses
#' 
#' The \pkg{ouch} package provides facilities for phylogenetic comparative analysis based on Ornstein-Uhlenbeck models of trait evolution along a phylogeny.
#' Multivariate data and complex adaptive hypotheses are supported.
#' 
#' @name ouch-package
#' @aliases ouch
#' @rdname package
#' @docType package
#' @section Classes:
#' The basic class, `ouchtree`, is provided to encode a phylogenetic tree.
#' Plot and print methods are provided.
#' 
#' The class `browntree` derives from class `ouchtree` and encodes the results of fitting a Brownian Motion model to data.
#' 
#' The class `hansentree` derives from class `ouchtree` and encodes the results of fitting a Hansen model to data.
#' @section Detailed Documentation:
#' \describe{
#'   \item{Phylogenies in \pkg{ouch} format}{[`ouchtree()`], [`ape2ouch()`]}
#'   \item{Brownian motion models}{[`brown()`]}
#'   \item{Ornstein-Uhlenbeck models}{[`hansen()`]}
#'   \item{Simulation of models}{[`simulate()`]}
#'   \item{Display of data}{[`plot()`]}
#'   \item{Examples}{[`anolis.ssd`], [`bimac`]}
#' }
#' @author Aaron A. King
#' @references
#' \Butler2004
#'
#' \Cressler2015
#' 
#' @useDynLib ouch, .registration = TRUE
#' @import methods
#' @keywords models
NULL

pStop <- function (fn, ...) {
  fn <- as.character(fn)
  stop("in ",sQuote(fn[1L]),": ",...,call.=FALSE)
}

pStop_ <- function (...) {
  stop(...,call.=FALSE)
}

pWarn <- function (fn, ...) {
  fn <- as.character(fn)
  warning("in ",sQuote(fn[1L]),": ",...,call.=FALSE)
}

pWarn_ <- function (...) {
  warning(...,call.=FALSE)
}

undef_method <- function (method, object) {
  o <- deparse(substitute(object))
  pStop_(sQuote(method)," is undefined for ",sQuote(o)," of class ",
    sQuote(class(object)),".")
}

reqd_arg <- function (method, object) {
  if (is.null(method) || length(method)==0)
    pStop_(sQuote(object)," is a required argument.")
  else
    pStop(method,sQuote(object)," is a required argument.")
}

undef_method <- function (method, object) {
  o <- deparse(substitute(object))
  stop(sQuote(method)," is undefined for ",sQuote(o)," of class ",
    sQuote(class(object)),".",call.=FALSE)
}

#' @importFrom stats runif
freeze <- function (seed = NULL) {
  if (!is.null(seed)) {
    if (!exists('.Random.seed',envir=.GlobalEnv)) runif(1)
    save.seed <- get('.Random.seed',envir=.GlobalEnv)
    set.seed(seed)
    invisible(save.seed)
  } else {
    invisible(NULL)
  }
}

thaw <- function (seed = NULL) {
  if (!is.null(seed)) {
    assign('.Random.seed',seed,envir=.GlobalEnv)
  }
  invisible(NULL)
}
