#' Ornstein-Uhlenbeck methods for comparative phylogenetic hypotheses
#' 
#' The \pkg{ouch} package provides facilities for phylogenetic comparative analysis based on Ornstein-Uhlenbeck models of trait evolution along a phylogeny.
#' Multivariate data and complex adaptive hypotheses are supported.
#' 
#' @name ouch-package
#' @rdname package
#' @docType package
#' @section Classes:
#' The basic class, \code{ouchtree}, is provided to encode a
#' phylogenetic tree.  Plot and print methods are provided.
#' 
#' The class \code{browntree} derives from class \code{ouchtree} and encodes
#' the results of fitting a Brownian Motion model to data.
#' 
#' The class \code{hansentree} derives from class \code{ouchtree} and encodes
#' the results of fitting a Hansen model to data.
#' @section Detailed Documentation:
#' \describe{
#'   \item{Phylogenies in \pkg{ouch} format}{\code{\link{ouchtree}}, \code{\link{ape2ouch}}}
#'   \item{Brownian motion models}{\code{\link{brown}}}
#'   \item{Ornstein-Uhlenbeck models}{\code{\link{hansen}}}
#'   \item{Simulation of models}{\code{\link[ouch:simulate]{simulate}}}
#'   \item{Display of data}{\code{\link[ouch:plot]{plot}}}
#'   \item{Examples}{\code{\link{anolis.ssd}}, \code{\link{bimac}}}
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

undef_method <- function (method, object) {
  o <- deparse(substitute(object))
  stop(sQuote(method)," is undefined for ",sQuote(o)," of class ",
    sQuote(class(object)),".",call.=FALSE)
}

reqd_arg <- function (method, object) {
  if (is.null(method) || length(method)==0)
    stop(sQuote(object)," is a required argument.",call.=FALSE)
  else
    stop(method,sQuote(object)," is a required argument.",call.=FALSE)
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
