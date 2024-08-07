##' Ornstein-Uhlenbeck methods for comparative phylogenetic hypotheses
##'
##' The \pkg{ouch} package provides facilities for phylogenetic comparative analysis based on Ornstein-Uhlenbeck models of trait evolution along a phylogeny.
##' Multivariate data and complex adaptive hypotheses are supported.
##'
##' @name ouch-package
##' @aliases ouch,package ouch-package
##' @rdname package
##' @family phylogenetic comparative models
##' @family methods for ouch trees
##' @family examples
##' @section Classes:
##' The basic class, `ouchtree`, is provided to encode a phylogenetic tree.
##' Plot and print methods are provided.
##'
##' The class `browntree` derives from class `ouchtree` and encodes the results of fitting a Brownian Motion model to data.
##'
##' The class `hansentree` derives from class `ouchtree` and encodes the results of fitting a Hansen model to data.
##' @section Detailed Documentation:
##'   - Phylogenies in \pkg{ouch} format: [`ouchtree()`], [`ape2ouch()`]
##'   - Brownian motion models: [`brown()`]
##'   - Ornstein-Uhlenbeck models: [`hansen()`], [`paint()`]
##'   - Simulation of models: [`simulate()`]
##'   - Display of data: [`plot()`]
##'   - Extraction of information from fitted models: [`summary()`], [`logLik()`], [`coef()`]
##'   - Example datasets: [`anolis.ssd`], [`bimac`]
##' @section Citing \pkg{ouch}:
##' Execute \code{citation("ouch")} to view the correct citation for publications.
##' @author Aaron A. King
##' @references
##' \Hansen1997
##'
##' \Butler2004
##'
##' \Cressler2015
##'
##' @useDynLib ouch, .registration = TRUE
##' @import methods
##' @keywords models
"_PACKAGE"

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

##' @importFrom stats runif
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
