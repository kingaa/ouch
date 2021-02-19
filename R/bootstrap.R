#' Bootstrapping for uncertainty quantification
#'
#' Generic bootstrapping for \pkg{ouch} models.
#'
#' \code{bootstrap} performs a parametric bootstrap for estimation of confidence intervals.
#'
#' @rdname bootstrap
#' @name bootstrap
#' @family methods
#'
#' @param object A fitted model object.
#' @param nboot integer; number of bootstrap replicates.
#' @param seed integer; setting \code{seed} to a non-\code{NULL} value allows one to fix the random seed (see \code{\link[ouch:simulate]{simulate}}).
#' @param ... Additional arguments are passed to \code{\link[ouch:update]{update}}.
NULL

setGeneric(
  "bootstrap",
  function (object, ...) {
    standardGeneric("bootstrap")
  }
)

#' @rdname bootstrap
#' @export
setMethod(
  "bootstrap",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("bootstrap","object")
  }
)

#' @rdname bootstrap
#' @export
setMethod(
  "bootstrap",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("bootstrap",object)
  }
)
