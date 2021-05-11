#' Bootstrapping for uncertainty quantification
#'
#' Generic bootstrapping for \pkg{ouch} models.
#'
#' `bootstrap` performs a parametric bootstrap for estimation of confidence intervals.
#'
#' @rdname bootstrap
#' @name bootstrap
#' @family methods for \pkg{ouch} trees
#' @example examples/bootstrap.R
#' @param object A fitted model object.
#' @param nboot integer; number of bootstrap replicates.
#' @param seed integer; setting `seed` to a non-`NULL` value allows one to fix the random seed (see [simulate]).
#' @param ... Additional arguments are passed to [`update`].
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
