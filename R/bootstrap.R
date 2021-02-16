#' Bootstrapping for uncertainty quantification
#'
#' Generic bootstrapping for \pkg{ouch} models.
#'
#' \code{bootstrap} performs a parametric bootstrap for estimation of confidence intervals.
#'
#' @rdname bootstrap
#' @name bootstrap
#'
#' @include brown.R hansen.R
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
  signature=signature(object="browntree"),
  function (object, nboot = 200, seed = NULL, ...) {
    simdata <- simulate(object,nsim=nboot,seed=seed)
    results <- vector(mode='list',length=nboot)
    toshow <- c("sigma.squared","theta","loglik","aic","aic.c","sic","dof")
    for (b in seq_len(nboot)) {
      results[[b]] <- summary(update(object,data=simdata[[b]],...))
    }
    as.data.frame(t(sapply(results,function(x)unlist(x[toshow]))))
  }
)

#' @rdname bootstrap
#' @export
setMethod(
  "bootstrap",
  signature=signature(object="hansentree"),
  function (object, nboot = 200, seed = NULL, ...) {
    simdata <- simulate(object,nsim=nboot,seed=seed)
    results <- vector(mode='list',length=nboot)
    toshow <- c("alpha","sigma.squared","optima","loglik","aic","aic.c","sic","dof")
    for (b in seq_len(nboot)) {
      results[[b]] <- summary(update(object,data=simdata[[b]],...))
    }
    as.data.frame(t(sapply(results,function(x)unlist(x[toshow]))))
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
