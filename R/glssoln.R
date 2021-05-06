#' Generalized least-squares solver
#'
#' Solves the generalized least squares problem.
#'
#' Given matrices \eqn{a}, \eqn{x}, \eqn{v}, `glssoln` computes \eqn{y} such that
#' \deqn{(x-ay)^T v^{-1} (x-ay)}
#' is minimized.
#' This is accomplished by first computing the Choleski decomposition of \eqn{v}:
#' \deqn{v=LL^T.}
#' One then solves for \eqn{y} in the equation
#' \deqn{L^{-1}ay=L^{-1}x.}
#' This is accomplished by means of a singular-value decomposition of \eqn{L^{-1} a}.
#'
#' The resulting \eqn{y} then satisfies
#' \deqn{x=ay+e,}
#' where the entries of \eqn{e} are the residuals.
#'
#' @keywords internal
#' @name glssoln
#' @rdname glssoln
#' @return
#' `glssoln` returns a list of two named components:
#' - `coeff` is \eqn{y} as above.
#' - `residuals` is \eqn{e} as above.
glssoln <- function (a, x, v, tol = sqrt(.Machine$double.eps)) {
  n <- length(x)
  vh <- tryCatch(
    chol(v),
    error=function (e) {
      pWarn(
        "ouch:::glssoln","Choleski decomposition ",
        "of variance-covariance matrix fails with error: ",
        dQuote(conditionMessage(e)),"."
      )
      NULL
    }
  )
  if (is.null(vh)) {
    y <- rep(NA,NCOL(a))
    e <- rep(NA,n)
    dim(y) <- NCOL(a)
    dim(e) <- n
  } else {
    s <- svd(forwardsolve(vh,a,upper.tri=TRUE,transpose=TRUE))
    inds <- s$d > tol*max(s$d)
    svals <- s$d[inds,drop=FALSE]
    r <- length(svals)
    svals <-  diag(1/svals,nrow=r,ncol=r)
    y <- (s$v[,inds,drop=FALSE]%*%
            (svals %*%
               t(s$u[,inds,drop=FALSE])))%*%
      forwardsolve(vh,x,upper.tri=TRUE,transpose=TRUE)
    e <- a%*%y-x
    dim(y) <- dim(y)[1]
    dim(e) <- n
  }
  list(
    coeff=y,
    residuals=e
  )
}
