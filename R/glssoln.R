glssoln <- function (a, x, v, tol = sqrt(.Machine$double.eps)) {
  n <- length(x)
  vh <- tryCatch(
    chol(v),
    error=function (e) {
      warning(
        "in ",sQuote("ouch:::glssoln"),": Choleski decomposition ",
        "of variance-covariance matrix fails with error: ",
        dQuote(conditionMessage(e)),
        call.=FALSE
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
