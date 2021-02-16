glssoln <- function (a, x, v, tol = sqrt(.Machine$double.eps)) {
  n <- length(x)
  vh <- try(
            chol(v),
            silent=FALSE
            )
  if (inherits(vh,'try-error')) {
    warning(
            "glssoln: Choleski decomposition of variance-covariance matrix fails",
            call.=FALSE
            )
    y <- rep(NA,ncol(a))
    e <- rep(NA,n)
    dim(y) <- ncol(a)
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
