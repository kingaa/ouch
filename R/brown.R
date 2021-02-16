#' Phylogenetic Brownian motion models
#' 
#' The function \code{brown} creates a \code{browntree} object by fitting a
#' Brownian-motion model to data.
#' 
#' @name brown
#' @rdname brown
#' @author Aaron A. King
#' @seealso \code{\link{ouchtree}}, \code{\link{hansen}}, \code{\link{bimac}}, \code{\link{anolis.ssd}}
#' @keywords models
#' @references
#' \Butler2004
#' @keywords models
#'
NULL

setClass(
  'browntree',
  contains='ouchtree',
  representation=representation(
    call='call',
    nchar='integer',
    data='list',
    theta='list',
    sigma='numeric',
    loglik='numeric'
  )
)

#' @rdname brown
#' @include ouchtree.R glssoln.R rmvnorm.R
#' 
#' @param data Phenotypic data for extant species, i.e., at the terminal ends
#' of the phylogenetic tree.  This can either be a numeric vector or a list.
#' If it is a numeric vector, there must be one entry for every node.  If it is
#' a list, it must consist entirely of numeric vectors, each of which has one
#' entry per node.  A data-frame is coerced to a list.
#' @param tree A phylogenetic tree, specified as an \code{ouchtree} object.
#'
#' @return \code{brown} returns an object of class \code{browntree}.
#' 
#' @export
brown <- function (data, tree) {
  
  if (!is(tree,'ouchtree'))
    stop(sQuote("tree")," must be an object of class ",sQuote("ouchtree"))
  
  if (is.data.frame(data)) {
    nm <- rownames(data)
    data <- lapply(as.list(data),function(x){names(x)<-nm;x})
  }
  if (is.numeric(data)) {
    nm <- deparse(substitute(data))[1]
    data <- list(data)
    names(data) <- nm
  }
  if (is.list(data)) {
    if (
      any(sapply(data,class)!='numeric') ||
        any(sapply(data,length)!=tree@nnodes)
    )
      stop(sQuote("data")," vector(s) must be numeric, with one entry per node of the tree")
    if (any(sapply(data,function(x)(is.null(names(x)))||(!setequal(names(x),tree@nodes)))))
      stop(sQuote("data"), " vector names (or data-frame row names) must match node names of ", sQuote("tree"))
    for (xx in data) {
      no.dats <- which(is.na(xx[tree@nodes[tree@term]]))
      if (length(no.dats)>0)
        stop("missing data on terminal node(s): ",paste(tree@nodes[tree@term[no.dats]],collapse=','))
    }
  } else
    stop(sQuote("data")," must be either a single numeric data set or a list of numeric data sets")

  nterm <- tree@nterm
  nchar <- length(data)
  if (is.null(names(data))) names(data) <- paste('char',1:nchar,sep='')

  data <- lapply(data,function(x)x[tree@nodes])
  
  dat <- lapply(data,function(y)y[tree@term])

  w <- matrix(data=1,nrow=nterm,ncol=1)
  b <- tree@branch.times
  
  sols <- lapply(dat,function(x)glssoln(w,x,b))
  theta <- lapply(sols,function(x)x$coeff) # optima
  e <- sapply(sols,function(x)x$residuals) # residuals
  q <- t(e)%*%solve(b,e)
  v <- q/nterm
  dev <- nchar*nterm*(1+log(2*pi))+nchar*log(det(b))+nterm*log(det(v))
  sigma <- t(chol(v))
  
  new(
    'browntree',
    as(tree,'ouchtree'),
    call=match.call(),
    nchar=nchar,
    data=as.list(data),
    theta=theta,
    sigma=sigma[lower.tri(sigma,diag=T)],
    loglik=-0.5*dev
  )
}

brown.deviate <- function(n = 1, object) {
  sigma <- matrix(0,object@nchar,object@nchar)
  sigma[lower.tri(sigma,diag=T)] <- object@sigma
  v <- kronecker(sigma%*%t(sigma),object@branch.times)
  m <- sapply(object@theta,rep,object@nterm)
  X <- array(
    data=NA,
    dim=c(object@nnodes,object@nchar,n),
    dimnames=list(
      object@nodes,
      names(object@data),
      paste('rep',seq(n),sep='.')
    )
  )
  X[object@term,,] <- array(
    rmvnorm(n=n,mean=m,var=v),
    dim=c(object@nterm,object@nchar,n)
  )
  apply(X,3,as.data.frame)
}

#' @rdname simulate
#' @include simulate.R
#' @importFrom stats runif
#' @param object fitted model object
#' @export
setMethod(
  'simulate',
  signature=signature(object='browntree'),
  function (object, nsim = 1, seed = NULL, ...) {
    if (!is.null(seed)) {
      if (!exists('.Random.seed',envir=.GlobalEnv)) runif(1)
      save.seed <- get('.Random.seed',envir=.GlobalEnv)
      set.seed(seed)
    }
    X <- brown.deviate(n=nsim,object)
    if (!is.null(seed)) {
      assign('.Random.seed',save.seed,envir=.GlobalEnv)
    }
    X
  }
)

#' @rdname update
#' @importFrom stats update
#' @export
setMethod(
  'update',
  signature=signature(object='browntree'),
  function (object, data, ...) {
    if (missing(data)) data <- object@data
    brown(
      data=data,
      tree=object,
      ...
    )
  }
)

#' @rdname coef
#' @importFrom stats coef
#' @return \code{coef} applied to a \code{browntree} object extracts a list with three elements:
#' \describe{
#'   \item{\code{sigma}}{the coefficients of the sigma matrix.}
#'   \item{\code{theta}}{a list of the estimated optima, one per character.}
#'   \item{\code{sigma..sq.matrix}}{the sigma-squared matrix itself.}
#' }
#' @export
setMethod(
  'coef',
  signature=signature(object="browntree"),
  function (object, ...) {
    list(
      sigma=object@sigma,
      theta=object@theta,
      sigma.sq.matrix=sym.par(object@sigma)
    )
  }
)

#' @rdname print
#' @export
setMethod(
  'show',
  signature=signature(object="browntree"),
  function (object) {
    print(as(object,'browntree'))
    invisible(NULL)
  }
)

#' @rdname print
#' @return \code{print} displays the tree as a table, along with the coefficients of the fitted model and diagnostic information.
#' @export
setMethod(
  'print',
  signature=signature(x='browntree'),
  function (x, ...) {
    cat("\ncall:\n")
    print(x@call)
    print(as(x,'data.frame'),...)
    sm <- summary(x)
    cat('\nsigma squared:\n')
    print(sm$sigma.squared)
    cat('\ntheta:\n')
    print(sm$optima)
    print(unlist(sm[c("loglik","deviance","aic","aic.c","sic","dof")]))
    invisible(x)
  }
)

#' @rdname loglik
#' @importFrom stats logLik
#' @export
setMethod(
  "logLik",
  signature=signature(object="browntree"),
  function(object)object@loglik
)

#' @rdname summary
#' @return \code{summary} applied to a \code{browntree} object returns information about the fitted model, including parameter estimates and quantities describing the goodness of fit.
#' @export
setMethod(
  "summary",
  signature=signature(object="browntree"),
  function (object, ...) {
    cf <- coef(object)
    dof <- object@nchar+length(object@sigma)
    deviance=-2*logLik(object)
    aic <- deviance+2*dof
    aic.c <- aic+2*dof*(dof+1)/(object@nterm*object@nchar-dof-1)
    sic <- deviance+log(object@nterm*object@nchar)*dof
    list(
      call=object@call,
      sigma.squared=cf$sigma.sq.matrix,
      theta=cf$theta,
      loglik=logLik(object),
      deviance=deviance,
      aic=aic,
      aic.c=aic.c,
      sic=sic,
      dof=dof
    )
  }
)

setAs(
  from='browntree',
  to='data.frame',
  def = function (from) {
    cbind(
      as(as(from,'ouchtree'),'data.frame'),
      as.data.frame(from@data)
    )
  }
)
