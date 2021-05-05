#' Ornstein-Uhlenbeck models of trait evolution
#' 
#' The function \code{hansen} fits an Ornstein-Uhlenbeck model to data.
#' The fitting is done using \code{optim} or \code{subplex}.
#' 
#' The Hansen model for the evolution of a multivariate trait \eqn{X} along a lineage can be written as a stochastic differential equation (Ito diffusion)
#' \deqn{dX=\alpha(\theta(t)-X(t))dt+\sigma dB(t),}{dX = alpha (theta(t)-X(t)) dt + sigma dB(t),}
#' where \eqn{t} is time along the lineage,
#' \eqn{\theta(t)}{theta(t)} is the optimum trait value, \eqn{B(t)} is a standard Wiener process (Brownian motion),
#' and \eqn{\alpha}{alpha} and \eqn{\sigma}{sigma} are matrices
#' quantifying, respectively, the strength of selection and random drift.
#' Without loss of generality, one can assume \eqn{\sigma}{sigma} is lower-triangular.
#' This is because only the infinitesimal variance-covariance matrix
#' \eqn{\sigma^2=\sigma\sigma^T}{sigma^2 = sigma\%*\%transpose(sigma)}
#' is identifiable, and for any admissible variance-covariance matrix, we can choose \eqn{\sigma}{sigma} to be lower-triangular.
#' Moreover, if we view the basic model as describing evolution on a fitness landscape, then \eqn{\alpha}{alpha} will be symmetric.
#' If we further restrict ourselves to the case of stabilizing selection, \eqn{\alpha}{alpha} will be positive definite as well.
#' We make these assumptions and therefore can assume that the matrix \eqn{\alpha}{alpha} has a lower-triangular square root.
#' 
#' The \code{hansen} code uses unconstrained numerical optimization to maximize the likelihood.
#' To do this, it parameterizes the \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2} matrices in a special way:
#' each matrix is parameterized by \code{nchar*(nchar+1)/2} parameters, where \code{nchar} is the number of quantitative characters.
#' Specifically, the parameters initialized by the \code{sqrt.alpha} argument of \code{hansen} are used
#' to fill the nonzero entries of a lower-triangular matrix (in column-major order),
#' which is then multiplied by its transpose to give the selection-strength matrix.
#' The parameters specified in \code{sigma} fill the nonzero entries in the lower triangular \eqn{\sigma}{sigma} matrix.
#' When \code{hansen} is executed, the numerical optimizer maximizes the likelihood over these parameters.
#'
#' @name hansen
#' @aliases hansentree-class
#' @rdname hansen
#' @family phylogenetic comparative models
#' @author Aaron A. King
#' @seealso
#' \code{\link[stats:optim]{optim}}, \code{\link[subplex:subplex]{subplex}}, \code{\link{bimac}}, \code{\link{anolis.ssd}}
#' @references
#' \Butler2004
#'
#' \Cressler2015
#' 
#' @keywords models
#' @example examples/hansen.R
#' 
NULL

setClass(
  'hansentree',
  contains='ouchtree',
  representation=representation(
    call='call',
    nchar='integer',
    optim.diagn='list',
    hessian='matrix',
    data='list',
    regimes='list',
    beta='list',
    theta='list',
    sigma='numeric',
    sqrt.alpha='numeric',
    loglik='numeric'
  )
)

setAs(
  from='hansentree',
  to='data.frame',
  def = function (from) {
    cbind(
      as(as(from,'ouchtree'),'data.frame'),
      as.data.frame(from@regimes),
      as.data.frame(from@data)
    )
  }
)

#' @rdname hansen
#' @include ouchtree.R glssoln.R rmvnorm.R
#' @importFrom stats optim
#' @importFrom subplex subplex
#'
#' @param data Phenotypic data for extant species, i.e., species at the terminal twigs of the phylogenetic tree.
#' This can either be a single named numeric vector, a list of \code{nchar} named vectors, or a data frame containing \code{nchar} data variables.
#' There must be an entry per variable for every node in the tree; use \code{NA} to represent missing data.
#' If the
#' data are supplied as one or more named vectors, the names attributes are taken to correspond to the node names specified when the \code{ouchtree} was constructed (see \code{\link{ouchtree}}).
#' If the data are supplied as a
#' data-frame, the rownames serve that purpose.
#' @param tree A phylogenetic tree, specified as an \code{ouchtree} object.
#' @param regimes A vector of codes, one for each node in the tree, specifying the selective regimes hypothesized to have been operative.
#' Corresponding to each node, enter the code of the regime hypothesized for the branch segment terminating in that node.
#' For the root node, because it has no branch segment terminating on it, the regime specification is irrelevant.
#' If there are \code{nchar} quantitative characters, then one can specify a single set of \code{regimes} for all characters or a list of \code{nchar} regime specifications, one for each character.
#' @param sqrt.alpha,sigma These are used to initialize the optimization algorithm.
#' The selection strength matrix \eqn{\alpha}{alpha} and the random drift variance-covariance matrix \eqn{\sigma^2}{sigma^2} are parameterized by their matrix square roots.
#' Specifically, these initial guesses are each packed into lower-triangular matrices (column by column).
#' The product of this matrix with its transpose is the \eqn{\alpha}{alpha} or \eqn{\sigma^2}{sigma^2} matrix.
#' See Details for more information.
#' @param fit If \code{fit=TRUE}, then the likelihood will be maximized.
#' If \code{fit=FALSE}, the likelihood will be evaluated at the specified values of \code{sqrt.alpha} and \code{sigma};
#' the optima \code{theta} will be returned as well.
#' @param method The method to be used by the optimization algorithm.
#' See \code{\link[subplex]{subplex}} and \code{\link{optim}} for information on the available options.
#' @param hessian If \code{hessian=TRUE}, then the Hessian matrix will be computed by \code{optim}.
#' @param \dots Additional arguments will be passed as \code{control} options to \code{optim} or \code{subplex}.
#' See \code{\link{optim}} and \code{\link[subplex:subplex]{subplex}} for information on the available options.
#' 
#' @return \code{hansen} returns an object of class \code{hansentree}.
#' @export
hansen <- function (data, tree, regimes, sqrt.alpha, sigma,
  fit = TRUE,
  method = c("Nelder-Mead","subplex","BFGS","L-BFGS-B"),
  hessian = FALSE,
  ...) {

  if (!is(tree,'ouchtree'))
    pStop("hansen",sQuote("tree")," must be an object of class ",sQuote("ouchtree"),".")

  if (missing(data)) {
    if (is(tree,"hansentree")) {
      data <- tree@data
    } else {
      pStop("hansen",sQuote("data")," must be specified.")
    }
  }
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
      pStop("hansen",sQuote("data")," vector(s) must be numeric, with one entry per node of the tree.")
    if (any(sapply(data,function(x)(is.null(names(x)))||(!setequal(names(x),tree@nodes)))))
      pStop("hansen",sQuote("data"), " vector names (or data-frame row names) must match node names of ", sQuote("tree"),".")
    for (xx in data) {
      no.dats <- which(is.na(xx[tree@nodes[tree@term]]))
      if (length(no.dats)>0)
        pStop("missing data on terminal node(s): ",
          paste(sQuote(tree@nodes[tree@term[no.dats]]),collapse=', '),".")
    }
  } else
    pStop("hansen",sQuote("data")," must be either a single numeric data set or a list of numeric data sets.")

  nchar <- length(data)
  if (is.null(names(data))) names(data) <- paste('char',seq_len(nchar),sep='')
  
  if (any(sapply(data,function(x)(is.null(names(x)))||(!setequal(names(x),tree@nodes)))))
    pStop("hansen","each data set must have names corresponding to the node names.")
  data <- lapply(data,function(x)x[tree@nodes])
  dat <- do.call(c,lapply(data,function(y)y[tree@term]))
  
  nsymargs <- nchar*(nchar+1)/2
  nalpha <- length(sqrt.alpha)
  nsigma <- length(sigma)

  if (nalpha!=nsymargs)
    pStop("hansen","the length of ",sQuote("sqrt.alpha")," must be a triangular number.")

  if (nsigma!=nsymargs)
    pStop("hansen","the length of ",sQuote("sigma")," must be a triangular number.")

  if (missing(regimes)) {
    if (is(tree,"hansentree")) {
      regimes <- tree@regimes
      beta <- tree@beta
    } else {
      pStop("hansen",sQuote("regimes")," must be specified.")
    }
  }
  if (is.data.frame(regimes)) {
    nm <- rownames(regimes)
    regimes <- lapply(as.list(regimes),function(x){names(x)<-nm;x})
  }
  if (is.list(regimes)) {
    if (any(sapply(regimes,length)!=tree@nnodes))
      pStop("hansen","each element in ",sQuote("regimes")," must be a vector with one entry per node of the tree.")
  } else {
    if (length(regimes)!=tree@nnodes)
      pStop("hansen","there must be one entry in ",sQuote("regimes")," per node of the tree.")
    nm <- deparse(substitute(regimes))[1]
    regimes <- list(regimes)
    names(regimes) <- nm
  }
  
  if (any(!sapply(regimes,is.factor)))
    pStop("hansen",sQuote("regimes")," must be of class ",sQuote("factor")," or a list of ",sQuote("factor")," objects.")

  if (length(regimes)==1)
    regimes <- rep(regimes,nchar)

  if (length(regimes) != nchar)
    pStop("hansen","you must supply a regime-specification vector for each character.")

  if (any(sapply(regimes,function(x)(is.null(names(x)))||(!setequal(names(x),tree@nodes)))))
    pStop("hansen","each regime specification must have names corresponding to the node names.")
  regimes <- lapply(regimes,function(x)x[tree@nodes])
  beta <- regime_spec(tree,regimes)

  optim.diagn <- vector(mode='list',length=0)

  if (fit) { ## maximize the likelihood

    method <- match.arg(method)

    if (method=='subplex') {
      opt <- subplex(
        par=c(sqrt.alpha,sigma),
        fn = function (par) {
          ou_lik_fn(
            tree=tree,
            alpha=sym_par(par[seq(nalpha)]),
            sigma=sym_par(par[nalpha+seq(nsigma)]),
            beta=beta,
            dat=dat
          )$deviance
        },
        hessian=hessian,
        control=list(...)
      )
      if (opt$convergence!=0) {
        message("unsuccessful convergence, code ",opt$convergence,", see documentation for ",sQuote("subplex"))
        if (!is.null(opt$message))
          message(sQuote("subplex")," message: ",opt$message)
        pWarn("hansen","unsuccessful convergence.")
      }
    } else {
      opt <- optim(
        par=c(sqrt.alpha,sigma),
        fn = function (par) {
          ou_lik_fn(
            tree=tree,
            alpha=sym_par(par[seq(nalpha)]),
            sigma=sym_par(par[nalpha+seq(nsigma)]),
            beta=beta,
            dat=dat
          )$deviance
        },
        gr=NULL,
        hessian=hessian,
        method=method,
        control=list(...)
      )
      if (opt$convergence!=0) {
        message("unsuccessful convergence, code ",opt$convergence,", see documentation for ",sQuote("optim"))
        if (!is.null(opt$message))
          message(sQuote("optim")," message: ",opt$message)
        pWarn("hansen","unsuccessful convergence.")
      }
    }

    sqrt.alpha <- opt$par[seq(nalpha)]
    sigma <- opt$par[nalpha+seq(nsigma)]
    optim.diagn <- list(convergence=opt$convergence,message=opt$message)

  }

  if (hessian) {
    hs <- opt$hessian
    ##     se <- sqrt(diag(solve(0.5*hs)))
    ##     se.alpha <- se[seq(nalpha)]
    ##     se.sigma <- se[nalpha+seq(nsigma)]
  } else {
    hs <- matrix(NA,0,0)
    ##     se.alpha <- rep(NA,nalpha)
    ##     se.sigma <- rep(NA,nalpha)
  }

  sol <- ou_lik_fn(
    tree=tree,
    alpha=sym_par(sqrt.alpha),
    sigma=sym_par(sigma),
    beta=beta,
    dat=dat
  )
  theta.x <- sol$coeff
  reg <- sets_of_regimes(tree,regimes)
  theta <- vector('list',nchar)
  names(theta) <- names(data)
  count <- 1
  for (n in seq_len(nchar)) {
    theta[[n]] <- theta.x[seq(from=count,length=length(reg[[n]]),by=1)]
    names(theta[[n]]) <- as.character(reg[[n]])
    count <- count+length(reg[[n]])
  }
  
  new(
    'hansentree',
    as(tree,'ouchtree'),
    call=match.call(),
    nchar=nchar,
    optim.diagn=optim.diagn,
    hessian=hs,
    data=as.list(data),
    regimes=as.list(regimes),
    beta=beta,
    theta=theta,
    sigma=sigma,
    sqrt.alpha=sqrt.alpha,
    loglik=-0.5*sol$deviance
  )
}

## note that, on input, alpha and sigma are full symmetric matrices
ou_lik_fn <- function (tree, alpha, sigma, beta, dat) {
  n <- length(dat)
  ev <- eigen(alpha,symmetric=TRUE)
  w <- .Call(ouch_weights,object=tree,lambda=ev$values,S=ev$vectors,beta=beta)
  v <- .Call(ouch_covar,object=tree,lambda=ev$values,S=ev$vectors,sigma.sq=sigma)
  gsol <- try(
    glssoln(w,dat,v),
    silent=FALSE
  )
  if (inherits(gsol,'try-error')) { # return Inf deviance (so that optimizer can keep trying)
    e <- rep(NA,n)
    theta <- rep(NA,ncol(w))
    dev <- Inf
  } else {                              # return finite deviance
    e <- gsol$residuals
    theta <- gsol$coeff
    q <- e%*%solve(v,e)
    det.v <- determinant(v,logarithm=TRUE)
    if (det.v$sign!=1)
      pStop("ou_lik_fn","non-positive determinant.")
    dev <- n*log(2*pi)+as.numeric(det.v$modulus)+q[1,1]
  }
  list(
    deviance=dev,
    coeff=theta,
    weight=w,
    vcov=v,
    resids=e
  )
}

sym_par <- function (x) {
  nchar <- floor(sqrt(2*length(x)))
  if (nchar*(nchar+1)!=2*length(x)) {
    pStop_("a symmetric matrix is parameterized by a triangular number of parameters.") #nocov
  }
  y <- matrix(0,nchar,nchar)
  y[lower.tri(y,diag=TRUE)] <- x
  y%*%t(y)
}

## sym.unpar <- function (x) {
##   y <- t(chol(x))
##   y[lower.tri(y,diag=TRUE)]
## }

sets_of_regimes <- function (object, regimes) {
  lapply(regimes,function(x)sort(unique(x)))
}

regime_spec <- function (object, regimes) {
  nterm <- object@nterm
  nchar <- length(regimes)
  reg <- sets_of_regimes(object,regimes)
  nreg <- sapply(reg,length)
  beta <- vector(mode='list',length=nterm)
  for (i in seq_len(nterm)) {
    p <- object@lineages[[object@term[i]]]
    np <- length(p)
    beta[[i]] <- vector(mode='list',length=nchar)
    for (n in seq_len(nchar)) {
      beta[[i]][[n]] <- matrix(data=NA,nrow=np,ncol=nreg[n])
      for (ell in seq_len(nreg[n])) {
        beta[[i]][[n]][,ell] <- ifelse(regimes[[n]][p]==reg[[n]][ell],1,0)
      }
    }
  }
  beta
}

## Solve the matrix equation
##   A . X + X . A = B
## for X, where we have assumed A = A'.
## 
## sym.solve <- function (a, b) {
##   n <- nrow(a)
##   d <- array(data=0,dim=c(n,n,n,n))
##   for (k in seq_len(n)) {
##     d[k,,k,] <- d[k,,k,] + a
##     d[,k,,k] <- d[,k,,k] + a
##   }
##   dim(b) <- n*n
##   dim(d) <- c(n*n,n*n)
##   x <- solve(d,b)
##   dim(x) <- c(n,n)
##   x
## }

hansen_deviate <- function (n = 1, object) {
  ev <- eigen(sym_par(object@sqrt.alpha),symmetric=TRUE)
  w <- .Call(ouch_weights,object=object,lambda=ev$values,S=ev$vectors,beta=object@beta)
  v <- .Call(ouch_covar,object=object,lambda=ev$values,S=ev$vectors,sigma.sq=sym_par(object@sigma))
  X <- array(
    data=NA,
    dim=c(object@nnodes,object@nchar,n),
    dimnames=list(
      object@nodes,
      names(object@data),
      paste('rep',seq(n),sep='.')
    )
  )

  theta <- do.call(c,object@theta)

  X[object@term,,] <- array(
    data=rmvnorm(
      n=n,
      mean=as.numeric(w%*%theta),
      var=v
    ),
    dim=c(object@nterm,object@nchar,n)
  )
  apply(X,3,as.data.frame)
}

#' @rdname coef
#' @include coef.R
#' @importFrom stats coef
#' @return \code{coef} applied to a \code{hansentree} object returns a named list containing the estimated \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2} matrices(given as the \code{alpha.matrix} and \code{sigma.sq.matrix} elements, respectively) but also the MLE returned by the optimizer
#' (as \code{sqrt.alpha} and \code{sigma}, respectively).
#' \strong{The latter elements should not be interpreted, but can be used to restart the algorithm, etc.}
#' @export
setMethod(
  'coef',
  signature=signature(object='hansentree'),
  function (object, ...) {
    list(
      sqrt.alpha=object@sqrt.alpha,
      sigma=object@sigma,
      theta=object@theta,
      alpha.matrix=sym_par(object@sqrt.alpha),
      sigma.sq.matrix=sym_par(object@sigma)
    )
  }
)

#' @rdname logLik
#' @include logLik.R
#' @importFrom stats logLik
#' @export
setMethod(
  "logLik",
  signature=signature(object='hansentree'),
  function (object) object@loglik
)

#' @rdname summary
#' @include summary.R
#' @return \code{summary} applied to a \code{hansentree} method displays the estimated \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2} matrices as well as various quantities describing the goodness of model fit.
#' @export
setMethod(
  "summary",
  signature=signature(object='hansentree'),
  function (object, ...) {
    cf <- coef(object)
    ## if (length(object@hessian)>0)
    ##   se <- sqrt(diag(solve(0.5*object@hessian)))
    dof <- length(object@sqrt.alpha)+length(object@sigma)+sum(sapply(object@theta,length))
    deviance=-2*logLik(object)
    aic <- deviance+2*dof
    aic.c <- aic+2*dof*(dof+1)/(object@nterm*object@nchar-dof-1)
    sic <- deviance+log(object@nterm*object@nchar)*dof
    list(
      call=object@call,
      conv.code=object@optim.diagn$convergence,
      optimizer.message=object@optim.diagn$message,
      alpha=cf$alpha.matrix,
      sigma.squared=cf$sigma.sq.matrix,
      optima=cf$theta,
      loglik=logLik(object),
      deviance=deviance,
      aic=aic,
      aic.c=aic.c,
      sic=sic,
      dof=dof
    )
  }
)

#' @rdname print
#' @include print.R
#' @return \code{print} displays the tree as a table, along with the coefficients of the fitted model and diagnostic information.
#' @export
setMethod(
  'print',
  signature=signature(x='hansentree'),
  function (x, ...) {
    cat("\ncall:\n")
    print(x@call)
    print(as(x,'data.frame'),...)
    if (length(x@optim.diagn)>0) {
      if (x@optim.diagn$convergence!=0)
        cat("\n",sQuote("optim")," convergence code: ",x@optim.diagn$convergence)
      if (!is.null(x@optim.diagn$message))
        cat("\n",sQuote("optim")," diagnostic message: ",x@optim.diagn$message)
    }
    sm <- summary(x)
    cat('\nalpha:\n')
    print(sm$alpha)
    cat('\nsigma squared:\n')
    print(sm$sigma.squared)
    cat('\ntheta:\n')
    print(sm$optima)
    print(unlist(sm[c("loglik","deviance","aic","aic.c","sic","dof")]))
    invisible(x)
  }
)

#' @rdname print
#' @include print.R
#' @export
setMethod(
  'show',
  signature=signature(object='hansentree'),
  function (object) {
    print(as(object,'hansentree'))
    invisible(NULL)
  }
)

#' @rdname plot
#' @include plot.R
#' @export
setMethod(
  "plot",
  signature=signature(x="hansentree"),
  function (x, ..., regimes, legend = TRUE) {
    if (missing(regimes)) regimes <- x@regimes
    f <- getMethod("plot","ouchtree")
    f(x=x,regimes=regimes,legend=legend,...)
  }
)

#' @rdname simulate
#' @include simulate.R package.R
#' @importFrom stats runif
#' @export
setMethod(
  'simulate',
  signature=signature(object='hansentree'),
  function (object, nsim = 1, seed = NULL, ...) {
    seed <- freeze(seed)
    X <- hansen_deviate(n=nsim,object)
    thaw(seed)
    X
  }
)

#' @rdname update
#' @include update.R
#' @importFrom stats update
#' @inheritParams hansen
#' @export
setMethod(
  'update',
  signature=signature(object='hansentree'),
  function (object, data, regimes, sqrt.alpha, sigma, ...) {
    if (missing(sqrt.alpha)) sqrt.alpha <- object@sqrt.alpha
    if (missing(sigma)) sigma <- object@sigma
    hansen(
      data=data,
      tree=object,
      regimes=regimes,
      sqrt.alpha=sqrt.alpha,
      sigma=sigma,
      ...
    )
  }
)

#' @rdname bootstrap
#' @include bootstrap.R
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
