brown.fit <- function (data, node, ancestor, times) {
  .Deprecated('brown',package='ouch')
  pt <- parse.tree(node,ancestor,times)
  n <- pt$N
  v <- pt$branch.times
  w <- matrix(data=1,nrow=pt$N,ncol=1)
  dat <- data[pt$term]
  no.dats <- which(is.na(dat))
  if (length(no.dats) > 0)
    stop("Missing data on terminal nodes: ",node[pt$term[no.dats]])
  g <- glssoln(w,dat,v)
  theta <- g$coeff
  e <- g$residuals
  sigma <- sqrt((e %*% solve(v,e))/n)
  dim(sigma) <- 1
  u = n * (1 + log(2*pi*sigma*sigma)) + log(det(v))
  dim(u) <- 1
  df <- 2
  list(
       sigma=sigma,
       theta=theta,
       loglik=-u/2,
       deviance=u,
       aic=u+2*df,
       sic=u+log(n)*df,
       df=df
       )
}

brown.dev <- function(n = 1, node, ancestor, times, sigma, theta) {
  .Deprecated('simulate',package='ouch')
  pt <- parse.tree(node,ancestor,times)
  v <- pt$branch.times
  x <- t(rmvnorm(n, rep(theta,dim(v)[1]), as.numeric(sigma^2)*v))
  do.call(
          c,
          apply(
                x,
                1,
                function (z) {
                  y <- rep(NA,length(node))
                  y[pt$term] <- z
                  list(y)
                }
                )
          )
}

hansen.fit <- function (data, node, ancestor, times,
                        regimes = NULL,
                        interval = c(0,100),
                        tol = 1e-12) {
  .Deprecated('hansen',package='ouch')
  pt <- parse.tree(node,ancestor,times,regimes)
  n <- pt$N
  dat <- data[pt$term]
  no.dats <- which(is.na(dat))
  if (length(no.dats) > 0)
    stop("Missing data on terminal nodes: ",node[pt$term[no.dats]])
  r <- optimize(
                f=badness,
                interval=interval,
                tol=tol,
                maximum=F,
                data=dat,
                parsed.tree=pt
                )
  alpha <- r$minimum
  w <- weight.matrix.deprec(alpha,pt)
  v <- scaled.covariance.matrix(alpha,pt)
  g <- glssoln(w,dat,v)
  theta <- g$coeff
  if (pt$R>0) {
    names(theta) <- c('0',as.character(pt$regime.set))
  }
  e <- g$residuals
  sigma <- sqrt((e %*% solve(v,e))/n)
  dim(sigma) <- 1
  u = r$objective
  dim(u) <- 1
  df <- pt$R+3
  list(
       alpha=alpha,
       sigma=sigma,
       theta=theta,
       loglik=-u/2,
       deviance=u,
       aic=u+2*df,
       sic=u+log(n)*df,
       df=df
       )
}

hansen.prof <- function (alpha,
                         data, node, ancestor, times,
                         regimes = NULL) {
  .Deprecated('hansen',package='ouch')
  pt <- parse.tree(node,ancestor,times,regimes)
  n <- pt$N
  dat <- data[pt$term]
  no.dats <- which(is.na(dat))
  if (length(no.dats) > 0)
    stop("Missing data on terminal nodes: ",node[pt$term[no.dats]])
  u <- sapply(alpha,badness,data=dat,parsed.tree=pt)
  df <- pt$R+3
  list(
       alpha=alpha,
       loglik=-u/2,
       deviance=u,
       aic=u+2*df,
       sic=u+log(n)*df
       )
}

hansen.dev <- function(n = 1, node, ancestor, times, regimes = NULL, alpha, sigma, theta) {
  .Deprecated('simulate',package='ouch')
  pt <- parse.tree(node,ancestor,times,regimes)
  w <- weight.matrix.deprec(alpha,pt)
  v <- scaled.covariance.matrix(alpha,pt)
  x <- t(rmvnorm(n,as.vector(w%*%theta),as.numeric(sigma^2)*v))
  do.call(
          c,
          apply(
                x,
                1,
                function (z) {
                  y <- rep(NA,length(node))
                  y[pt$term] <- z
                  list(y)
                }
                )
          )
}

badness <- function (alpha, data, parsed.tree) {
  a <- alpha
  n <- length(data)
  w <- weight.matrix.deprec(a,parsed.tree)
  v <- scaled.covariance.matrix(a,parsed.tree)
  g <- glssoln(w,data,v)
  e <- g$residuals
  sigmasq <- (e %*% solve(v,e)) / n
  dim(sigmasq) <- 1
  u <- n * (1 + log(2*pi*sigmasq)) + log(det(v))
  dim(u) <- 1
  u
}

weight.matrix.deprec <- function (alpha, parsed.tree) {
  N <- parsed.tree$N
  R <- parsed.tree$R
  if (R > 0) {
    tree.depth <- parsed.tree$tree.depth
    ep <- parsed.tree$epochs
    beta <- parsed.tree$beta
    w <- matrix(data=0,nrow=N,ncol=R+1)
    w[,1] <- exp(-alpha*tree.depth)
    for (i in 1:N) {
      delta <- diff(exp(alpha*(ep[[i]]-tree.depth)))
      for (k in 1:R) {
        w[i,k+1] <- -sum(delta * beta[[i+N*(k-1)]])
      }
    }
  } else {
    w <- matrix(data=1,nrow=parsed.tree$N,ncol=1)
  }
  w
}

scaled.covariance.matrix <- function (alpha, parsed.tree) {
  tree.depth <- parsed.tree$tree.depth
  bt <- parsed.tree$branch.times
  if (alpha == 0) {
    v <- bt
  } else {
    a <- 2*alpha
    if (parsed.tree$R > 0) {
      v <- exp(-a*tree.depth)*expm1(a*bt)/a
    } else {
      v <- exp(-a*(tree.depth-bt))/a
    }
  }
  v
}

set.of.regimes <- function (ancestors, regime.specs) {
  unique(regime.specs[!is.root.node(ancestors)])
}

regimes <- function (ancestors, times, regime.specs, term) {
  N <- length(term)
  reg <- set.of.regimes(ancestors,regime.specs)
  R <- length(reg)
  beta <- vector(R*N, mode="list")
  for (i in 1:N) {
    for (k in 1:R) {
      p <- pedigree.deprec(ancestors, term[i])
      n <- length(p)
      beta[[i + N*(k-1)]] <- as.integer(regime.specs[p[1:(n-1)]] == reg[k])
    }
  }    
  beta
}

parse.tree <- function (nodenames, ancestors, times, regime.specs=NULL) {
  nodenames <- as.character(nodenames)
  ancestors <- as.character(ancestors)
  if (!is.valid.ouch.tree(nodenames,ancestors,times,regime.specs))
    stop('the specified tree is not in valid ouch format')
  term <- terminal.twigs(nodenames,ancestors) # get rownumbers of terminal nodes
  N <- length(term)                     # number of terminal nodes
  anc <- ancestor.numbers(nodenames,ancestors)
  bt <- branch.times(anc,times,term)   # absolute times of branch points
  e <- epochs(anc,times,term)
  if (is.null(regime.specs)) {          # useful for BM models
    pt <- list(
               N=N,
               R=0,
               tree.depth=max(times),
               term=term,
               branch.times=bt,
               epochs=e
               )
  } else {                              # useful for Hansen models
    reg <- set.of.regimes(anc,as.factor(regime.specs))
    pt <- list(
               N=N,
               R=length(reg),
               tree.depth = max(times),
               term=term,
               branch.times=bt,
               epochs=e,
               regime.set=reg,
               beta=regimes(anc,times,as.factor(regime.specs),term)
               )
  }
  return(pt)
}

tree.plot <- function (node, ancestor, times, names = NULL, regimes = NULL) {
  .Deprecated('plot',package='ouch')

  node <- as.character(node)
  ancestor <- as.character(ancestor)
  if (!is.valid.ouch.tree(node,ancestor,times,regimes))
    stop("the tree is not in valid ouch format");

  rx <- range(times,na.rm=T)
  rxd <- 0.1*diff(rx)

  anc <- ancestor.numbers(node,ancestor)

  if (is.null(regimes))
    regimes <- factor(rep(1,length(anc)))

  levs <- levels(as.factor(regimes))
  palette <- rainbow(length(levs))

  for (r in 1:length(levs)) {
    root <- which(is.root.node(anc))
    y <- arrange.tree(root,anc)
    x <- times
    f <- which(!is.root.node(anc) & regimes == levs[r])
    pp <- anc[f]
    X <- array(data=c(x[f], x[pp], rep(NA,length(f))),dim=c(length(f),3))
    Y <- array(data=c(y[f], y[pp], rep(NA,length(f))),dim=c(length(f),3))
    oz <- array(data=1,dim=c(2,1))
    X <- kronecker(t(X),oz)
    Y <- kronecker(t(Y),oz)
    X <- X[2:length(X)]
    Y <- Y[1:(length(Y)-1)]
    C <- rep(palette[r],length(X))
    if (r > 1) par(new=T)
    par(yaxt='n')
    plot(X,Y,type='l',col=C,xlab='time',ylab='',xlim = rx + c(-rxd,rxd),ylim=c(0,1))
    if (!is.null(names))
      text(X[seq(1,length(X),6)],Y[seq(1,length(Y),6)],names[f],pos=4)
  }
}

is.valid.ouch.tree <- function (node, ancestor, times, regimes = NULL) {
  .Deprecated('ouchtree',package='ouch')
  valid <- TRUE
  node <- as.character(node)
  ancestor <- as.character(ancestor)
  n <- length(node)
  if (length(unique(node)) != n) {
    warning('node names must be unique')
    valid <- FALSE
  }
  if (
      (length(ancestor) != n) ||
      (length(times) != n)
      ) {
    warning('invalid tree: all columns must be of the same length')
    valid <- FALSE
  }
  if (!is.null(regimes) && (length(regimes) != n)) {
    warning('regimes must be the same length as the other columns')
    valid <- FALSE
  }
  root <- which(is.root.node(ancestor))
  if (length(root) != 1) {
    warning('the tree must have a unique root node, designated by its having ancestor = NA')
    return(FALSE)
  }
  term <- as.list(terminal.twigs(node,ancestor))
  if (length(term) <= 0) {
    warning("there ought to be at least one terminal node, don't you think?")
    valid <- FALSE
  }
  outs <- which((!is.root.node(ancestor) & !(ancestor %in% node)))
  if (length(outs) > 0) {
    str <- sprintf("the ancestor of node %s is not in the tree\n", node[outs])
    warning(str,call.=F)
    valid <- FALSE
  }
  anc <- ancestor.numbers(node,ancestor)
  ck <- all(
            sapply(
                   1:n,
                   function(x) {
                     good <- root %in% pedigree.deprec(anc,x)
                     if (!good) {
                       str <- sprintf("node %s is disconnected", node[x])
                       warning(str,call.=F)
                     }
                     good
                   }
                   )
            )
  valid && ck
}

pedigree.deprec <- function (anc, k) {
  p <- k
  k <- anc[k]
  while (!is.root.node(k)) {
    if (k %in% p)
      stop('this is no tree: circularity detected at node ',k,call.=FALSE)
    p <- c(p,k)
    k <- anc[k]
  }
  p
}

