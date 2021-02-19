#' Plotting functions.
#'
#' Plot phylogenies.
#'
#' @name plot
#' @rdname plot
#' @family methods
#' @param x object to plot.
#' @param y ignored.
#' @param regimes factor or character; a vector of regime paintings.
#' @param node.names logical; should node names be displayed?
#' @param legend logical; display a legend?
#' @param labels character; taxon labels.
#' @param ... additional arguments, passed to \code{text}.
NULL

tree.plot.internal <- function (x, regimes = NULL, labels = x@nodelabels, legend = TRUE, ...) {
  rx <- range(x@times,na.rm=T)
  rxd <- 0.1*diff(rx)
  anc <- x@anc.numbers
  term <- x@term
  root <- which(is.root.node(anc))
  if (is.null(regimes)) {
    regimes <- factor(rep('unspec',length(anc)))
    names(regimes) <- x@nodes
  } else if (length(regimes)!=x@nnodes)
    stop(sQuote("regimes")," must be a vector with one entry for each node")
  if ((is.null(names(regimes)))||(!setequal(names(regimes),x@nodes)))
    stop("regime specifications must have names corresponding to the node names")
  regimes <- regimes[x@nodes]
  levs <- levels(as.factor(regimes))
  ## if the root is the only one with a certain regime, toss that regime out
  if (sum(regimes%in%regimes[root])==1)
    levs <- setdiff(levs,regimes[root])
  palette <- rainbow(length(levs))
  for (r in seq_along(levs)) {
    yy <- arrange.tree(root,anc,term)
    xx <- x@times
    f <- which(!is.root.node(anc) & regimes==levs[r])
    pp <- anc[f]
    X <- array(data=c(xx[f],xx[pp],rep(NA,length(f))),dim=c(length(f),3))
    Y <- array(data=c(yy[f],yy[pp],rep(NA,length(f))),dim=c(length(f),3))
    oz <- array(data=1,dim=c(2,1))
    X <- kronecker(t(X),oz)
    Y <- kronecker(t(Y),oz)
    X <- X[-1]
    Y <- Y[-length(Y)]
    C <- rep(palette[r],length(X))
    if (r > 1) par(new=T)
    par(yaxt='n')
    plot(X,Y,type='l',col=C,xlab='time',ylab='',xlim=rx+c(-rxd,rxd),ylim=c(0,1),...)
    if (!is.null(labels))
      text(X[seq(1,length(X),6)],Y[seq(1,length(Y),6)],labels[f],pos=4,...)
  }
  if (legend)
    legend('topleft',levs,lwd=1,col=palette,bty='n')
  invisible(NULL)
}

## The following function lays out the tree.
## Contribution by M. Butler: vertical spacing of terminal leaves is even
arrange.tree <- function (root, anc, term) {
  k <- which(anc==root) # nodes immediately descended from 'root'
  reltree <- numeric(length(anc))
  reltree[term] <- seq(from=0,to=1,length.out=length(term)) ## space terminal taxa
  if (root %in% term) { ## terminate recursion
    setNames(reltree[root],root)
  } else { ## recurse
    p <- list()
    for (j in seq_along(k)) {
      p[[j]] <- arrange.tree(k[j],anc,term)
    }
    y <- unlist(p)
    yy <- setNames(mean(y[as.character(k)]),root) # x is mean y-coordinate of descendants 
    y <- c(y,yy)
    oo <- as.character(sort(as.numeric(names(y))))  # sorted in node order
    y[oo]
  }
}
