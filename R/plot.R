##' \pkg{ouch} plotting functions
##'
##' Plot phylogenetic trees, with or without regime paintings.
##'
##' @name plot
##' @rdname plot
##' @family methods for ouch trees
##' @param x object to plot.
##' @param regimes factor or character; a vector of regime paintings.
##' @param node.names logical; should node names be displayed?
##' @param ladderize logical; should the tree be ladderized?
##' @param legend logical; display a legend?
##' @param palette function or character; specifies the colors to be used for the several regimes on the tree.
##' Specified as a function, when given an integer, `n`, the function should create a vector of `n` colors.
##' See, for example [`rainbow`][grDevices::rainbow()].
##' One can also specify the `n` colors as a vector of color codes.
##' There must be at least as many colors as levels in the `regimes`.
##' @param labels character; taxon labels.
##' @param margin numeric; width of the right margin (as a fraction of the plot width).
##' Adjust this if labels are clipped (see Examples below).
##' One can also adjust the width of the left margin (for example to aid in the formatting of the figure legend).
##' To do this, furnish `margin=c(L, R)`, where `L` and `R` are the widths of the right and left margins, respectively, as fractions of the plot width.
##' Obviously, in this case, we must have `L+R<1`.
##' @param text_opts options for the labels; passed to [`text`][graphics::text()].
##' @param legend_opts options for the the legend; passed to [`legend`][graphics::legend()].
##' @param ... additional arguments, passed to [`plot`][graphics::plot()].
##' 
##' @inheritParams graphics::plot.default
##' @importFrom graphics plot text legend par
##' @importFrom grDevices rainbow
NULL

setGeneric("plot")

tree.plot.internal <- function (
  x, ...,
  regimes,
  ladderize,
  palette,
  labels,
  legend,
  margin,
  xlab = "", ylab = "",
  xaxp = NULL,
  text_opts,
  legend_opts
) {
  ladderize <- as.logical(ladderize)
  legend <- as.logical(legend)
  rx <- range(x@times,na.rm=T)
  margin <- as.numeric(margin)
  if (length(margin)>2L||length(margin)<1L)
    pStop("plot",sQuote("margin")," should be one or two numbers.")
  if (!isTRUE(all(margin>=0 & margin<1)))
    pStop("plot",sQuote("margin")," should be between 0 and 1.")
  if (!isTRUE(sum(margin)<1))
    pStop("plot",sQuote("margin")," sum >= 1!")
  if (length(margin)>1L) {
    lmargin <- margin[1L]
    margin <- margin[2L]
  } else {
    lmargin <- 0
  }
  rxd <- margin*diff(rx)/(1-margin)
  lxd <- lmargin*diff(rx)/(1-lmargin)
  anc <- x@anc.numbers
  root <- which(is.root.node(anc))
  if (is.null(regimes)) {
    regimes <- factor(rep('unspec',length(anc)))
    names(regimes) <- x@nodes
  } else if (length(regimes)!=x@nnodes)
    pStop("plot","there must be one entry in ",sQuote("regimes")," per node of the tree.")
  if (is.null(names(regimes)) || !setequal(names(regimes),x@nodes))
    pStop("plot","regime specifications must have names corresponding to the node names.")
  regimes <- regimes[x@nodes]
  levs <- levels(as.factor(regimes))
  ## if the root is the only one with a certain regime, toss that regime out
  if (sum(regimes%in%regimes[root])==1)
    levs <- setdiff(levs,regimes[root])
  if (is.function(palette))
    palette <- palette(length(levs))
  else if (!(is.character(palette) && length(palette)>=length(levs)))
    pStop("plot",sQuote("palette")," must be either a function or a character vector of length >= ",length(levs),".")
  if (ladderize) {
    cs <- clade_size(root,anc)
  } else {
    cs <- NULL
  }
  xx <- x@times
  yy <- arrange_tree(root,anc,cs)/(length(x@term)+1)
  for (r in seq_along(levs)) {
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
    if (r > 1) par(new=TRUE)
    graphics::plot.default(
            X,Y,
            type='l',col=C,
            xaxp=if (is.null(xaxp)) c(rx,1) else xaxp,
            yaxt='n',
            xlab=xlab,ylab=ylab,
            xlim=rx+c(-lxd,rxd),ylim=c(0,1),
            ...
          )
    if (!is.null(labels)) {
      do.call(
        graphics::text,
        c(
          list(
            x=X[seq.int(1,length(X),6)],
            y=Y[seq.int(1,length(Y),6)],
            labels=labels[f],
            pos=4
            ),
          text_opts
        )
      )
    }
  }
  if (legend) {
    do.call(
      graphics::legend,
      c(
        list(x="topleft",y=levs,col=palette,lwd=1,bty="n"),
        legend_opts
      )
    )
  }
  invisible(NULL)
}

clade_size <- function (root, anc, n = integer(length(anc))) {
  children <- which(anc==root)
  if (length(children) > 0) {
    for (child in children) {
      n <- clade_size(child,anc,n)
    }
    n[root] <- sum(n[children])
  } else {
    n[root]=1
  }
  n
}

arrange_tree <- function (root, anc, cs, ypos = numeric(length(anc))) {
  children <- which(anc==root)
  if (length(children) > 0) {
    if (!is.null(cs)) {
      children <- children[order(cs[children])]
    }
    for (child in children) {
      ypos <- arrange_tree(child,anc,cs,ypos)
    }
    ypos[root] <- mean(ypos[children])
  } else {
    ypos[root] <- max(ypos)+1
  }
  ypos
}

##' @rdname plot
##' @importFrom graphics par
##' @importFrom grDevices rainbow
##' @importFrom stats setNames
##' @example examples/bimac2.R
##' @export
setMethod(
  "plot",
  signature=signature(x="ouchtree"),
  function (
    x, ..., regimes = NULL, ladderize = TRUE,
    node.names = FALSE,
    legend = !is.null(regimes), labels, frame.plot = FALSE,
    palette = rainbow,
    margin = 0.1,
    text_opts = list(),
    legend_opts = list()
  ) {
    if (missing(labels)) labels <- x@nodelabels
    if (node.names) {
      lbld <- !is.na(labels)
      labels[lbld] <- paste(x@nodes[lbld],labels[lbld])
      labels[!lbld] <- x@nodes[!lbld]
    }
    if (is.data.frame(regimes)) {
      nm <- rownames(regimes)
      regimes <- lapply(as.list(regimes),function(x){names(x)<-nm;x})
    }
    if (!(is.list(regimes) || is.null(regimes))) {
      if (length(regimes)!=x@nnodes)
        pStop("plot","there must be one entry in ",sQuote("regimes")," per node of the tree.")
      nm <- deparse(substitute(regimes))[1]
      regimes <- list(regimes)
      names(regimes) <- nm
    }
    if (is.null(regimes)) {
      tree.plot.internal(
        x,regimes=NULL,ladderize=ladderize,
        palette=palette,labels=labels,legend=legend,frame.plot=frame.plot,
        ...,margin=margin,text_opts=text_opts,legend_opts=legend_opts
      )
    } else {
      oldpar <- par(mfrow=c(1,length(regimes)))
      on.exit(par(oldpar))
      for (r in regimes) {
        tree.plot.internal(
          x,regimes=r,ladderize=ladderize,
          palette=palette,labels=labels,legend=legend,frame.plot=frame.plot,
          ...,margin=margin,text_opts=text_opts,legend_opts=legend_opts
        )
      }
    }
    invisible(NULL)
  }
)
