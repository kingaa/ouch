#' Phylogenetic tree object in 'ouch' format.
#' 
#' An object containing a phylogenetic tree in a form suitable for using \pkg{ouch} methods.
#' 
#' \code{ouchtree} creates an \code{ouchtree} object.
#' This contains the topology, branch times, and epochs.
#' It also (optionally) holds names of taxa for display purposes.
#' 
#' @name ouchtree
#' @rdname ouchtree
#' @aliases ouchtree ouchtree-class
#' 
#' @author Aaron A. King
#' @seealso \code{ouchtree}, \code{ape2ouch}, \code{brown}, \code{hansen}
#' @keywords models
#' @example examples/bimac1.R
#' 
NULL

setClass(
  'ouchtree',
  representation=representation(
    nnodes = 'integer',
    nodes = 'character',
    ancestors = 'character',
    nodelabels = 'character',
    times = 'numeric',
    root = 'integer',
    nterm = 'integer',
    term = 'integer',
    anc.numbers = 'integer',
    lineages = 'list',
    epochs = 'list',
    branch.times = 'matrix',
    depth = 'numeric'
  )
)

#' @rdname ouchtree
#' @param nodes A character vector giving the name of each node.
#' These are used internally and must be unique.
#' @param ancestors Specification of the topology of the phylogenetic tree.
#' This is in the form of a character vector specifying the name
#' (as given in the \code{nodes} argument)
#' of the immediate ancestor of each node.
#' In particular, the i-th name is that of the ancestor of the i-th node.
#' The root node is distinguished by having no ancestor (i.e., \code{NA}).
#' @param times A vector of nonnegative numbers, one per node in the tree,
#' specifying the time at which each node is located.
#' Time should be increasing from the root node to the terminal twigs.
#' @param labels Optional vector of node labels.
#' These will be used in plots to label nodes.
#' It is not necessary that these be unique.
#'
#' @include package.R
#' @export ouchtree
ouchtree <- function (nodes, ancestors, times, labels = as.character(nodes)) {

  nodes <- as.character(nodes)
  ancestors <- as.character(ancestors)

  n <- length(nodes)
  if (length(unique(nodes)) != n) stop("node names must be unique")
  if (length(ancestors) != n)
    stop("invalid tree: ",sQuote("ancestors")," must have the same length as ",sQuote("nodes"))
  if (length(times) != n) 
    stop("invalid tree: ",sQuote("times")," must have the same length as ",sQuote("nodes"))
  if (length(labels) != n)
    stop("invalid tree: ",sQuote("labels")," must be the same length as ",sQuote("nodes"))

  root <- which(is.root.node(ancestors))
  if (length(root) != 1)
    stop("invalid tree: there must be a unique root node, designated by its having ancestor = NA")
  if (times[root] != 0)
    stop("the algorithms assume that the root node is at time=0")
  
  term <- terminal.twigs(nodes,ancestors)
  if (length(term) <= 0) 
    stop("invalid tree: there ought to be at least one terminal node, don't you think?")

  outs <- which((!is.root.node(ancestors) & !(ancestors %in% nodes)))
  if (length(outs) > 0) {
    for (out in outs) {
      warning(
        sprintf(
          "the ancestor of node %s is not in the tree",
          nodes[out]
        ),
        call.=FALSE)
    }
    stop("invalid tree")
  }
  
  anc <- ancestor.numbers(nodes,ancestors)

  if (any(anc==seq(along=anc),na.rm=TRUE)) {
    w <- which(anc==seq(along=anc))
    stop("this is no tree: node ",nodes[w[1]]," is its own ancestor",call.=FALSE)
  }

###  lineages <- build.lineages(anc)
  lineages <- vector(mode='list',length=n)
  todo <- root
  k <- 1
  while (k <= length(todo)) {
    todo <- c(todo,which(anc==todo[k]))
    a <- anc[todo[k]]
    lineages[[todo[k]]] <- c(todo[k],lineages[[a]])
    if (todo[k] %in% lineages[[a]]) 
      stop("this is no tree: circularity detected at node ",nodes[todo[k]]," in ",sQuote("ouchtree"),call.=FALSE)
    k <- k+1
  }

  for (k in 1:n) {
    if (!(root %in% lineages[[k]]))
      stop("node ",nodes[k]," is disconnected",call.=FALSE)
  }

  new(
    'ouchtree',
    nnodes = length(nodes),
    nodes=nodes,
    ancestors=ancestors,
    nodelabels=as.character(labels),
    times=times,
    root=root,
    nterm=length(term),
    term=term,
    anc.numbers=anc,
    lineages=lineages,
    epochs=epochs(lineages,times,term), # absolute times of epochs
    branch.times=branch.times(lineages,times,term), # absolute times of branch points
    depth=max(times)
  )
}

## map ancestor names to row numbers
ancestor.numbers <- function (nodenames, ancestors) { 
  sapply(ancestors,function(x)charmatch(x,nodenames),USE.NAMES=FALSE)
}

## nodenames of terminal twigs (terminal nodes are not ancestors)
terminal.twigs <- function (nodenames, ancestors) {
  which(nodenames %in% setdiff(nodenames,unique(ancestors)))
}

build.lineages <- function (ancestors) {
  n <- length(ancestors)
  lineages <- vector(mode='list',length=n)
  pedigree <- function (k) {
    if (is.null(lineages[[k]])) {
      a <- ancestors[k]
      if (is.root.node(a)) {
        lineages[[k]] <<- k
      } else {
        if (is.null(lineages[[a]])) Recall(a)
        if (k %in% lineages[[a]]) 
          stop('this is no tree: circularity detected at node ',k,call.=FALSE)
        lineages[[k]] <<- c(k,lineages[[a]])
      }
    }
    NULL
  }
  for (k in 1:n)
    pedigree(k)
  lineages
}

branch.times <- function (lineages, times, term) {
  N <- length(term)
  bt <- matrix(data=0,nrow=N,ncol=N)
  for (i in 1:N) {
    pedi <- lineages[[term[i]]]
    for (j in seq_len(i-1)) {
      pedj <- lineages[[term[j]]]
      for (k in 1:length(pedi)) {
        if (any(pedj == pedi[k])) break
      }
      bt[j,i] <- bt[i,j] <- times[pedi[k]]
    }
    bt[i,i] <- times[term[i]]
  }
  bt
}

epochs <- function (lineages, times, term) {
  N <- length(term)
  e <- vector(mode='list',length=N)
  for (i in 1:N) {
    pedi <- lineages[[term[i]]]
    e[[i]] <- times[pedi]
  }
  e
}

is.root.node <- function (anc) {
  is.na(anc)
}

#' @rdname print
#' @include print.R
#' @export
setMethod(
  'print',
  signature=signature(x='ouchtree'),
  function (x, ...) {
    print(as(x,'data.frame'),...)
    invisible(x)
  }
)

#' @rdname print
#' @include print.R
#' @export
setMethod(
  'show',
  signature=signature(object='ouchtree'),
  function (object) {
    print(as(object,'ouchtree'))
    invisible(NULL)
  }
)

setAs(
  from='ouchtree',
  to='data.frame',
  def = function (from) {
    df <- data.frame(
      nodes=from@nodes,
      ancestors=from@ancestors,
      times=from@times,
      labels=from@nodelabels,
      row.names=from@nodes
    )
    rownames(df) <- from@nodes
    df
  }
)

#' @rdname plot
#' @include plot.R
#' @importFrom grDevices rainbow
#' @importFrom graphics par text
#' @importFrom stats setNames
#' @export
setMethod(
  "plot",
  signature=signature(x="ouchtree"),
  function (x, y, regimes = NULL, node.names = FALSE, legend = TRUE, ..., labels) {
    if (!missing(y)) warning(sQuote("y")," is ignored.")
    labels <- x@nodelabels
    if (node.names) {
      lbld <- !is.na(labels)
      labels[lbld] <- paste(x@nodes[lbld],labels[lbld])
      labels[!lbld] <- x@nodes[!lbld]
    }
    if (is.data.frame(regimes)) {
      nm <- rownames(regimes)
      regimes <- lapply(as.list(regimes),function(x){names(x)<-nm;x})
    }
    if (is.list(regimes)) {
      if (any(sapply(regimes,length)!=x@nnodes))
        stop("each element in ",sQuote("regimes")," must be a vector with one entry per node of the tree")
    } else if (!is.null(regimes)) {
      if (length(regimes)!=x@nnodes)
        stop("there must be one entry in ",sQuote("regimes")," per node of the tree")
      nm <- deparse(substitute(regimes))[1]
      regimes <- list(regimes)
      names(regimes) <- nm
    }
    if (is.null(regimes)) {
      invisible(tree.plot.internal(x,regimes=NULL,labels=labels,legend=legend,...))
    } else {
      oldpar <- par(mfrow=c(1,length(regimes)))
      on.exit(par(oldpar))
      retval <- lapply(
        regimes,
        function (r) tree.plot.internal(x,regimes=r,labels=labels,...)
      )
      invisible(retval)
    }
  }
)
