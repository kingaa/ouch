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
