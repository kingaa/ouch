##' Painting regimes on a phylogenetic tree
##'
##' Function to paint selective regimes on a phylogenetic tree.
##'
##' The names of `subtree` and `branch` must be the names of nodes of `tree`.
##' The painting proceeds in a particular order:
##' one can overpaint a branch.
##' The subtrees indicated by the elements of `subtree` are painted first, in order.
##' Then the branches indicated by `branch` are painted.
##' If `tree` is of class [`hansentree`][hansen], then `paint` begins with the regimes specified in the `regimes` slot of `tree`.
##' Otherwise, `paint` begins with a blank canvas,
##' i.e., a tree painted with the single regime "nonspec".
##' Note that, if `tree` is a multivariate `hansentree`, then there are multiple regime specifications contained in `tree`.
##' In this case, the argument `which` lets you pick which one you wish to begin with;
##' by default, the first is used.
##'
##' @keywords models
##' @family methods for ouch trees
##' @family phylogenetic comparative models
##' @include ouchtree.R
##'
##' @param tree An object of class [`ouchtree`].
##' @param subtree An optional named vector specifying the root nodes of subtrees.
##' Each branch that descends from this node will be painted with the specified regime.
##' @param branch An optional named vector specifying the end nodes of branches.
##' The unique branch that terminates at the named node will be painted with the specified regime.
##' @param which integer;
##' if `tree` is a [`hansentree`][hansen], start not with a blank canvas but with the regime specifications `tree` contains for the character indicated by `which`.
##' @return A vector of class 'factor' with names corresponding to the nodes in `tree`, specifying selective regimes.
##' @author Aaron A. King
##' @example examples/paint.R
##' @importFrom utils head tail
##' @export
paint <- function (tree, subtree, branch, which = 1) {
  if (!is(tree,'ouchtree'))
    pStop("paint",sQuote("tree")," must be of class ",sQuote("ouchtree"),".")
  if (is(tree,'hansentree')) {
    regimes <- try(tree@regimes[[which]],silent=FALSE)
    if (inherits(regimes,'try-error'))
      pStop("paint",sQuote("paint")," error: invalid ",sQuote("which"),".")
  } else {
    regimes <- rep('unspec',length(tree@nodes))
    names(regimes) <- tree@nodes
  }
  if (!missing(subtree)) {
    st.nm <- names(subtree)
    if (is.null(st.nm))
      pStop("paint",sQuote("subtree")," must be a named vector.")
    if (!all(st.nm%in%tree@nodes))
      pStop("paint","all names of ",sQuote("subtree")," must be names of nodes of ",sQuote("tree"),".")
    subtree <- as.character(subtree)
  } else {
    subtree <- character(0)
    st.nm <- character(0)
  }
  if (!missing(branch)) {
    br.nm <- names(branch)
    if(is.null(br.nm))
      pStop("paint",sQuote("branch")," must be a named vector.")
    if (!all(br.nm%in%tree@nodes))
      pStop("paint","all names of ",sQuote("branch")," must be names of nodes of ",sQuote("tree"),".")
    branch <- as.character(branch)
  } else {
    branch <- character(0)
    br.nm <- character(0)
  }
  tog <- as.factor(c(as.character(subtree),as.character(branch)))
  subtree <- head(tog,length(subtree))
  branch <- tail(tog,length(branch))
  names(subtree) <- st.nm
  names(branch) <- br.nm
  for (k in seq(along=subtree)) {
    st <- sapply(tree@lineages,function(x,y)y%in%x[-1],y=st.nm[k])
    regimes[st] <- as.character(subtree[k])
  }
  for (k in seq(along=branch)) {
    br <- tree@nodes%in%br.nm[k]
    regimes[br] <- as.character(branch[k])
  }
  as.factor(regimes)
}
