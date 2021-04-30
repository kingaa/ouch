#' Painting regimes on a phylogenetic tree
#' 
#' Function to paint selective regimes on a phylogenetic tree.
#' 
#' The names of \code{subtree} and \code{branch} must be the names of nodes of \code{tree}.
#' The painting proceeds in a particular order:
#' one can overpaint a branch.
#' The subtrees indicated by the elements of \code{subtree} are painted first, in order.
#' Then the branches indicated by \code{branch} are painted.
#' If \code{tree} is a simple \code{ouchtree} object, then \code{paint} begins with a blank canvas,
#' i.e., a tree painted with the single regime \dQuote{nonspec}.
#' If \code{tree} is of class \code{hansentree}, then \code{paint} begins with the regimes specified in the \code{regimes} slot of \code{tree}.
#' Note that, if \code{tree} is a multivariate \code{hansentree}, then there are multiple regime specifications contained in \code{tree}.
#' In this case, the argument \code{which} lets you pick which one you wish to begin with;
#' by default, the first is used.
#' 
#' @param tree An object of class \code{ouchtree}.
#' @param subtree An optional named vector specifying the root nodes of subtrees.
#' Each branch that descends from this node will be painted with the
#' specified regime.
#' @param branch An optional named vector specifying the end nodes of branches.
#' The unique branch that terminates at the named node will be painted with the specified regime.
#' @param which integer;
#' if \code{tree} is a \code{hansentree}, start not with a blank canvas but with the regime specifications \code{tree} contains for the character indicated by \code{which}.
#' @return A vector of class \sQuote{factor} with names corresponding to the nodes in \code{tree}, specifying selective regimes.
#' @author Aaron A. King
#' @seealso \code{ouchtree}, \code{hansen}
#' @keywords models
#' @example examples/paint.R
#' @include ouchtree.R
#' @importFrom utils head tail
#' @export paint
paint <- function (tree, subtree, branch, which = 1) {
  if (!is(tree,'ouchtree'))
    stop(sQuote("tree")," must be of class ",sQuote("ouchtree"))
  if (is(tree,'hansentree')) {
    regimes <- try(tree@regimes[[which]],silent=FALSE)
    if (inherits(regimes,'try-error'))
      stop(sQuote("paint")," error: invalid ",sQuote("which"))
  } else {
    regimes <- rep('unspec',length(tree@nodes))
    names(regimes) <- tree@nodes
  }
  if (!missing(subtree)) {
    st.nm <- names(subtree)
    if (is.null(st.nm))
      stop(sQuote("subtree")," must be a named vector")
    if (!all(st.nm%in%tree@nodes))
      stop("all names of ",sQuote("subtree")," must be names of nodes of ",sQuote("tree"))
    subtree <- as.character(subtree)
  } else {
    subtree <- character(0)
    st.nm <- character(0)
  }
  if (!missing(branch)) {
    br.nm <- names(branch)
    if(is.null(br.nm))
      stop(sQuote("branch")," must be a named vector")
    if (!all(br.nm%in%tree@nodes))
      stop("all names of ",sQuote("branch")," must be names of nodes of ",sQuote("tree"))
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
