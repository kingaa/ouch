\name{paint}
\alias{paint}
\title{Painting selective regimes on a phylogenetic tree.}
\description{
  Function to paint selective regimes on a phylogenetic tree.
}
\usage{
paint(tree, subtree, branch)
}
\arguments{
  \item{tree}{
    An object of class \code{ouchtree}.
  }
  \item{subtree}{
    An optional named vector specifying the root nodes of subtrees.
    Each branch that descends from this node will be painted with the specified regime.
  }
  \item{branch}{
    An optional named vector specifying the end nodes of branches.
    The unique branch that terminates at the named node will be painted with the specified regime.
  }
}
\details{
  The names of \code{subtree} and \code{branch} must be the names of nodes of \code{tree}.
  The painting proceeds in a particular order: one can overpaint a branch.
  The subtrees indicated by the elements of \code{subtree} are painted first, in order.
  Then the branches indicated by \code{branch} are painted.
  If \code{tree} is a simple \code{ouchtree} object, then \code{paint} begins with a blank canvas, i.e., a tree painted with the single regime "nonspec".
  If \code{tree} inherits class \code{hansentree}, then \code{paint} begins with the regimes specified in the \code{regimes} slot of \code{tree}.
}
\value{
  A vector of class 'factor' with names corresponding to the nodes in \code{tree}, specifying selective regimes.
}
\author{Aaron A. King kingaa at umich dot edu}
\seealso{\code{ouchtree}, \code{hansen}}
\examples{
data(bimac)
attach(bimac)
x <- ouchtree(nodes=node,times=time/max(time),ancestors=ancestor,labels=species)
r <- paint(x,subtree=c("1"="medium","9"="large","2"="small"),branch=c("38"="large","2"="medium"))
# compare to bimac['OU.LP']
plot(x,regimes=r,node.names=TRUE)
detach(bimac)
}
\keyword{models}
