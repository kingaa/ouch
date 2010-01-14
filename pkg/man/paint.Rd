\name{paint}
\alias{paint}
\title{Painting selective regimes on a phylogenetic tree.}
\description{
  Function to paint selective regimes on a phylogenetic tree.
}
\usage{
paint(tree, subtree, branch, which = 1)
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
  \item{which}{
    integer;
    if \code{tree} is a \code{hansentree}, start not with a blank canvas but with the regime specifications \code{tree} contains for the character indicated by \code{which}.
  }
}
\details{
  The names of \code{subtree} and \code{branch} must be the names of nodes of \code{tree}.
  The painting proceeds in a particular order: one can overpaint a branch.
  The subtrees indicated by the elements of \code{subtree} are painted first, in order.
  Then the branches indicated by \code{branch} are painted.
  If \code{tree} is a simple \code{ouchtree} object, then \code{paint} begins with a blank canvas, i.e., a tree painted with the single regime "nonspec".
  If \code{tree} inherits class \code{hansentree}, then \code{paint} begins with the regimes specified in the \code{regimes} slot of \code{tree}.
  Note that, if \code{tree} is a multivariate \code{hansentree}, then there are multiple regime specifications contained in \code{tree}.
  In this case, the argument \code{which} lets you pick which one you wish to begin with;
  by default, the first is used.
}
\value{
  A vector of class 'factor' with names corresponding to the nodes in \code{tree}, specifying selective regimes.
}
\author{Aaron A. King kingaa at umich dot edu}
\seealso{\code{ouchtree}, \code{hansen}}
\examples{
data(bimac)
x <- with(bimac,ouchtree(nodes=node,times=time/max(time),ancestors=ancestor,labels=species))
r <- paint(x,subtree=c("1"="medium","9"="large","2"="small"),branch=c("38"="large","2"="medium"))
plot(x,regimes=r,node.names=TRUE)
# compare to bimac['OU.LP']
h5 <- hansen(data=log(bimac['size']),tree=x,regimes=bimac['OU.LP'],sqrt.alpha=1,sigma=1,reltol=1e-5)
r <- paint(h5,branch=c("18"="large"),subtree=c("9"="small"))
plot(x,regimes=r,node.names=TRUE)
}
\keyword{models}
