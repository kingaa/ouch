\name{ouchtree}
\docType{class}
\alias{ouchtree}
\alias{ouchtree-class}
\alias{show,ouchtree-method}
\alias{show-ouchtree}
\alias{plot,ouchtree-method}
\alias{plot-ouchtree}
\alias{print,ouchtree-method}
\alias{print-ouchtree}
\alias{coerce,ouchtree,data.frame-method}
\title{Phylogenetic tree object in 'ouch' format.}
\description{
  An object containing a phylogenetic tree in a form suitable for using \pkg{ouch} methods.
}
\usage{
ouchtree(nodes, ancestors, times, labels = as.character(nodes))
}
\arguments{
  \item{nodes}{
    A character vector giving the name of each node.  These are used internally and must be unique.
  }
  \item{ancestors}{
    Specification of the topology of the phylogenetic tree.
    This is in the form of a character vector specifying the name (as given in the \code{nodes} argument) of the immediate ancestor of each node.
    In particular, the i-th name is that of the ancestor of the i-th node.
    The root node is distinguished by having no ancestor (i.e., NA).
  }
  \item{times}{
    A vector of nonnegative numbers, one per node in the tree, specifying the time at which each node is located.
    Time should be increasing from the root node to the terminal twigs.
  }
  \item{labels}{
    Optional vector of node labels.
    These will be used in plots to label nodes.
    It is not necessary that these be unique.
  }
}
\details{
  \code{ouchtree} creates an \code{ouchtree} object.
  This contains the topology, branch times, and epochs.
  It also holds (optionally) names of taxa for display purposes.
}
\section{Methods}{
  \describe{
    \item{\code{plot(tree,regimes=NULL,node.names=FALSE,legend=TRUE,...)}}{
      displays the phylogenetic tree graphically.
    }
    \item{\code{print()}}{
      displays the tree in table form.
    }
    \item{\code{show()}}{
      displays the tree.
    }
    \item{coerce}{
      An \code{ouchtree} object can be coerced to a data-frame via \code{as(object,"data.frame")}.
    }
  }
}
\author{Aaron A. King kingaa at umich dot edu}
\seealso{\code{ouchtree}, \code{ape2ouch}, \code{brown}, \code{hansen}}
\examples{
data(bimac)
tree <- with(bimac,ouchtree(nodes=node,ancestors=ancestor,times=time,labels=species))
plot(tree)
plot(tree,node.names=TRUE)
print(tree)
}
\keyword{models}
