\name{ape2ouch}
\alias{ape2ouch}
\title{Convert an "ape" tree to an "ouch" tree.}
\description{
  \code{ape2ouch} translates \pkg{ape}'s \code{phylo} representation of a phylogenetic tree into \pkg{ouch}'s \code{ouchtree} representation.
  The user can change the branch lengths while preserving the topology.
}
\usage{
ape2ouch(tree, scale = TRUE, branch.lengths = tree$edge.length)
}
\arguments{
  \item{tree}{
    a tree of class \code{phylo} created in package \pkg{ape}.
  }
  \item{scale}{
    if \code{scale=TRUE}, the tree's depth will be scaled to 1.
    If \code{scale} is a number, then the branch lengths will be scaled by this number.
  }
  \item{branch.lengths}{
    optional vector of branch lengths.
  }
}
\author{Aaron A. King <kingaa at umich dot edu>}
\keyword{models}
