\name{brown}
\alias{brown}
\title{Fit a Brownian-motion model of evolution along a phylogenetic tree}
\description{
  \code{brown} fits a Brownian motion model to the given data.
}
\usage{
brown(data, tree)
}
\arguments{
  \item{data}{Phenotypic data for extant species, i.e., at the terminal
    ends of the phylogenetic tree.  This can either be a numeric vector
    or a list.  If it is a numeric vector, there must be one entry for
    every node.  If it is a list, it must consist entirely of numeric
    vectors, each of which has one entry per node.  A data-frame is
    coerced to a list.}
  \item{tree}{A phylogenetic tree, specified as an \code{ouchtree}
    object.}
}
\value{
  \code{brown} returns an object of class \code{browntree}.  See
  \code{\link{browntree-class}} for information on the methods of this
  class.
}
\seealso{\code{\link{ouchtree}}, \code{\link{browntree}}}
\references{
  Butler, M.A. and A.A. King (2004) Phylogenetic comparative analysis: a
  modeling approach for adaptive evolution. American Naturalist
  164:683-695.
}
\author{Aaron A. King <kingaa at umich dot edu>}
\keyword{models}
