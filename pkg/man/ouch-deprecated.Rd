\name{ouch-deprecated}
\alias{ouch-deprecated}
\alias{is.valid.ouch.tree}
\alias{tree.plot}
\alias{brown.fit}
\alias{brown.dev}
\alias{hansen.fit}
\alias{hansen.dev}
\alias{hansen.prof}
\title{Deprecated functions in the "ouch" package}
\description{
  These functions are deprecated.
}
\usage{
is.valid.ouch.tree(node, ancestor, times, regimes = NULL)
tree.plot(node, ancestor, times, names = NULL, regimes = NULL)
brown.fit(data, node, ancestor, times)
brown.dev(n = 1, node, ancestor, times, sigma, theta)
hansen.fit(data, node, ancestor, times, regimes = NULL,
           interval = c(0, 100), tol = 1e-12)
hansen.prof(alpha, data, node, ancestor, times, regimes = NULL)
hansen.dev(n = 1, node, ancestor, times, regimes = NULL,
           alpha, sigma, theta)
}
\arguments{
  \item{data}{Phenotypic data for extant species, i.e., at the
    terminal ends of the phylogenetic tree.}
  \item{node}{Specification of the names of the nodes.}
  \item{ancestor}{Specification of the topology of the phylogenetic
    tree.  This is in the form of a character vector of node names, one
    for each node in the tree.  The i-th name is that of the ancestor
    of the i-th node.  The root node is distinguished by having no
    ancestor (i.e., NA).}
  \item{times}{A vector of nonnegative numbers, one per node in the
    tree, specifying the time at which each node is located.  The root node
    should be assigned time 0.}
  \item{regimes}{A vector of codes, one for each node in the
    tree, specifying the selective regimes hypothesized to have been
    operative.  Corresponding to each node, enter the code of the regime
    hypothesized for the branch segment terminating in that node.
    For the root node, because it has no branch segment terminating on
    it, the regime specification is irrelevant.}
  \item{names}{Optional vector of species names.}
  \item{n}{Number of pseudorandom data sets to generate.}
  \item{sigma}{the value of \eqn{\sigma} to be used in the simulations.}
  \item{alpha}{Value of selection-strength parameter \eqn{\alpha}.}
  \item{theta}{the value of \eqn{\theta} to be used in the simulations.}
  \item{interval}{The interval which will be searched for the optimal
    \eqn{\alpha}.  By default, \eqn{0 < \alpha < 100}.}
  \item{tol}{Convergence tolerance.}
}
\details{
  Phylogenetic trees in \pkg{ouch} are now implemented using the S4
  class \code{ouchtree}.  The validity checks formerly implemented by
  \code{is.valid.ouch.tree} are now performed when an \code{ouchtree} is
  contructed; see \code{\link{ouchtree}}, \code{\link{ouchtree-class}}.
  The \code{ouchtree} class has a \code{plot} method, obviating the need
  for \code{tree.plot}; see \code{\link{ouchtree-class}}.

  Brownian motion (BM) models are now built around the S4 class
  \code{browntree}.  Fitting BM models is now accomplished by means of
  \code{\link{brown}}; generating random deviates from a fitted model,
  by the \code{simulate} method for \code{browntree}-class objects.  See
  \code{\link{browntree-class}}.

  The Ornstein-Uhlenbeck-based Hansen models are now implemented using
  the S4 class \code{hansentree}.  They are fit using
  \code{\link{hansen}}; see \code{\link{hansentree-class}} for
  documentation of the \code{simulate} method for this class.
}
\seealso{\code{\link{ouch-package}}}
\author{Aaron A. King <kingaa at umich dot edu>}
\keyword{models}
