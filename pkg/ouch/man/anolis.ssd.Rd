\name{anolis.ssd}
\alias{anolis.ssd}
\docType{data}
\title{Greater Antillean anolis lizard sexual size dimorphism data.}
\description{
  The dataset consists of sexual size-dimorphism data for 38 species of anoles from Cuba, Hispaniola, Jamaica, and Puerto Rico (Butler, Schoener, and Losos 2000).
  Each of these species belongs to one of six microhabitat types, or ``ecomorphs'' (sensu Williams, 1972):
  trunk-ground, grass-bush, trunk, trunk-crown, twig, and crown-giant.
  The data were used to demonstrate an evolutionary association between habitat type and degree of sexual size dimorphism.
}
\usage{data(anolis.ssd)}
\format{
  A data frame with 38 observations on the following 6 variables.
  \describe{
    \item{node}{Labels for the nodes.}
    \item{species}{Names of extant species.}
    \item{log.SSD}{Log sexual size dimorphism of extant species.}
    \item{ancestor}{Ancestor node.}
    \item{time}{Time of node.}
    \item{OU.1}{a factor with levels \code{ns}}
    \item{OU.7}{a factor with levels corresponding to ecomorph
      (\code{tg} \code{tc} \code{gb} \code{cg} \code{tw} \code{tr}
      \code{anc})}
  }
}
\details{
  Size dimorphism was calcuated as the log-ratio of male snout-to-vent length to female snout-to-vent length (males are larger).

  In this example, we tested three models of evolution:
  Brownian motion, Ornstein-Uhlenbeck with one global optimum, and Ornstein-Uhlenbeck with 7 optima (one for each ecomorph type plus an additional one for an ``unknown'' type).

  For the 7-optima model, we assigned each terminal branch to an optimum according to the ecomorph type of the extant species.
  Because we had no information to help guide hypotheses about internal branches, we assigned internal branches to the ``unknown'' selective regime.
  The phylogeny of these species is consistent with and adaptive radiation, with a burst of speciation events early in the evolutionary history of this clade (see phylogeny in Butler & King (2004) or example below.
}
\source{
  Butler, M.A. and A.A. King. 2004.
  Phylogenetic comparative analysis: a modeling approach for adaptive evolution.
  American Naturalist 164:683-695.
}
\references{
  Butler, M. A., T. W. Schoener, and J. B. Losos. 2000.
  The relationship between sexual size dimorphism and habitat use in Greater Antillean Anolis lizards.
  Evolution, 54:259-272.

  Williams, E. E. 1972.
  The origin of faunas. Evolution of lizard congeners in a complex island fauna: a trial analysis.
  Evol. Biol., 6:47-89.
}
\examples{
data(anolis.ssd)
tree <- with(anolis.ssd,ouchtree(node,ancestor,time/max(time),species))
plot(tree,node.names=TRUE)
print(h1 <- brown(anolis.ssd['log.SSD'],tree))
plot(h1)
print(h2 <- hansen(anolis.ssd['log.SSD'],tree,anolis.ssd['OU.1'],sqrt.alpha=1,sigma=1))
plot(h2)
print(h3 <- hansen(anolis.ssd['log.SSD'],tree,anolis.ssd['OU.7'],sqrt.alpha=1,sigma=1))
plot(h3)
}
\author{Marguerite A. Butler <mbutler at hawaii dot edu> and Aaron A. King <kingaa at umich dot edu>}
\keyword{models}

