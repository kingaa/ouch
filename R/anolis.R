#' Greater Antillean anolis lizard sexual size dimorphism data
#' 
#' The dataset consists of sexual size-dimorphism data for 38 species of anoles from Cuba, Hispaniola, Jamaica, and Puerto Rico (Butler, Schoener, and Losos 2000).
#' Each of these species belongs to one of six microhabitat types, or \dQuote{ecomorphs} (sensu Williams, 1972):
#' trunk-ground, grass-bush, trunk, trunk-crown, twig, and crown-giant.
#' The data were used to demonstrate an evolutionary association between habitat type and degree of sexual size dimorphism.
#' 
#' Size dimorphism was calcuated as the log-ratio of male snout-to-vent length to female snout-to-vent length (males are larger).
#' 
#' In this example, we tested three models of evolution:
#' Brownian motion, Ornstein-Uhlenbeck with one global optimum, and Ornstein-Uhlenbeck with 7 optima (one for each ecomorph type plus an additional one for an \dQuote{unknown} type).
#' 
#' For the 7-optima model, we assigned each terminal branch to an optimum according to the ecomorph type of the extant species.
#' Because we had no information to help guide hypotheses about internal branches, we assigned
#' internal branches to the \dQuote{unknown} selective regime.
#' The phylogeny of these species is consistent with and adaptive radiation, with a burst of speciation events early in the evolutionary history of this clade (see phylogeny in Butler & King (2004) or example below).
#' 
#' @name anolis.ssd
#' @rdname anolis_ssd
#' @docType data
#' @format
#' A data frame with 38 observations on the following 6 variables.
#' \describe{
#'   \item{node}{Labels for the nodes.}
#'   \item{species}{Names of extant species.}
#'   \item{log.SSD}{Log sexual size dimorphism of extant species.}
#'   \item{ancestor}{Ancestor node.}
#'   \item{time}{Time of node.}
#'   \item{OU.1}{a factor with levels \code{ns}}
#'   \item{OU.7}{a factor with levels corresponding to ecomorph
#'     (\code{tg} \code{tc} \code{gb} \code{cg} \code{tw} \code{tr} \code{anc})}
#' }
#' @author Marguerite A. Butler, Aaron A. King
#' @references
#' \Butler2000
#'
#' \Williams1972
#' 
#' @source \Butler2004
#' @keywords models
#' @example examples/anolis.R
#' 
NULL
