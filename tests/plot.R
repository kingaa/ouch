library(ouch)
set.seed(1106216184L)

png(filename="plot%02d.png")

tree <- with(bimac,ouchtree(node,ancestor,time/max(time),spcode))
plot(tree,node.names=TRUE)

plot(tree,regimes=bimac["OU.3"],node.names=TRUE,ladderize=TRUE)
plot(tree,regimes=bimac["OU.3"],node.names=FALSE,ladderize=FALSE)

bimac <- bimac[order(runif(n=nrow(bimac))),]
tree <- with(bimac,ouchtree(node,ancestor,time,spcode))
plot(tree,regimes=bimac["OU.3"],node.names=FALSE,ladderize=TRUE)
plot(tree,regimes=bimac["OU.3"],node.names=TRUE,ladderize=FALSE)

tree1 <- with(bimac,ouchtree(node,ancestor,time,species))
plot(tree1)
plot(tree1,regimes=bimac["OU.LP"],margin=0.7)
try(plot(tree,margin=c(-0.1,0.8),palette=hcl.colors))
try(plot(tree,margin=c(-0.1,0.8,0.5),palette=hcl.colors))
try(plot(tree,margin=c(0.3,0.8),palette=hcl.colors))
try(plot(tree,margin=c(NA,0.8),palette=hcl.colors))
try(plot(tree,margin=2,palette=hcl.colors))
try(plot(tree,margin=-1,palette=hcl.colors))
plot(tree,margin=0.4,palette=hcl.colors(1))
try(plot(tree,margin=0.2,regimes=bimac["OU.3"],palette=hcl.colors(1)))
plot(tree,margin=0.4,regimes=bimac["OU.3"],palette=rainbow(5))
plot(tree,margin=c(0.2,0.4),labels=bimac$species,palette=hcl.colors)

try(plot(tree,regimes=1:3))
try(plot(tree,regimes=bimac$OU.LP))

dev.off()
