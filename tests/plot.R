library(ouch)
set.seed(1106216184L)

png(filename="plot%02d.png")

tree <- with(bimac,ouchtree(node,ancestor,time/max(time),species))
plot(tree,node.names=TRUE)

plot(tree,regimes=bimac["OU.3"],node.names=TRUE)

bimac <- bimac[order(runif(n=nrow(bimac))),]
tree <- with(bimac,ouchtree(node,ancestor,time/max(time),species))
plot(tree,regimes=bimac["OU.3"],node.names=TRUE)

dev.off()
