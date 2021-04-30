library(ouch)
set.seed(1106216184L)

png(filename="plot%02d.png")

tree <- with(bimac,ouchtree(node,ancestor,time/max(time),species))
plot(tree,node.names=TRUE)

plot(tree,regimes=bimac["OU.3"],node.names=TRUE,ladderize=TRUE)
plot(tree,regimes=bimac["OU.3"],node.names=FALSE,ladderize=FALSE)

bimac <- bimac[order(runif(n=nrow(bimac))),]
tree <- with(bimac,ouchtree(node,ancestor,time,species))
plot(tree,regimes=bimac["OU.3"],node.names=FALSE,ladderize=TRUE)
plot(tree,regimes=bimac["OU.3"],node.names=TRUE,ladderize=FALSE)

bimac$silly <- ifelse(is.na(bimac$species),NA_character_,"anamagrophilexidicenappendocitemiousness")
tree1 <- with(bimac,ouchtree(node,ancestor,time/max(time),silly))
plot(tree1)
plot(tree1,regimes=bimac["OU.LP"],margin=0.7)
try(plot(tree,margin=c(0.01,0.8),palette=hcl.colors))
try(plot(tree,margin=2,palette=hcl.colors))
try(plot(tree,margin=-1,palette=hcl.colors))
plot(tree,margin=0.2,palette=hcl.colors(1))
try(plot(tree,margin=0.2,regimes=bimac["OU.3"],palette=hcl.colors(1)))
plot(tree,margin=0.2,regimes=bimac["OU.3"],palette=rainbow(5))

dev.off()
