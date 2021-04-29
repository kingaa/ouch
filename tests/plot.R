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
tree <- with(bimac,ouchtree(node,ancestor,time/max(time),silly))
plot(tree)
plot(tree,margin=c(0.01,0.8),palette=hcl.colors)

dev.off()
