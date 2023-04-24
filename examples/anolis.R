## Analysis of sexual size dimorphism data
\donttest{ ## Save time for CRAN
tree <- with(anolis.ssd,ouchtree(node,ancestor,time/max(time),species))
plot(tree,node.names=TRUE)

h1 <- brown(anolis.ssd['log.SSD'],tree)
h1
plot(h1)

h2 <- hansen(anolis.ssd['log.SSD'],tree,anolis.ssd['OU.1'],sqrt.alpha=1,sigma=1)
h2
plot(h2)

h3 <- hansen(anolis.ssd['log.SSD'],tree,anolis.ssd['OU.7'],sqrt.alpha=1,sigma=1)
h3
plot(h3)
}
