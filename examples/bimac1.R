tree <- with(
  bimac,
  ouchtree(nodes=node,ancestors=ancestor,times=time,labels=species)
)
tree

plot(tree)
plot(tree,node.names=TRUE)
