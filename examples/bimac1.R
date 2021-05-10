tree <- with(
  bimac,
  ouchtree(nodes=node,ancestors=ancestor,times=time,labels=spcode)
)
tree

plot(tree)
plot(tree, node.names=TRUE)    # display node names
