tree <- with(
  bimac,
  ouchtree(nodes=node,ancestors=ancestor,times=time,labels=spcode)
)
tree

plot(tree)
plot(tree, node.names=TRUE)    # display node names

## When taxon names are long, they are cut off when the
## default settings are used.  For example:
tree2 <- with(
  bimac,
  ouchtree(nodes=node,ancestors=ancestor,times=time,
    labels=ifelse(is.na(species),NA,paste(species,island,sep=", "))
  )
)
plot(tree2) # long species names are cut off
## This is fixed by increasing right margin and font size:
plot(tree2,margin=0.5)
