library(ouch)
set.seed(1277742405L)

tree <- with(bimac,ouchtree(node,ancestor,time/max(time),species))

try(brown(bimac))
try(brown(tree=bimac))
try(brown(tree=tree))
try(brown(tree=tree,data=bimac$size[1:10]))
try(brown(tree=tree,data=bimac$size))
try(brown(tree=tree,data=setNames(c(bimac$size[-c(44,45)],NA,NA),bimac$node)))
try(
  brown(
    tree=tree,
    data=as.character(bimac$size)
  )
)
brown(
  tree=tree,
  data=list(
    setNames(bimac$size,bimac$node),
    setNames(rnorm(n=45),bimac$node)
  )
) -> m
brown(
  tree=tree,
  data=list(
    A=setNames(bimac$size,bimac$node),
    B=setNames(rnorm(n=45),bimac$node)
  )
) -> m
update(m,data=setNames(2*bimac$size,bimac$node)) -> m
