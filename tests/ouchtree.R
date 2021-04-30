library(ouch)

x <- anolis.ssd
x$node[15] <- x$node[12]
try(with(x,ouchtree(node,ancestor,time,species)))

x <- anolis.ssd
try(with(x,ouchtree(node,ancestor[1:10],time,species)))
try(with(x,ouchtree(node,ancestor,time[1:20],species)))
try(with(x,ouchtree(node,ancestor,time,species[1:20])))

x$time[1] <- 10
try(with(x,ouchtree(node,ancestor,time,species)))

x <- anolis.ssd
x$ancestor[15] <- 15
try(with(x,ouchtree(node,ancestor,time,species)))

x <- anolis.ssd
x$ancestor[15] <- NA
try(with(x,ouchtree(node,ancestor,time,species)))

x <- anolis.ssd
x$ancestor[3] <- 15
try(with(x,ouchtree(node,ancestor,time,species)))

x <- anolis.ssd
x$ancestor[15] <- 888
try(with(x,ouchtree(node,ancestor,time,species)))

x <- anolis.ssd
x$ancestor[15:18] <- 888
try(with(x,ouchtree(node,ancestor,time,species)))
