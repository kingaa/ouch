library(ouch)

treedat <- data.frame(node=1:7,anc=NA,time=NA)
treedat$anc[-1] <- c(1,2,2,1,5,5)
treedat$time <- c(0,1,3,3,2,3,3)
treedat$reg <- as.factor(c("a","a","a","b","b","a","b"))

treedat$data1[c(3,4,6,7)] <- c(-0.00133741562273541, -1.01339260146802, -0.0022452348823774, -1.00140366716313)
treedat$data2[c(3,4,6,7)] <- c(0.0161586802663828, 1.01327463666717, 0.0125379175716403, 0.991540176709494)
treedat$data3[c(3,4,6,7)] <- c(-4,-1,0,-2)
treedat$data4[c(3,4,6,7)] <- c(0.5,1,-3,1)
treedat$data5[c(3,4,6,7)] <- c(2,0,1,-2)

tree <- with(treedat,ouchtree(node,anc,time))
x <- treedat[c(6,4,5,7,3,2,1),c("data1","data2")]
btree <- brown(data=x,tree)
print(btree)
htree <- hansen(data=x,tree,regimes=treedat["reg"],sqrt.alpha=c(1,0,1),sigma=c(0.1,0,0.1),fit=F)
print(htree)
htree <- update(htree,sqrt.alpha=c(1,-0.1,1),sigma=c(0.5,0.5,1),fit=F)
print(htree)

htree <- hansen(
                tree=tree,
                data=treedat[c('data3','data4','data5')],
                regimes=treedat['reg'],
                sqrt.alpha=c(1,0.5,0.5,1,0,3),
                sigma=c(0.5,-0.1,-0.25,1,0.5,1),
                fit=F
                )
print(htree)

htree <- hansen(
                tree=tree,
                data=treedat[c('data3','data4','data5')],
                regimes=treedat['reg'],
                sqrt.alpha=c(1,0.5,0.5,1,0,3),
                sigma=c(0.5,-0.1,-0.25,1,0.5,1),
                fit=F
                )
