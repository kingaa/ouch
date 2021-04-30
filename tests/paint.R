library(ouch)

try(paint(bimac))

with(
  bimac,
  ouchtree(nodes=node,times=time/max(time),ancestors=ancestor,labels=species)
) -> t

h <- hansen(log(bimac['size']),t,bimac['OU.3'],sqrt.alpha=1,sigma=1)

paint(t,subtree=c("1"="medium","9"="large","2"="small"),
  branch=c("38"="large","2"="medium"))
paint(h,subtree=c("1"="medium","9"="large","2"="small"),
  branch=c("38"="large","2"="medium"))
try(
  paint(t,subtree=c("medium","large","small"),
    branch=c("38"="large","2"="medium"))
)
try(
  paint(t,subtree=c(`1`="medium",`9`="large","small"),
    branch=c("38"="large","2"="medium"))
)
try(
  paint(t,subtree=c(`1`="medium",`9`="large",`bob`="small"),
    branch=c("38"="large","2"="medium"))
)
try(
  paint(t,branch=c("large","medium"))
)
paint(t,branch=c('38'="large",'2'="medium"))
try(paint(t,branch=c("large",'2'="medium")))
try(paint(t,branch=c('bob'="large",'2'="medium")))
try(
  paint(h,branch=c("38"="large","2"="medium"),which=3)
)
