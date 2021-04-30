library(ouch)
options(digits=3)

try(ouch:::glssoln(0,1,0))
try(ouch:::glssoln(matrix(c(1,1,0,1),2,2),c(1,1),matrix(c(1,1,1,1),2,2)))
ouch:::glssoln(matrix(c(1,1,0,1),2,2),c(1,1),matrix(c(2,1,1,2),2,2))
