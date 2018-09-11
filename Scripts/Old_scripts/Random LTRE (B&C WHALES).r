#  bookkeeping
rm(list=ls())
options(digits=4)
library(popbio)

 pods<-c(
+ 0,0.0067,0.1632,0,0.9535,0.8827,0,0,0,0.0802,0.9586,0,0,0,0.0414,0.9752,
+ 0,0.0062,0.1737,0,1,0.9020,0,0,0,0.0694,0.9582,0,0,0,0.0418,0.9855,
+ 0,0.0037,0.0988,0,0.9562,0.9030,0,0,0,0.0722,0.9530,0,0,0,0.0406,0.9798,
+ 0,0.0043,0.1148,0,1,0.9015,0,0,0,0.0727,0.9515,0,0,0,0.0485,0.9667,
+ 0,0.0042,0.1054,0,0.8165,0.8903,0,0,0,0.0774,0.9515,0,0,0,0.0485,0.9810,
+ 0,0.0027,0.0732,0,1,0.9123,0,0,0,0.0730,0.9515,0,0,0,0.0485,0.9545,
+ 0,0.0025,0.0651,0,1,0.9254,0,0,0,0.0746,0.9515,0,0,0,0.0485,0.9810,
+ 0,0.0047,0.1159,0,1,0.9200,0,0,0,0.0800,0.9706,0,0,0,0.0294,0.9608,
+ 0,0.0068,0.1761,0,1,0.9241,0,0,0,0.0759,0.9562,0,0,0,0.0438,1,
+ 0,0.0061,0.1418,0,1,0.9167,0,0,0,0.0833,0.9286,0,0,0,0.0714,1,
+ 0,0.0050,0.1251,0,1,0.9216,0,0,0,0.0784,0.9515,0,0,0,0.0485,0.9810,
+ 0,0.0021,0.0542,0,1,0.9254,0,0,0,0.0746,0.9515,0,0,0,0.0485,0.9810,
+ 0,0.0027,0.0732,0,1,0.9286,0,0,0,0.0714,0.9515,0,0,0,0.0485,0.9810,
+ 0,0.0045,0.1220,0,1,0.9286,0,0,0,0.0714,0.9515,0,0,0,0.0485,1,
+ 0,0.0052,0.1428,0,1,0.9286,0,0,0,0.0714,0.9515,0,0,0,0.0485,0.9810,
+ 0,0.0037,0.0998,0,1,0.9286,0,0,0,0.0714,0.9515,0,0,0,0.0485,0.9810,
+ 0,0.0047,0.1273,0,1,0.9286,0,0,0,0.0714,0.9515,0,0,0,0.0485,0.9810,
+ 0,0.0024,0.0797,0,1,0.8929,0,0,0,0.0595,0.9515,0,0,0,0.0485,1)

## Covariance matrix
p1<-matrix(pods, nrow=18, byrow=TRUE)

# addcolumn names
colnames(p1)<- paste("a", rep(1:4,each=4), 1:4, sep="")

## re-order columns to plot matrix by columns
x<- order(paste("a", 1:4, rep(1:4,each=4), sep=""))
p1<-p1[,x]

covmat<- cov(p1)

## plots matching figure 3 in Brault & Caswell 1993 or figure 10.10 in Caswell 2001 (p 272)
persp(covmat, theta = 45, phi = 15, box = FALSE, main = expression("Fig 10.10a Killer Whale Covariances"))

w1<- matrix(apply(p1, 2, mean), nrow=4)

wS<-sensitivity(w1)

## V matrix of contributions
contmat <- covmat * c(wS) %*% t(c(wS))

persp(contmat, theta=45, phi=15, box=FALSE, main=expression(paste("Fig 10.10b Contributions to V(", lambda, ")")) )

## contributions of V associated with aij  (matrix on page 271 in Caswell)
A<-matrix(apply(contmat, 2, mean), nrow=4)

dimnames(A)<-dimnames(whale)

# matrix on page 271
round(A/sum(A),3)
           yearling juvenile mature postreprod
yearling      0.000    0.029  0.705          0
juvenile      0.032    0.092  0.000          0
mature        0.000    0.071  0.070          0
postreprod    0.000    0.000  0.000          0


