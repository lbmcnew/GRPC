library(popbio)

#Demographic Parameters
Fy1 = 0.06
Fa1 = 0.06
Sj1 = 0.5
Sy1 = 0.64
Sa1 = 0.71

Fy2 = 0.09
Fa2 = 0.11
Sj2 = 0.5
Sy2 = 0.42
Sa2 = 0.5

Fy3 = 0.21
Fa3 = 0.23
Sj3 = 0.38
Sy3 = 0.34
Sa3 = 0.42

FySj1 = Fy1*Sj1
FaSj1 = Fa1*Sj1

FySj2 = Fy2*Sj2
FaSj2 = Fa2*Sj2

FySj3 = Fy3*Sj3
FaSj3 = Fa3*Sj3

 gpc<-c(
+ FySj1, FaSj1, Sy1, Sa1,
+ FySj2, FaSj2, Sy2, Sa2,
+ FySj3, FaSj3, Sy3, Sa3)

## Covariance matrix
#combine demographic rates into matrix where rows are populations
p1<-matrix(gpc, nrow=3, byrow=TRUE)

# addcolumn names
colnames(p1)<- paste("a", rep(1:2,each=2), 1:2, sep="")

## re-order columns to plot matrix by columns
x<- order(paste("a", 1:2, rep(1:2,each=2), sep=""))
p1<-p1[,x]

covmat<- cov(p1)

## plots matching figure 3 in Brault & Caswell 1993 or figure 10.10 in Caswell 2001 (p 272)
persp(covmat, theta = 60, phi = 15, box = FALSE, main = expression("Matrix Element Covariances GPC"))

# Define a mean matrix
w1<- matrix(apply(p1, 2, mean), nrow=2)

#calculate sensitivities for mean matrix
wS<-sensitivity(w1)

wE<-elasticity(w1)

## V matrix of contributions
contmat <- covmat * c(wS) %*% t(c(wS))

persp(contmat, theta=45, phi=25, box=TRUE, main=expression(paste("Contributions to V(", lambda, ")")) )

## contributions of V associated with aij  (matrix on page 271 in Caswell)
A<-matrix(apply(contmat, 2, mean), nrow=2)

dimnames(A)<-dimnames(whale)

# matrix on page 271
round(A/sum(A),3)



