rm(list=ls())
options(digits=4)
library(popbio)

#Vital Rates
S_vitalrates<- list(Snest = 0.88,
  Srenest = 0.67,
  Stcl1   = 11.8,
  Stcl2   = 11.6,
  Snsurv1y = 0.04,
  Snsurv2y = 0.10,
  Snsurv1a = 0.05,
  Snsurv2a = 0.12,
  Scpere   = 0.77,
  Sbsurv  = 0.29,
  Sfperc   = 0.41,
  Ssj      = 0.5,
  Ssy      = 0.64,
  Ssa      = 0.71); unlist(S_vitalrates)

N_vitalrates<- list(Nnest = 0.81,
  Nrenest  =  0.4,
  Ntcl1    =  11.6,
  Ntcl2    =  10.9,
  Nnsurv1y =  0.08,
  Nnsurv2y =  0.16,
  Nnsurv1a =  0.09,
  Nnsurv2a =  0.18,
  Ncpere   =  0.87,
  Nbsurv  =   0.27,
  Nfperc   =  0.56,
  Nsj      =  0.5,
  Nsy      =  0.42,
  Nsa      =  0.5); unlist(N_vitalrates)
  
SM_vitalrates<- list(SMnest = 0.86,
  SMrenest  =  0.4,
  SMtcl1    =  13.4,
  SMtcl2    =  10.8,
  SMnsurv1y =  0.16,
  SMnsurv2y =  0.28,
  SMnsurv1a =  0.18,
  SMnsurv2a =  0.31,
  SMcpere   =  0.83,
  SMbsurv  =   0.34,
  SMfperc   =  0.51,
  SMsj      =  0.38,
  SMsy      =  0.34,
  SMsa      =  0.42); unlist(SM_vitalrates)

S_elements <- expression(((Snest*Stcl1*Snsurv1y)+((1-Snsurv1y)*Srenest*Stcl2*Snsurv2y))*Scpere*Sbsurv*Sfperc*.5*Ssj,
                   ((Snest*Stcl1*Snsurv1a)+((1-Snsurv1a)*Srenest*Stcl2*Snsurv2a))*Scpere*Sbsurv*Sfperc*.5*Ssj, Ssy, Ssa)
                   
N_elements <- expression(((Nnest*Ntcl1*Nnsurv1y)+((1-Nnsurv1y)*Nrenest*Ntcl2*Nnsurv2y))*Ncpere*Nbsurv*Nfperc*.5*Nsj,
                   ((Nnest*Ntcl1*Nnsurv1a)+((1-Nnsurv1a)*Nrenest*Ntcl2*Nnsurv2a))*Ncpere*Nbsurv*Nfperc*.5*Nsj, Nsy, Nsa)
                   
SM_elements<- expression(((SMnest*SMtcl1*SMnsurv1y)+((1-SMnsurv1y)*SMrenest*SMtcl2*SMnsurv2y))*SMcpere*SMbsurv*SMfperc*.5*SMsj,
                   ((SMnest*SMtcl1*SMnsurv1a)+((1-SMnsurv1a)*SMrenest*SMtcl2*SMnsurv2a))*SMcpere*SMbsurv*SMfperc*.5*SMsj, SMsy, SMsa)
                   

#Substitute vital rates into all of the elements
S_vrs <- sapply(S_elements, eval, S_vitalrates, NULL); S_vrs
N_vrs <- sapply(N_elements, eval, N_vitalrates, NULL); N_vrs
SM_vrs <- sapply(SM_elements, eval, SM_vitalrates, NULL); SM_vrs

#Calculate dimensions of matrix
n <- sqrt(length(S_elements)); n


#Recast vector of elements into a matrix by row
S_mat <- matrix(S_vrs, nrow = n, byrow = TRUE); S_mat
N_mat <- matrix(N_vrs, nrow = n, byrow = TRUE); N_mat
SM_mat <- matrix(SM_vrs, nrow = n, byrow = TRUE); SM_mat

#Eigenvalues and eigenvectors of matrix
S_ev<-eigen(S_mat)
N_ev<-eigen(N_mat)
SM_ev<-eigen(SM_mat)

#Lambda
S_imax <- which.max(S_ev$values); S_imax
S_lambda <- Re(S_ev$values[S_imax]); S_lambda

N_imax <- which.max(N_ev$values); N_imax
N_lambda <- Re(N_ev$values[N_imax]); N_lambda

SM_imax <- which.max(SM_ev$values); SM_imax
SM_lambda <- Re(SM_ev$values[SM_imax]); SM_lambda

# SAD

S_w.tmp <- S_ev$vectors[, S_imax]
S_w <- Re(S_w.tmp/sum(S_w.tmp)); S_w

N_w.tmp <- N_ev$vectors[, N_imax]
N_w <- Re(N_w.tmp/sum(N_w.tmp)); N_w

SM_w.tmp <- SM_ev$vectors[, SM_imax]
SM_w <- Re(SM_w.tmp/sum(SM_w.tmp)); SM_w

# Reprod values
S_v.tmp <- eigen(t(S_mat))$vectors[,S_imax]
S_v <- Re(S_v.tmp/S_v.tmp[1]); S_v

N_v.tmp <- eigen(t(N_mat))$vectors[,N_imax]
N_v <- Re(N_v.tmp/N_v.tmp[1]); N_v

SM_v.tmp <- eigen(t(SM_mat))$vectors[,SM_imax]
SM_v <- Re(SM_v.tmp/SM_v.tmp[1]); SM_v

#Sensitivities
S_sens <- (S_v%*%t(S_w))/as.numeric(S_v %*% S_w); S_sens
N_sens <- (N_v%*%t(N_w))/as.numeric(N_v %*% N_w); N_sens
SM_sens <- (SM_v%*%t(SM_w))/as.numeric(SM_v %*% SM_w); SM_sens

#Elasticities
S_elas <- (S_mat/S_lambda)* S_sens; S_elas
N_elas <- (N_mat/N_lambda)* N_sens; N_elas
SM_elas <- (SM_mat/SM_lambda)* SM_sens; SM_elas

##Calculate the sensitivities, elasticities, and VSSs of Lower Leverl vital rates
#Take the derivative of each elsment with respoect to each lower level vital rates
S_deriv.funcs <- sapply(S_elements, deriv, namevec = names(S_vitalrates), function.arg = TRUE)
N_deriv.funcs <- sapply(N_elements, deriv, namevec = names(N_vitalrates), function.arg = TRUE)
SM_deriv.funcs <- sapply(SM_elements, deriv, namevec = names(SM_vitalrates), function.arg = TRUE)

#Substitute vital rates into the derivative functions
S_devs <- lapply(S_deriv.funcs, function(x) do.call(x,S_vitalrates))
N_devs <- lapply(N_deriv.funcs, function(x) do.call(x,N_vitalrates))
SM_devs <- lapply(SM_deriv.funcs, function(x) do.call(x,SM_vitalrates))

#Create a data frame for teh output for ll rates
S_lower <-data.frame(Estimate = unlist(S_vitalrates), Sensitivity = 0, Elasticity = 0, VSS = 0)
N_lower <-data.frame(Estimate = unlist(N_vitalrates), Sensitivity = 0, Elasticity = 0, VSS = 0)
SM_lower <-data.frame(Estimate = unlist(SM_vitalrates), Sensitivity = 0, Elasticity = 0, VSS = 0)

#For each vital rate, multiply the derivatives by sens
#loop for South population
for (i in 1:length(S_vitalrates)){
  S_derivs <- matrix(as.numeric(lapply(S_devs, function(x)
    attr(x,"gradient")[i])), nrow = n, byrow = TRUE)
  S_lower[i,2] <- sum(S_derivs * S_sens)
  S_lower[i,3] <- S_vitalrates[[i]]/S_lambda*sum(S_derivs*S_sens)
 S_num <- (1-S_vitalrates[[i]])*(S_vitalrates[[i]])
 S_lower[i,4] <- sqrt(S_num)/S_lambda*sum(S_derivs*S_sens)
 }
 S_lower
 
 #loop for North population
for (i in 1:length(N_vitalrates)){
  N_derivs <- matrix(as.numeric(lapply(N_devs, function(x)
    attr(x,"gradient")[i])), nrow = n, byrow = TRUE)
  N_lower[i,2] <- sum(N_derivs * N_sens)
  N_lower[i,3] <- N_vitalrates[[i]]/N_lambda*sum(N_derivs*N_sens)
 N_num <- (1-N_vitalrates[[i]])*(N_vitalrates[[i]])
 N_lower[i,4] <- sqrt(N_num)/N_lambda*sum(N_derivs*N_sens)
 }
 N_lower
 
#loop for Smoky population
for (i in 1:length(SM_vitalrates)){
  SM_derivs <- matrix(as.numeric(lapply(SM_devs, function(x)
    attr(x,"gradient")[i])), nrow = n, byrow = TRUE)
  SM_lower[i,2] <- sum(SM_derivs * SM_sens)
  SM_lower[i,3] <- SM_vitalrates[[i]]/SM_lambda*sum(SM_derivs*SM_sens)
 SM_num <- (1-SM_vitalrates[[i]])*(SM_vitalrates[[i]])
 SM_lower[i,4] <- sqrt(SM_num)/SM_lambda*sum(SM_derivs*SM_sens)
 }
 SM_lower
 

############ Random effects LTRE###################

# enter matrices
South =c(unlist(S_vitalrates))
North =c(unlist(N_vitalrates))
Smoky =c(unlist(SM_vitalrates))

allmat <- rbind(South, North, Smoky)

#Name columns
colnames(allmat) <- c("NEST", "RENEST", "TCL1", "TCL2", "NSURV1y", "NSURV2y", "NSURV1a", "NSURV2a", "CperE", "BSURV", "FperC", "Sj", "Sy", "Sa")

#Calculate covariance matrix
(covmat <- cov(allmat))

#Calculate a Mean Vital Rates across populations
meanvitals <- apply(allmat,2,mean);meanvitals

#Calculate Lower Level Sensitivities for the Mean Matrix

vitalrates<- list(NEST=0.85, RENEST=0.49, TCL1=12.26, TCL2=11.1, NSURV1y=0.093, NSURV2y=0.180, NSURV1a=0.1067, NSURV2a=0.203, CperE=0.823, BSURV=0.300, FperC=0.493, Sj=0.46, Sy=0.467, Sa=0.543); unlist(vitalrates)

elements <- expression(((NEST*TCL1*NSURV1y)+((1-NSURV1y)*RENEST*TCL2*NSURV2y))*CperE*BSURV*FperC*.5*Sj,
                   ((NEST*TCL1*NSURV1a)+((1-NSURV1a)*RENEST*TCL2*NSURV2a))*CperE*BSURV*FperC*.5*Sj, Sy, Sa); elements

#Substitute vital rates into all of the elements
vrs <- sapply(elements, eval, vitalrates, NULL); vrs

#Calculate dimensions of matrix
n <- sqrt(length(elements)); n

#Recast vector of elements into a matrix by row
mat <- matrix(vrs, nrow = n, byrow = TRUE); mat

#Eigenvalues and eigenvectors of matrix
ev<-eigen(mat)

#Lambda
imax <- which.max(ev$values); imax
lambda <- Re(ev$values[imax]); lambda

w.tmp <- ev$vectors[, imax]
w <- Re(w.tmp/sum(w.tmp)); w

# Reprod values
v.tmp <- eigen(t(mat))$vectors[,imax]
v <- Re(v.tmp/v.tmp[1]); v

#Sensitivities
sens <- (v%*%t(w))/as.numeric(v %*% w); sens

#Elasticities
elas <- (mat/lambda)* sens; elas

#Take the derivative of each elsment with respoect to each lower level vital rates
deriv.funcs <- sapply(elements, deriv, namevec = names(vitalrates), function.arg = TRUE)

#Substitute vital rates into the derivative functions
devs <- lapply(deriv.funcs, function(x) do.call(x,vitalrates))

#Create a data frame for teh output for ll rates
lower <-data.frame(Estimate = unlist(vitalrates), Sensitivity = 0, Elasticity = 0, VSS = 0)

#For each vital rate, multiply the derivatives by sens
for (i in 1:length(vitalrates)){
  derivs <- matrix(as.numeric(lapply(devs, function(x)
    attr(x,"gradient")[i])), nrow = n, byrow = TRUE)
  lower[i,2] <- sum(derivs * sens)
  lower[i,3] <- vitalrates[[i]]/lambda*sum(derivs*sens)
 num <- (1-vitalrates[[i]])*(vitalrates[[i]])
 lower[i,4] <- sqrt(num)/lambda*sum(derivs*sens)
 }
 lower

sens<- lower[,2]


# matrix of contributions
# c(sens) converts a matrix to a vector by *columns*
contmat <-  c(sens) %*% t(c(sens)) * covmat

# contributions to Var(lambda) by matrix elements
# use byrow=F because columns were reordered
contrib <- matrix(apply(contmat, 2, sum), nrow=1, byrow=F)
(varlamb <- sum(contrib))
(contrib1 <- contrib/varlamb)

#Variance of the three lambdas (to compare to the approximation above)
alllambda <- c(S_lambda, N_lambda, SM_lambda)
(varlambda<- var(alllambda))

#contributions sum to 0.0169, actual variance is 0.017

write.table(contmat, file = "C:/Users/LB McNew/Documents/GPC/WORD PROCESSING/MANUSCRIPTS/4. Comparative Demography/Analysis/Contmat.txt", sep = " ", quote = FALSE, append = FALSE)


## Surface plots of covariance and contributions
persp(covmat, theta = 45, phi = 25, box = FALSE, main = expression("Matrix Element Covariances GPC"))

persp(contmat, theta=45, phi=25, box=TRUE, main=expression(paste("Contributions to V(", lambda, ")")) )


 
# 2-panel plot, fig. 3 of Brault and Caswell 1993, fig. 10.10 of Caswell 2001:272
par(mfrow=c(1,2), mar=c(1,1,1,1), lwd=2, las=1, cex.lab=2)
persp(covmat, theta=45, phi=18, xlab="covmat by col", ylab="covmat by col", zlab="var-cov values", 
main=expression("Covariances"))
persp(contmat, theta=45, phi=18, xlab="contmat by col", ylab="contmat by col", zlab="contributions", 
main=expression(paste("Contributions to Var(",lambda,")")))
# fig. 3a of Brault and Caswell has an error
# element labels should be a11, a21, a31...a44 for columns (not rows)

# 2-panel plot for matrix image
windows()
par(mfrow=c(1,2))
image2(covmat*1e3, round=1)
image2(contmat*1e4, round=1)
# fig. 3a of Brault and Caswell has an error, peak B is labeled cov(P1,P2)
# fig. 10.10a of Caswell has an error, peak B is labeled cov(G1,G2)
# peak B is cov(G1,P2))