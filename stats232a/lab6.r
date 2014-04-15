#Lab 6
######

print(paste("Adaptive PLS Fits to Motorcycle Data"),quote=F)
print(paste(""),quote=F)
library(MASS)

#Input motorcycle data
######################
#y is vector of accelerations
#t.rep is vector of replicated observed times
#t.obs is vector of observed distinct times
#t.all is vector of all distinct times 
#t.miss is vector of times with no data

mot <- matrix(scan("motor.dat"),ncol=2,byrow=T)
t.rep <- mot[,1]
y <- mot[,2]
n <- length(y)
p <- t.rep[n]
t.all <- 1:p

#Construct n x p design matrix X
################################
X <- matrix(rep(0,n*p),nrow=n) 
for(k in 1:n){
        X[k,t.rep[k]]=1
}

#Construct p x 1 replication vector, n x q matrix C, q x p matrix D: X=CD
#########################################################################
rep <- diag(t(X) %*% X)
t.obs <- t.all[rep > .001] 
t.miss <- t.all[rep < .001]
q <- length(t.obs)

print(paste("Row dimension n of X is", n),quote=F)
print(paste("Column dimension p of X is", p),quote=F)
print(paste("Rank of X is", q),quote=F)
print(paste(""),quote=F)

print(paste("The replication vector of length",length(rep),":"),quote=F)
print(rep,quote=F)
print(paste(""),quote=F)

print(paste("The",length(t.obs),"times t_obs with observations:"),quote=F)
print(t.obs,quote=F)
print(paste(""),quote=F)

print(paste("The",length(t.miss),"times t_miss without observations:"),
quote=F)
print(t.miss,quote=F)
print(paste(""),quote=F)

#C <- X[,rep > .001]
#Cplus <- solve(t(C)%*%C)%*%t(C)
#D <- Cplus%*%X

#print(paste("Row dimension n of C is", n),quote=F)
#print(paste("Column dimension q of C is", q),quote=F)
#print(paste("Rank of C is", q),quote=F)
#print(paste(""),quote=F)

#print(paste("Row dimension q of D is", q),quote=F)
#print(paste("Column dimension p of D is", p),quote=F)
#print(paste("Rank of D is", q),quote=F)
#print(paste(""),quote=F)

#Frobenius Norm function
########################
normsq <- function(A){
z <- as.vector(A)
sum(z^2)
} 

#First-difference matrix of dimension (n -1 ) times n
##################################################### 
del <- function(n){
A <- matrix(rep(0, (n-1)*n), ncol=n)
for(k in 1:(n-1)){
        A[k,k] <- -1
        A[k,k+1] <- 1    
    }
   A 
}

#Higher-order difference matrices
#################################
D1 <- del(p)
D2 <- del(p-1)%*%D1
#D3 <- del(p-2)%*%D2
#D4 <- del(p-3)%*%D3

#PLS linear operator
####################
Dif <- D2
XX <- t(X) %*% X
A.pls <- function(lam){
        X %*% ginv(XX + lam*t(Dif)%*%Dif) %*% t(X)    
}

#Minimum relevant eigenvalue
############################
elmin <- function(lam){
        B <- XX + lam*t(Dif)%*%Dif
        el <- eigen(B,symmetric=T,only.values=T)$values
        min(el)        
}

#PLS estimator function (Dmhat = D%*%mhat)
##########################################
fit.pls <- function(lam){
        etahat <- A.pls(lam) %*% y
        mhat <- ginv(XX + lam*t(Dif)%*%Dif) %*% t(X) %*% y       
        list(mhat=mhat,etahat=etahat)   
}

#Estimated risk function
########################
rhat <- function(A){
rh <- normsq(y - A%*%y) + (2*sum(diag(A)) - n)*ss
rh/q
}

risk.pls <- function(lam){
        rhat(A.pls(lam))   
} 


#LS estimator and its risk
##########################
fit.ls <- fit.pls(0)
etahat.ls <- fit.ls$etahat 
ss <- normsq(y - etahat.ls)/(n-q)
risk.ls <- risk.pls(0)

print(paste("LS variance estimate is", round(ss,3)),quote=F)
print(paste("Estimated risk of LS fit is", round(risk.ls,3)),quote=F)
print(paste(""),quote=F)

#Several PLS estimators and their estimated risks
#################################################
risk.vec <- c(risk.pls(1000),risk.pls(2000),risk.pls(3000),risk.pls(4000),
risk.pls(5000))
print(paste("PLS lambda and estimated risks"),quote=F) 
print(c(1000,2000,3000,4000,5000),quote=F)
print(round(risk.vec,3),quote=F)
print(paste(""),quote=F)

#Minimizing estimated PLS risks
###############################
minout <- optimize(risk.pls, interval=c(0,10^6))
lam.opt <- minout$minimum
rk.opt <- minout$objective
fit.opt <- fit.pls(lam.opt)
etahat.opt <- fit.opt$etahat
Dmhat.opt <- fit.opt$Dmhat
mhat.opt <- fit.opt$mhat

print(paste("Minimum estimated PLS risk =",round(rk.opt,3)),quote=F)
print(paste("Optimal PLS lambda =",round(lam.opt,3)),quote=F)
print(paste(""),quote=F)

print(paste("Min eigenvalue of A(lambda.opt) =",round(elmin(lam.opt),4)),
quote=F)
print(paste(""),quote=F)

#Plots
######
res <- y - etahat.opt
mhat.obs <- mhat.opt[t.obs]
mhat.miss <- mhat.opt[t.miss]

postscript(file="fig1.eps", width=6, height=8, horizontal=F, pointsize=10)
par(mfrow=c(2,1))

plot(t.rep, etahat.opt, xlab="Time", ylab="Acceleration",
ylim=c(-140,80), main="PLS Fit and Data", pch="o", cex=.4)
points(t.rep, y, pch="x",cex=.4)

plot(t.obs, mhat.obs, xlab="Time",ylab="Acceleration",
main="Estimated and Extrapolated Means",ylim = c(-140,80),pch="o",cex=.4)
points(t.miss, mhat.miss, pch="x", cex=.4)

dev.off()

postscript(file="fig2.eps", width=6, height=8, horizontal=F, pointsize=10)
par(mfrow=c(2,1))

plot(t.rep,res,xlab="Time",ylab="Residual",
main="Residuals versus Time",pch="o",cex=.4)

#qqnorm(res,cex=.4)
#qqline(res)

dev.off()

