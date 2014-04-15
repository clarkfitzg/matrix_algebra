#Lab 3
######

#Read in rat litter data 
########################
p1 <- 4
p2 <- 4
data <- matrix(scan("litter.dat"),byrow=T,ncol=3)
y <- data[,1]
n <- length(y)
p <- p1*p2

#Incidence matrix
#################
num <- NULL
for(j in 1:p2){
  for(i in 1:p1){
    tem <- y[data[,2]==i & data[,3]==j]
    num <- c(num, length(tem))
  }
}

C <- matrix(rep(0,n*p), nrow=n)
count <- 0
for(j in 1:p){
  for(k in 1:num[j]){
    count <- count + 1
    C[count,j] <- 1
  }
}

#Projection function of dimension n 
###################################
flat <- function(n){
   u <- rep(1,n)/sqrt(n)
   J <- outer(u,u)         
   H <- diag(n) - J
#return(J,H)
list(J=J,H=H)
}

#Two-way projections
####################
flat.1 <- flat(p1)
H1 <- flat.1$H
J1 <- flat.1$J

flat.2 <- flat(p2)
H2 <- flat.2$H
J2 <- flat.2$J

P0 <- kronecker(J2,J1)
P1 <- kronecker(J2,H1)
P2 <- kronecker(H2,J1)
P12 <- kronecker(H2,H1)

#Norm squared
#############
normsq <- function(x){
  y <- sum(x^2)
  return(y)
}

#Invoke Moore-Penrose inverse
#############################
library(MASS)

#Competing fits
###############

#Full model
D <- solve(t(C)%*%C) %*% t(C)
mhat <- D %*% y
etahat <- C %*% mhat
df.ss <- n-p
sshat <- normsq(y - etahat)/df.ss
print(paste("Full model estimated means are:"))
print(round(mhat,3))
print(paste("Estimated sigma^2 is", round(sshat,3)))
print(paste(""))

#Model AB
C.AB <- C %*% (P0 + P1 + P2)
mhat.AB <- ginv(C.AB) %*% y
etahat.AB <- C %*% mhat.AB
df.AB <- (p1-1)*(p2-1)
top.AB <- (normsq(etahat - etahat.AB))/df.AB
T.AB <- top.AB/sshat
pval.AB <- 1 - pf(T.AB,df.AB,df.ss)
print(paste("Model AB p-value is", round(pval.AB,4)))
print(paste(""))

#Model A
C.A <- C %*% (P0 + P1)
mhat.A <- ginv(C.A) %*% y
etahat.A <- C %*% mhat.A
df.A <- p1*p2 - p1
top.A <- (normsq(etahat - etahat.A))/df.A
T.A <- top.A/sshat
pval.A <- 1 - pf(T.A,df.A,df.ss)
print(paste("Model A p-value is", round(pval.A,4)))
print(paste(""))

#Model B
C.B <- C %*% (P0 + P2)
mhat.B <- ginv(C.B) %*% y
etahat.B <- C %*% mhat.B
df.B <- p1*p2 - p2
top.B <- (normsq(etahat - etahat.B))/df.B
T.B <- top.B/sshat
pval.B <- 1 - pf(T.B,df.B,df.ss)
print(paste("Model B p-value is", round(pval.B,4)))
print(paste(""))

#Model 0
C.0 <- C %*% P0
mhat.0 <- ginv(C.0) %*% y
etahat.0 <- C %*% mhat.0
df.0 <- p1*p2 - 1
top.0 <- (normsq(etahat - etahat.0))/df.0
T.0 <- top.0/sshat
pval.0 <- 1 - pf(T.0,df.0,df.ss)
print(paste("Model 0 p-value is", round(pval.0,4)))
print(paste(""))


#Estimated risk function
########################
rhat <- function(eta,moddim){
  rhat <- (normsq(y -eta) + (2*moddim - n)*sshat)/p 
  return(rhat)
}


#Estimated risk calculations
############################

#Full model
Rhat <- rhat(etahat,p)
print(paste("Full model estimated risk is", round(Rhat,4)))
print(paste(""))

#Model AB
Rhat.AB <- rhat(etahat.AB,p1+p2-1)
print(paste("Model AB estimated risk is", round(Rhat.AB,4)))
print(paste(""))

#Model A
Rhat.A <- rhat(etahat.A,p1)
print(paste("Model A estimated risk is", round(Rhat.A,4)))
print(paste(""))

#Model B
Rhat.B <- rhat(etahat.B,p2)
print(paste("Model B estimated risk is", round(Rhat.B,4)))
print(paste(""))

#Model 0
Rhat.0 <- rhat(etahat.0,1)
print(paste("Model 0 estimated risk is", round(Rhat.0,4)))
print(paste(""))

#Residual plot
##############
Residuals <- y - C %*% mhat
postscript(file="lab3-1.eps",height=9,width=6.5,horizontal=F,pointsize=12)
par(mfrow=c(2,1))
qqnorm(Residuals, main="Residuals from Full Model Fit")
qqline(Residuals)
dev.off()
