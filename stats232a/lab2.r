#Lab 2
######

#Input vineyard data
####################
vine.mat <- matrix(scan("vineyard.dat"),ncol=4,byrow=T)
y <- vine.mat[,4]
n <- length(y)
library(MASS)


######################################################################
#Degree 2 Fits
######################################################################
deg <- 2
d1 <- deg + 1
x <- 1:n
X <- matrix(rep(0, n*d1),ncol=d1,byrow=T)
for(c in 0:deg){
   X[,c+1] <- x^c
}   

print(paste("Degree of the polynomial fit is:",deg),quote=F)
print(paste(""),quote=F)

#Classical method
#################
eta1 <- X %*% solve(t(X)%*%X) %*% t(X) %*% y

#Pseudoinverse method
#####################
eta2 <- X %*% ginv(X) %*% y

#SVD method
###########
U <- svd(X)$u
eta3 <- U %*% t(U) %*% y

#QR method
##########
Q <- qr.Q(qr(X))
eta4 <- Q %*% t(Q) %*% y

#Orthogonal Polynomial Method
#############################
x0 <- rep(1,n)/sqrt(n)
W <- cbind(x0,poly(x,deg))
eta5 <- W %*% t(W) %*% y

#Plots
######
postscript(file="lab2-1.eps",width=7,height=10,horizontal=F,
pointsize=14)
par(mfcol=c(3,2))

plot(1:n,y,xlab="Vineyard Row",ylab="Yield",type="p",ylim =
c(0,30))
title(main="Classical Degree 2 Fit")
lines(1:n,eta1)

plot(1:n,y,xlab="Vineyard Row",ylab="Yield",type="p",ylim =
c(0,30))
title(main="Pseudoinverse Degree 2 Fit")
lines(1:n,eta2)

plot(1:n,y,xlab="Vineyard Row",ylab="Yield",type="p",ylim =
c(0,30))
title(main="SVD Degree 2 Fit")
lines(1:n,eta3)

plot(1:n,y,xlab="Vineyard Row",ylab="Yield",type="p",ylim =
c(0,30))
title(main="QR Degree 2 Fit")
lines(1:n,eta4)

plot(1:n,y,xlab="Vineyard Row",ylab="Yield",type="p",ylim =
c(0,30))
title(main="Orth. Poly. Degree 2 Fit")
lines(1:n,eta5)

dev.off()

######################################################################
#Degree 6 Fits
######################################################################
deg <- 6
d1 <- deg + 1
x <- 1:n
X <- matrix(rep(0, n*d1),ncol=d1,byrow=T)
for(c in 0:deg){
   X[,c+1] <- x^c
}   

print(paste("Degree of the polynomial fit is:",deg),quote=F)
print(paste(""),quote=F)

#Classical method
#################
#eta1 <- X %*% solve(t(X)%*%X) %*% t(X) %*% y
print(paste("Classical method fails"),quote=F)
print(paste(""),quote=F)

#Pseudoinverse method
#####################
eta2 <- X %*% ginv(X) %*% y

#SVD method
###########
U <- svd(X)$u
eta3 <- U %*% t(U) %*% y

#QR method
##########
Q <- qr.Q(qr(X))
eta4 <- Q %*% t(Q) %*% y

#Orthogonal Polynomial Method
#############################
x0 <- rep(1,n)/sqrt(n)
W <- cbind(x0,poly(x,deg))
eta5 <- W %*% t(W) %*% y

#Plots
######
postscript(file="lab2-2.eps",width=7,height=10,horizontal=F,
pointsize=14)
par(mfcol=c(3,2))

plot(1:n,y,xlab="Vineyard Row",ylab="Yield",type="p",ylim =
c(0,30))
title(main="Classical Degree 6 Fit")
#lines(1:n,eta1)

plot(1:n,y,xlab="Vineyard Row",ylab="Yield",type="p",ylim =
c(0,30))
title(main="Pseudoinverse Degree 6 Fit")
lines(1:n,eta2)

plot(1:n,y,xlab="Vineyard Row",ylab="Yield",type="p",ylim =
c(0,30))
title(main="SVD Degree 6 Fit")
lines(1:n,eta3)

plot(1:n,y,xlab="Vineyard Row",ylab="Yield",type="p",ylim =
c(0,30))
title(main="QR Degree 6 Fit")
lines(1:n,eta4)

plot(1:n,y,xlab="Vineyard Row",ylab="Yield",type="p",ylim =
c(0,30))
title(main="Orth. Poly. Degree 6 Fit")
lines(1:n,eta5)

dev.off()


######################################################################
#Degree 14 Fits
######################################################################
deg <- 14
d1 <- deg + 1
x <- 1:n
X <- matrix(rep(0, n*d1),ncol=d1,byrow=T)
for(c in 0:deg){
   X[,c+1] <- x^c
}   

print(paste("Degree of the polynomial fit is:",deg),quote=F)
print(paste(""),quote=F)

#Classical method
#################
#eta1 <- X %*% solve(t(X)%*%X) %*% t(X) %*% y
print(paste("Classical method fails"),quote=F)
print(paste(""),quote=F)

#Pseudoinverse method
#####################
eta2 <- X %*% ginv(X) %*% y

#SVD method
###########
U <- svd(X)$u
eta3 <- U %*% t(U) %*% y

#QR method
##########
Q <- qr.Q(qr(X))
eta4 <- Q %*% t(Q) %*% y

#Orthogonal Polynomial Method
#############################
x0 <- rep(1,n)/sqrt(n)
W <- cbind(x0,poly(x,deg))
eta5 <- W %*% t(W) %*% y

#Plots
######
postscript(file="lab 2-3.eps",width=7,height=10,horizontal=F,
pointsize=14)
par(mfcol=c(3,2))

plot(1:n,y,xlab="Vineyard Row",ylab="Yield",type="p",ylim =
c(0,30))
title(main="Classical Degree 14 Fit")
#lines(1:n,eta1)

plot(1:n,y,xlab="Vineyard Row",ylab="Yield",type="p",ylim =
c(0,30))
title(main="Pseudoinverse Degree 14 Fit")
lines(1:n,eta2)

plot(1:n,y,xlab="Vineyard Row",ylab="Yield",type="p",ylim =
c(0,30))
title(main="SVD Degree 14 Fit")
lines(1:n,eta3)

plot(1:n,y,xlab="Vineyard Row",ylab="Yield",type="p",ylim =
c(0,30))
title(main="QR Degree 14 Fit")
lines(1:n,eta4)

plot(1:n,y,xlab="Vineyard Row",ylab="Yield",type="p",ylim =
c(0,30))
title(main="Orth. Poly. Degree 14 Fit")
lines(1:n,eta5)

dev.off()




