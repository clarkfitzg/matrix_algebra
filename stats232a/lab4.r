#Lab 4
######

#Input monkey data
#File monkey0.dat is monkey.dat with text header line removed
##############################################################
data <- matrix(scan("monkey0.dat"),ncol=4,byrow=T)
fac1 <- data[,2]
fac2 <- data[,3]
fac3 <- data[,1]
y <- data[,4]
n <- length(y)

print(paste("Latin Square: Monkey Data"))
print(paste(""))

print(paste("No of Factor Levels in each Factor is 5"))
print(paste(""))

#Construct means deletion matrix
################################
ind <- fac1 + 5*(fac2-1) + 25*(fac3-1)
D <- matrix(rep(0,25*125),nrow=25)
for(i in 1:25){
   D[i,ind[i]] <- 1
}

#Univariate projection function
###############################
flat <- function(n){
   u <- rep(1,n)/sqrt(n)
   J <- outer(u,u)         
   H <- diag(n) - J
   list(J=J,H=H)
}

#Construct projection factors
#############################
flat.out <- flat(5)
H <- flat.out$H
J <- flat.out$J

#Three-fold Kronecker product
#############################
kprod <- function(A,B,C){
        kprod <- kronecker(B,C)
        kprod <- kronecker(A,kprod)
}

#Construct projections
######################
P0 <- kprod(J,J,J)
P1 <- kprod(J,J,H)
P2 <- kprod(J,H,J)
P3 <- kprod(H,J,J)

#Norm squared
#############
normsq <- function(x){
  y <- sum(x^2)
}

#Moore-Penrose inverse and rank function, using the svd
####################################################### 
pinv <- function(A){
  svdA <- svd(A)
  U <- svdA$u
  V <- svdA$v
  d <- svdA$d
  dr <- d[d>.000001]
  r <- length(dr)
  Ur <- as.matrix(U[,1:r])
  Vr <- as.matrix(V[,1:r])
  Linv <- matrix(rep(0,r*r),nrow=r)
  diag(Linv) <- 1/dr  
  inv <- Vr %*% Linv %*% t(Ur)
  list(inv=inv,r=r)
}

#Submodel fitting function
#########################
subfit <- function(P){
DP <- D %*% P 
mpinv <- pinv(DP)
etahat <- D %*% mpinv$inv %*% y
r <- mpinv$r
etahat <- as.vector(etahat)
list(etahat=etahat,r=r)
}

#LSE's under additive model
###########################
fit.123 <- subfit(P0 + P1 + P2 + P3)
etahat.123 <- fit.123$etahat
r.123 <- fit.123$r
rm(fit.123)

sshat <- normsq(y - etahat.123)/(n-r.123)

print(paste("Estimated additive model eta is:"),quote=F)
print(round(etahat.123,4))
print("",quote=F)


print(paste("Estimated additive model variance is",
round(sshat,4)),quote=F)
print("",quote=F)
df.123 <- n - r.123
print(paste("Degrees of freedom are", df.123),quote=F)
print("",quote=F)

#Residual plot
##############
Residuals <- y - etahat.123
postscript(file="resid.eps",height=9,width=6.5,horizontal=F,pointsize=12)
par(mfrow=c(2,1))
qqnorm(Residuals, main="Residuals from Additive Model Fit")
qqline(Residuals)
dev.off()

#Estimability
#############
P <- P0 + P1 + P2 + P3
DP <- D%*%P
DPplus <- pinv(DP)$inv
A <- DPplus %*% DP

dif <- rep(1,4)
dif[1] <- normsq(P0%*%A -P0)
dif[2] <- normsq(P1%*%A -P1)
dif[3] <- normsq(P2%*%A -P2)
dif[4] <- normsq(P3%*%A -P3)
print(paste("Estimability check for P_0m, P_1m, P_2m, P_3m:",
round(dif[1],4), round(dif[2],4),round(dif[3],4),round(dif[3],4)),quote=F)
print("", quote=F)

#Least squares estimates of parametric functions
################################################
m0 <- P0 %*% DPplus %*% y
m1 <- P1 %*% DPplus %*% y
m2 <- P2 %*% DPplus %*% y
m3 <- P3 %*% DPplus %*% y

print(paste("muhat is", round(m0[1],4)),quote=F)
print("", quote=F)

print(paste("alphahat is"),quote=F)
print(round(m1[1:5],4),quote=F)
print("", quote=F)

print(paste("betahat is"),quote=F)
print(round(c(m2[1],m2[6],m2[11],m2[16],m2[21]),4),quote=F)
print("", quote=F)

print(paste("gammahat is"),quote=F)
print(round(c(m3[1],m3[26],m3[51],m3[76],m3[101]),4),quote=F)
print("", quote=F)
