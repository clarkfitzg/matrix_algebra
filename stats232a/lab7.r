#Lab 7
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
  sum(x^2)
}

#Trace
######
tr <- function(A){
   sum(diag(A))
}

#Argmin
#######
argmin <- function(x){
      names(x) <- 1:length(x)
      label <- (names(x[x == min(x)]))
      as.integer(label)
      }

#Competing fits in Hypercube setup
##################################
Cplus <- solve(t(C)%*%C) %*% t(C)
H <- t(C)%*%C
rep <- diag(H)
reproot <- sqrt(rep)
Hroot <- diag(reproot)
Hrootinv <- diag(1/reproot)
U <- C%*%Hrootinv
eye <- diag(p)
N <- H - eye
z <- t(U)%*%y

Q <- function(t){
   t[1]*P0 + t[2]*P1 + t[3]*P2 + t[4]*P12
}

#Canonical shrinkage matrix
###########################
SS <- function(d){
   cen <- solve(Q(d)%*%N%*%Q(d) + eye)
   Hroot%*%Q(d)%*%cen%*%Q(d)%*%Hroot
}
#Hypercube estimator function
#############################
fit <- function(d){
   eta <- U %*% SS(d) %*% z
   m <- Cplus %*% eta 
   M <- array(m,dim=c(4,4))
   list(eta=eta,m=m,M=M)
}

#Estimated risk function
########################
risk <- function(d){
   (normsq(z - SS(d)%*%z) + (2*tr(SS(d))-p)*ss.ls)/p
}   

#Competing fits
###############

#General Model LS fit 
#####################
d <- c(1,1,1,1)
fit.ls <- fit(d)
eta.ls <- fit.ls$eta
m.ls <- fit.ls$m
M.ls <- fit.ls$M
ss.ls <- normsq(y - eta.ls)/(n-p)

print(paste("General model estimated means are:"),quote=F)
print(round(M.ls,2),quote=F)
print(paste(""),quote=F)
print(paste("Estimated sigma^2 is", round(ss.ls,3)),quote=F)
print(paste(""),quote=F)

print(paste("General model estimated risk is", round(risk(d),4)),quote=F)
print(paste(""),quote=F)


#ANOVA Submodel LS fits
#######################
rk.mat <- NULL
for(d1 in 0:1){
  for(d2 in 0:1){
    for(d3 in 0:1){
      for(d4 in 0:1){ 
        d <- c(d1,d2,d3,d4)      
        rk.mat <- rbind(rk.mat,c(risk(d),d))   
      }
    }
  }
}

print(paste("ANOVA submodels risk matrix is:"),quote=F)
print(round(rk.mat,4),quote=F)
print(paste(""),quote=F)

ind.opt <- argmin(rk.mat[,1])
row.opt <- rk.mat[ind.opt,]
d.sub <- c(as.integer(row.opt[2:5]))
rk.sub <- row.opt[1]

print(paste("Best ANOVA submodel est'd risk =",round(rk.sub,5)),quote=F)
print(paste("Best ANOVA submodel d is:"),quote=F)
print(round(d.sub,5),quote=F)
print(paste(""),quote=F)

fit.sub <- fit(d.sub)
eta.sub <- fit.sub$eta
m.sub <- fit.sub$m
M.sub <- fit.sub$M

print(paste("Optimal ANOVA submodel means are:"),quote=F)
print(round(M.sub,2),quote=F)
print(paste(""),quote=F)


#Risk minimization
##################
d.start <- d.sub
minout <- nlm(risk,d.start)
d.opt <- minout$estimate
rk.opt <- risk(d.opt)

print(paste("Optimal estimated risk =",round(rk.opt,5)),quote=F)
print(paste("Optimal d is"),quote=F)
print(round(d.opt,5),quote=F)
print(paste(""),quote=F)

#Best hypercube fit
S.opt <- SS(d.opt)
eta.opt <- U%*%S.opt%*%z
m.opt <- Cplus%*%eta.opt
M.opt <- array(m.opt,dim=c(4,4))

print(paste("Optimal HPLS means are:"),quote=F)
print(round(M.opt,2),quote=F)
print(paste(""),quote=F) 


#Mean and residual plots for full LS and HPLS fits
##################################################
postscript(file="lab7.eps",height=9,width=6.5,horizontal=F,pointsize=12)
par(mfcol=c(3,2))

plot(1:p1,M.ls[1,],xlab="Rat Litter Genotype",
ylab="Mean Litter Weight",ylim = c(45,65),type="n",cex=.5)
title(main="General Model Least Squares Fit")
points(1:p1,M.ls[1,],cex=.5)
lines(1:p1,M.ls[1,])

points(1:p1,M.ls[2,],cex=.5)
lines(1:p1,M.ls[2,],lty=2)

points(1:p1,M.ls[3,],cex=.5)
lines(1:p1,M.ls[3,],lty=3)

points(1:p1,M.ls[4,],cex=.5)
lines(1:p1,M.ls[4,],lty=4)


plot(1:p1,M.sub[1,],xlab="Rat Litter Genotype",
ylab="Mean Litter Weight",ylim = c(45,65),type="n")
title(main="Optimal ANOVA Submodel Fit")
points(1:p1,M.sub[1,],cex=.5)
lines(1:p1,M.sub[1,])

points(1:p1,M.sub[2,],cex=.5)
lines(1:p1,M.sub[2,],lty=2)

points(1:p1,M.sub[3,],cex=.5)
lines(1:p1,M.sub[3,],lty=3)

points(1:p1,M.sub[4,],cex=.5)
lines(1:p1,M.sub[4,],lty=4)


plot(1:p1,M.opt[1,],xlab="Rat Litter Genotype",
ylab="Mean Litter Weight",ylim = c(45,65),type="n")
title(main="Optimal HPLS Fit")
points(1:p1,M.opt[1,],cex=.5)
lines(1:p1,M.opt[1,])

points(1:p1,M.opt[2,],cex=.5)
lines(1:p1,M.opt[2,],lty=2)

points(1:p1,M.opt[3,],cex=.5)
lines(1:p1,M.opt[3,],lty=3)

points(1:p1,M.opt[4,],cex=.5)
lines(1:p1,M.opt[4,],lty=4)

dev.off()

