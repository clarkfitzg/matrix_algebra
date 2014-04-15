#Lab 5
######

#Input dental data as dental(G,C,D)
###################################
y <- scan("dental.dat")
y.arr <- array(y, dim=c(8,3,5))
dimen <- dim(y.arr)
pG <- dimen[1]
pC <- dimen[2]
pD <- dimen[3]
p <- pG*pC*pD

print(paste("Three Nominal Factors: Dental Data"),quote=F)
print(paste(""),quote=F)

print(paste("No. of Gold Alloys is", pG),quote=F)
print(paste("No. of Condensation Methods is", pC),quote=F)
print(paste("No. of Dentists is", pD),quote=F)
print(paste(""),quote=F)

print(paste("Total No. of Factor Level Triples is",p),quote=F)
print(paste(""),quote=F)

#Basic projection function of dimension n 
#########################################
flat <- function(n){
   u <- rep(1,n)/sqrt(n)
   J <- outer(u,u)         
   H <- diag(n) - J
   list(J=J,H=H)
}

#Construct basic factor projections
###################################
flat.G <- flat(pG)
HG <- flat.G$H
JG <- flat.G$J

flat.C <- flat(pC)
HC <- flat.C$H
JC <- flat.C$J


flat.D <- flat(pD)
HD <- flat.D$H
JD <- flat.D$J

#Three-fold Kronecker product
#############################
kprod <- function(A,B,C){
        kprod <- kronecker(B,C)
        kronecker(A,kprod)
        }

#Trace function
###############
tr <- function(A){
   sum(diag(A))      
   }

#Norm squared
#############
normsq <- function(x){
   sum(as.vector(x)^2)
   }

#Construct ANOVA projections
############################
P0 <- kprod(JD,JC,JG)
PG <- kprod(JD,JC,HG)
PC <- kprod(JD,HC,JG)
PD <- kprod(HD,JC,JG)
PGC <- kprod(JD,HC,HG)
PCD <- kprod(HD,HC,JG)
PGD <- kprod(HD,JC,HG)
PGCD <- kprod(HD,HC,HG)

#Pooled variance estimators
###########################
ssGD <- normsq(PGD%*%y)/tr(PGD)
ssGC <- normsq(PGC%*%y)/tr(PGC)
ssCD <- normsq(PCD%*%y)/tr(PCD)
ssGCD <- normsq(PGCD%*%y)/tr(PGCD)

print(paste("Variance estimator GD is", round(ssGD,2)),quote=F)
print(paste("Variance estimator GC is", round(ssGC,2)),quote=F)
print(paste("Variance estimator CD is", round(ssCD,2)),quote=F)
print(paste("Variance estimator GCD is", round(ssGCD,2)),quote=F)
print(paste(""),quote=F)

#Chosen variance estimator
##########################
ss <- ssGD
print(paste("Sigmahat^2 is", round(ss,2)),quote=F)
print(paste(""),quote=F)

print(paste("Estimated LSE risk is", round(ss,2)),quote=F)
print(paste(""),quote=F)

#Adaptive shrinkage constants
#############################
print(paste("Adaptive shrinkage constants"),quote=F)

tauhat0 <- ss*tr(P0)/p
what0 <- normsq(P0%*%y)/p - tauhat0
what0p <- max(what0,0)
what0m <- min(what0,0)
ahat0 <- what0p/(what0p + tauhat0)
print(paste("ahat0 is", round(ahat0,3)),quote=F)

tauhatG <- ss*tr(PG)/p
whatG <- normsq(PG%*%y)/p - tauhatG
whatGp <- max(whatG,0)
whatGm <- min(whatG,0)
ahatG <- whatGp/(whatGp + tauhatG)
print(paste("ahatG is", round(ahatG,3)),quote=F)

tauhatC <- ss*tr(PC)/p
whatC <- normsq(PC%*%y)/p - tauhatC
whatCp <- max(whatC,0)
whatCm <- min(whatC,0)
ahatC <- whatCp/(whatCp + tauhatC)
print(paste("ahatC is", round(ahatC,3)),quote=F)

tauhatD <- ss*tr(PD)/p
whatD <- normsq(PD%*%y)/p - tauhatD
whatDp <- max(whatD,0)
whatDm <- min(whatD,0)
ahatD <- whatDp/(whatDp + tauhatD)
print(paste("ahatD is", round(ahatD,3)),quote=F)

tauhatGC <- ss*tr(PGC)/p
whatGC <- normsq(PGC%*%y)/p - tauhatGC
whatGCp <- max(whatGC,0)
whatGCm <- min(whatGC,0)
ahatGC <- whatGCp/(whatGCp + tauhatGC)
print(paste("ahatGC is", round(ahatGC,3)),quote=F)

tauhatCD <- ss*tr(PCD)/p
whatCD <- normsq(PCD%*%y)/p - tauhatCD
whatCDp <- max(whatCD,0)
whatCDm <- min(whatCD,0)
ahatCD <- whatCDp/(whatCDp + tauhatCD)
print(paste("ahatCD is", round(ahatCD,3)),quote=F)

tauhatGD <- ss*tr(PGD)/p
whatGD <- normsq(PGD%*%y)/p - tauhatGD
whatGDp <- max(whatGD,0)
whatGDm <- min(whatGD,0)
ahatGD <- whatGDp/(whatGDp + tauhatGD)
print(paste("ahatGD is", round(ahatGD,3)),quote=F)

tauhatGCD <- ss*tr(PGCD)/p
whatGCD <- normsq(PGCD%*%y)/p - tauhatGCD
whatGCDp <- max(whatGCD,0)
whatGCDm <- min(whatGCD,0)
ahatGCD <- whatGCDp/(whatGCDp + tauhatGCD)
print(paste("ahatGCD is", round(ahatGCD,3)),quote=F)
print(paste(""),quote=F)

#Adaptive shrinkage estimator
############################
mhat.shr <-
ahat0*P0%*%y+ahatG*PG%*%y+ahatC*PC%*%y+ahatD*PD%*%y+ahatGC*PGC%*%y
mhat.shr <- mhat.shr+ahatCD*PCD%*%y+ahatGD*PGD%*%y+ahatGCD*PGCD%*%y

#Residuals from adaptive shrinkage estimator
############################################
res.shr <- y - mhat.shr

postscript(file="fig1.eps",height=6,width=6,horizontal=F,pointsize=10)
#par(mfrow=c(1,2))
qqnorm(res.shr, main="Adaptive Shrinkage Residuals")
qqline(res.shr)
dev.off()

#Estimated risk of adaptive shrinkage estimator
###############################################
first.shr <-
tauhat0*ahat0+tauhatG*ahatG+tauhatC*ahatC+tauhatD*ahatD+tauhatGC*ahatGC
first.shr <- first.shr+tauhatCD*ahatCD+tauhatGD*ahatGD+tauhatGCD*ahatGCD

second.shr <- what0m+whatGm+whatCm+whatDm+whatGCm+whatCDm+whatGDm+whatGCDm
Rhat.shr <- first.shr + second.shr 

print(paste("Estimated shrinkage risk Rhat.shr is",
round(Rhat.shr,2)),quote=F)
print(paste(""),quote=F)


#Projection function
####################
project <- function(x){
   if(x>1/2) pr <- 1 else pr <- 0
#   return(pr)
} 

#Adaptive projection constants
##############################
print(paste("Adaptive projection constants"),quote=F)

bhat0 <- project(ahat0)
print(paste("bhat0 is", round(bhat0,3)),quote=F)

bhatG <- project(ahatG)
print(paste("bhatG is", round(bhatG,3)),quote=F)

bhatC <- project(ahatC)
print(paste("bhatC is", round(bhatC,3)),quote=F)

bhatD <- project(ahatD)
print(paste("bhatD is", round(bhatD,3)),quote=F)

bhatGC <- project(ahatGC)
print(paste("bhatGC is", round(bhatGC,3)),quote=F)

bhatCD <- project(ahatCD)
print(paste("bhatCD is", round(bhatCD,3)),quote=F)

bhatGD <- project(ahatGD)
print(paste("bhatGD is", round(bhatGD,3)),quote=F)

bhatGCD <- project(ahatGCD)
print(paste("bhatGCD is", round(bhatGCD,3)),quote=F)
print(paste(""),quote=F)

#Adaptive projection estimator
##############################
mhat.pro <-
bhat0*P0%*%y+bhatG*PG%*%y+bhatC*PC%*%y+bhatD*PD%*%y+bhatGC*PGC%*%y
mhat.pro <- mhat.pro+bhatCD*PCD%*%y+bhatGD*PGD%*%y+bhatGCD*PGCD%*%y

#Estimated risk of adaptive projection estimator
################################################
Rhat.pro <- min(tauhat0,what0) + min(tauhatG,whatG) + min(tauhatC,whatC)
Rhat.pro <- Rhat.pro + min(tauhatD,whatD) + min(tauhatGC,whatGC)
Rhat.pro <- Rhat.pro + min(tauhatCD,whatCD) + min(tauhatGD,whatGD)
Rhat.pro <- Rhat.pro + min(tauhatGCD,whatGCD)
print(paste("Estimated projection risk Rhat.pro is",
round(Rhat.pro,2)),quote=F)
print(paste(""),quote=F)

#Two-way marginal fits
######################
mshr.CD <- ahat0*P0%*%y + ahatC*PC%*%y + ahatD*PD%*%y + ahatCD*PCD%*%y
msarr.CD <- array(mshr.CD, dim=c(pG,pC,pD))
ms.CD <- msarr.CD[1,,]

mshr.GC <- ahat0*P0%*%y + ahatG*PG%*%y + ahatC*PC%*%y + ahatGC*PGC%*%y
msarr.GC <- array(mshr.GC, dim=c(pG,pC,pD))
ms.GC <- msarr.GC[,,1]

mshr.GD <- ahat0*P0%*%y + ahatG*PG%*%y + ahatD*PD%*%y + ahatGD*PGD%*%y
msarr.GD <- array(mshr.GD, dim=c(pG,pC,pD))
ms.GD <- msarr.GD[,1,]


mpro.CD <- bhat0*P0%*%y + bhatC*PC%*%y + bhatD*PD%*%y + bhatCD*PCD%*%y
mparr.CD <- array(mpro.CD, dim=c(pG,pC,pD))
mp.CD <- mparr.CD[1,,]

mpro.GC <- bhat0*P0%*%y + bhatG*PG%*%y + bhatC*PC%*%y + bhatGC*PGC%*%y
mparr.GC <- array(mpro.GC, dim=c(pG,pC,pD))
mp.GC <- mparr.GC[,,1]

mpro.GD <- bhat0*P0%*%y + bhatG*PG%*%y + bhatD*PD%*%y + bhatGD*PGD%*%y
mparr.GD <- array(mpro.GD, dim=c(pG,pC,pD))
mp.GD <- mparr.GD[,1,]


#Plots of two-way marginal fits
###############################
postscript(file="fig2.eps",width=7.5,height=10,horizontal=F,
pointsize=10)
par(mfrow=c(3,2))

persp(ms.CD,xlab="Condensation",ylab="Dentist",zlab="Hardness",axes=T,
box=T,zlim=c(500,1000),ticktype="simple",theta=125,phi=30,expand=.5) 
title(main="SHR Fit to CD Means")

persp(mp.CD,xlab="Condensation",ylab="Dentist",zlab="Hardness",axes=T,
box=T, zlim=c(500,1000),ticktype="simple",theta=125,phi=30,expand=.5) 
title(main="PRO Fit to CD Means")

persp(ms.GC, xlab="Gold",ylab="Condensation",zlab="Hardness",axes=T,
box=T, zlim=c(500,1000),ticktype="simple",theta=205,phi=35,expand=.5)  
title(main="SHR Fit to GC Means")

persp(mp.GC, xlab="Gold",ylab="Condensation",zlab="Hardness",axes=T,
box=T, zlim=c(500,1000),ticktype="simple",theta=205,phi=35,expand=.5)  
title(main="PRO Fit to GC Means")

persp(ms.GD, xlab="Gold",ylab="Dentist",zlab="Hardness",axes=T,
box=T, zlim=c(500,1000),ticktype="simple",theta=150,phi=30,expand=.5)  
title(main="SHR Fit to GD Means")

persp(mp.GD, xlab="Gold",ylab="Dentist",zlab="Hardness",axes=T,
box=T, zlim=c(500,1000),ticktype="simple",theta=150,phi=30,expand=.5)  
title(main="PRO Fit to GD Means")

dev.off()


