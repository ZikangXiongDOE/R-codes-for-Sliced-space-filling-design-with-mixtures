source("Gibbssamplingformixtureswithlinearconstraints.R")
source("slicedHEDParallel.R")
source("SeqslicedHEDParallel.R")
source("InverseTransformMethod.R")
source("sliceddesignformixturesbasedonRSPD.R")



library(vcd)
library(SLHD)
nk <- 10
K <- 4
p <- 2
nv <- rep(nk,K)
slicedindex <- rep(1:K,nv)

## ITM
a=rep(0,p+1)
b=rep(1,p+1)

library(vcd)
#ternaryplot(mixlhd,col=rep(1:2,times=nv))

# SLHD package
library(SLHD)
nk <- 10
K <- 4
p <- 2
nv <- rep(nk,K)
slicedindex <- rep(1:K,nv)

## ITM
a=rep(0,p+1)
b=rep(1,p+1)
DSLHD <- maximinSLHD(K,nk,p)
plot(DSLHD$StandDesign[,2:3],col=DSLHD$StandDesign[,1],
     xaxs="i",yaxs="i",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1),pch=16,cex=1.5)
grid(nk,nk,lwd=2,col=4)
ternaryplot(ITM(DSLHD$StandDesign[,2:3],a,b),
            col=DSLHD$StandDesign[,1],pch=DSLHD$StandDesign[,1]+14,main="")
ternaryplot(ITM(DSLHD$StandDesign[slicedindex==1,2:3],a,b),
            col=DSLHD$StandDesign[slicedindex==1,1],
            pch=DSLHD$StandDesign[slicedindex==1,1]+14,main="")
ternaryplot(ITM(DSLHD$StandDesign[slicedindex==2,2:3],a,b),
            col=DSLHD$StandDesign[slicedindex==2,1],
            pch=DSLHD$StandDesign[slicedindex==2,1]+14,main="")
ternaryplot(ITM(DSLHD$StandDesign[slicedindex==3,2:3],a,b),
            col=DSLHD$StandDesign[slicedindex==3,1]
            ,pch=DSLHD$StandDesign[slicedindex==3,1]+14,main="")

## CDM
ddd=SSPDMFB(sum(nv),p+1,K,c(1,1))
table(ddd$Slice)
ternaryplot(ddd$Design,col=ddd$Slice+1,pch=ddd$Slice+15,grid_color = "black",cex=0.7,main="",labels_color = "black")

slicedindex <- rep(1:K,nv)
ts <- Gibbs_l(a,b,n = 100000)

# trainning sample for new methods
ts <- Gibbs_l(a,b,n = 100000)


## ParM lambda=1 means that we only care about uniformity of the overall design
overallD <- slicedHEDParallel(nv = rep(nk,K),ts = ts$sample,lambda = 1,ifparallel = TRUE)
#-- Divide the overall point set using the functions in the R package
library(twinning)
DivideIndex <- multiplet(overallD$Pn,K)
ternaryplot(overallD$Pn,col=DivideIndex+1,pch=DivideIndex+14)
ternaryplot(overallD$Pn[DivideIndex==1,],col=2,pch=15)
ternaryplot(overallD$Pn[DivideIndex==2,],col=3,pch=16)
ternaryplot(overallD$Pn[DivideIndex==3,],col=4,pch=17)

ParMD <- cbind(overallD$Pn,DivideIndex)

## SeqM lambda=1 means that we sequentially generate subdesigns and only care about uniformity of the overall design

SeqMD<- SeqslicedHEDParallel(nv = rep(nk,K),ts = ts$sample,lambda = 1,ifparallel = TRUE)
ternaryplot(SeqMD$Pn,col=SeqMD$slicenumber+1,pch=SeqMD$slicenumber+15)
ternaryplot(SeqMD$Pn[SeqMD$slicenumber==1,],col=1,pch=15,main="")
ternaryplot(SeqMD$Pn[SeqMD$slicenumber==2,],col=2,pch=16,main="")
ternaryplot(SeqMD$Pn[SeqMD$slicenumber==3,],col=3,pch=17,main="")


## MHED

ddnew <- slicedHEDParallel(nv = rep(nk,K),ts = ts$sample,ifparallel = TRUE)
ternaryplot(ddnew$Pn,col=ddnew$slicenumber,pch=ddnew$slicenumber+14,main="")
ternaryplot(ddnew$Pn[ddnew$slicenumber==1,],col=1,pch=15,main="")
ternaryplot(ddnew$Pn[ddnew$slicenumber==2,],col=2,pch=16,main="")
ternaryplot(ddnew$Pn[ddnew$slicenumber==3,],col=3,pch=17,main="")


## SeqHED
Seqddnew <- SeqslicedHEDParallel(nv = rep(nk,K),ts = ts$sample,ifparallel = TRUE)
ternaryplot(Seqddnew$Pn,col=Seqddnew$slicenumber,pch=Seqddnew$slicenumber+15,main="")
ternaryplot(Seqddnew$Pn[Seqddnew$slicenumber==1,],col=1,pch=15,main="")
ternaryplot(Seqddnew$Pn[Seqddnew$slicenumber==2,],col=2,pch=16,main="")
ternaryplot(Seqddnew$Pn[Seqddnew$slicenumber==3,],col=3,pch=17,main="")
