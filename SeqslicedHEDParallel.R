# Sliced representative points based on "Squential" Hybrid energy distance criterion

SeqslicedHEDParallel=function(nv,K=length(nv),ts,lambda=NULL,T=200,ifparallel=FALSE)
{ # Confirm the number of slices and the number of points in each slice
  if(length(nv)==1)
  {if(K==1){stop("please specify the number slice")}
    if(nv%%K!=0){stop("If no subset points are specified,nv should be a multiple of K")}
    nv <- rep(nv/K,K)}else{
      if(length(nv)!=K){stop("The number of slices maybe wrong")}}
  if(is.null(lambda)){
    lambda <- rep(1,K) 
    nall <- nv[1]
  for(k in 2:K){ 
    nall<- nall+nv[k]
    lambda[k] <- (nall+nv[k])/(2*nall)  # lambdak=(nc+nk)/(2*nc)
    }}else{    
  if(lambda<0|lambda>1){stop("lambda should in [0,1]")}
      lambda <- rep(lambda,K)}
  print(lambda)
  # Initial points
  N <- nrow(ts)
  p <- ncol(ts)
  trts=t(ts)  #The transposition of training samples is helpful for subsequent calculations
  nv <- sort(nv) # Give priority to generating point sets with small sample sizes
  slicenumber <- rep(1:K,times=nv)# The slice where each point is located(fixed)
  tau <- min(apply(ts,2,sd))/max(N,10000)
  if(ifparallel==TRUE){
    if (!require("foreach")) install.packages("foreach")
    if (!require("doParallel")) install.packages("doParallel")
    cl <- makeCluster(detectCores()-1)
    registerDoParallel(cl)
    cat(paste0("Parallel computing has been enabled, using ", getDoParWorkers(), " worker threads\n"))
  }
  nc=0         # the number of Current points 
  Pnc<- NULL   # Current points
  for(k in 1:K){
    Indexinitial <- sample(1:N,nv[k],replace = FALSE)
    Randomperturbation <- matrix(runif(nv[k]*p,-1,1)*tau,nv[k],p)
    Pnk <- ts[Indexinitial,]+Randomperturbation
    
    #MM algorithm iteration
    Pnknew <- Pnk   
   
    for(t in 1:T){
      if(ifparallel==TRUE){
        Pnknew_list <- foreach(i = 1:nv[k], .combine = rbind) %dopar% {
          xdy <- Pnk[i,]-trts      # Result are saved by row
          normxy <- sqrt(colSums(xdy*xdy))
          normxy[normxy==0] <- 1  # Handle zero norms
          
          xdxk <- t(Pnk[i,]-t(Pnk))
          normxxk <- sqrt(rowSums(xdxk*xdxk))
          normxxk[normxxk==0] <- 1  # Handle zero norms
          xdxk=xdxk/normxxk  # Normalization
          if(nc!=0){
            xdxc <- t(Pnk[i,]-t(Pnc))
            normxxc <- sqrt(rowSums(xdxc*xdxc))
            normxxc[normxxc==0] <- 1  # Handle zero norms
            xdxc <- xdxc/normxxc  # Normalization
            lambdatr <- lambda[k]*(nc)/(nc+nv[k])
          }else{xdxc <- matrix(0,2,2)
          lambdatr <- 0}
          
          # Update point
       
          update <- 1/mean(1/normxy)*(colMeans(ts/normxy)+(1-lambdatr)*colMeans(xdxk)+lambdatr*colMeans(xdxc))
          as.numeric(update)  
        }
        Pnknew <- as.matrix(Pnknew_list)  
      } else {
        for(i in 1:nv[k]){
          xdy <- Pnk[i,]-trts      # Result are saved by row
          normxy <- sqrt(colSums(xdy*xdy))
          normxy[normxy==0] <- 1  # Handle zero norms
          
          xdxk <- t(Pnk[i,]-t(Pnk))
          normxxk <- sqrt(rowSums(xdxk*xdxk))
          normxxk[normxxk==0] <- 1  # Handle zero norms
          xdxk=xdxk/normxxk  # Normalization
          if(nc!=0){
            xdxc <- t(Pnk[i,]-t(Pnc))
            normxxc <- sqrt(rowSums(xdxc*xdxc))
            normxxc[normxxc==0] <- 1  # Handle zero norms
            xdxc <- xdxc/normxxc  # Normalization
            lambdatr <- lambda[k]*(nc)/(nc+nv[k])
          }else{xdxc <- matrix(0,2,2)
          lambdatr <- 0}
          
          # Update point
        
          Pnknew[i,] <- 1/mean(1/normxy)*(colMeans(ts/normxy)+(1-lambdatr)*colMeans(xdxk)+lambdatr*colMeans(xdxc))        
          }
      }
      Pnk <- Pnknew 
  }
 nc <- nc+nv[k]
 Pnc <-rbind(Pnc,Pnk) 
    cat(sprintf("Iteration %d/%d complete\n", k, K))
  }
  
  # Close parallel cluster
  if(ifparallel==TRUE){
    stopCluster(cl)
    registerDoSEQ()  # Restore to serial computing
  }
  
  return(list(Pn=Pnc,slicenumber=slicenumber))
}

library(randtoolbox)

D=SeqslicedHEDParallel(c(10,20,30),ts=sobol(100000,2,scrambling = 3),lambda = NULL,ifparallel = T)
plot(D$Pn[D$slicenumber==1,],col=2,pch=16,cex=2,xlim=c(0,1),ylim=c(0,1))
plot(D$Pn[D$slicenumber==2,],col=3,pch=16,cex=2,xlim=c(0,1),ylim=c(0,1))
plot(D$Pn[D$slicenumber==3,],col=4,pch=16,cex=2,xlim=c(0,1),ylim=c(0,1))
plot(D$Pn,col=D$slicenumber,pch=D$slicenumber+15,cex=2,xlim=c(0,1),ylim=c(0,1))




