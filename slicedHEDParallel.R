# Sliced representative points based on Hybrid energy distance criterion
# $$E^{\lambda}\left(F,F_{\mathcal{P}_n}\right)$$
# nv: nv=c(n_1,...,n_k),if(length(nv)==1), number of points in each slice is n/K
# K: length(nv) if(length(nv)==1&K==1){please specify the number slice}
# ts: necessary training sample y_1,...y_N
# lambda: balanced parameter
# T: maximum number of iterations
# ifparallel: TRUE means parallel update of the points
slicedHEDParallel=function(nv,K=length(nv),ts,lambda=1/2,T=200,ifparallel=FALSE)
{ # Confirm the number of slices and the number of points in each slice
  if(length(nv)==1)
  {if(K==1){stop("please specify the number slice")}
    if(nv%%K!=0){stop("If no subset points are specified,nv should be a multiple of K")}
    nv <- rep(nv/K,K)}else{
      if(length(nv)!=K){stop("The number of slices maybe wrong")}}
  if(lambda<0|lambda>1){stop("lambda should in [0,1]")}

  # Initial points
  N <- nrow(ts)
  p <- ncol(ts)
  n <- sum(nv)
  Indexinitial <- sample(1:N,n,replace = FALSE)
  tau <- min(apply(ts,2,sd))/max(N,10000)
  
  Randomperturbation <- matrix(runif(n*p,-1,1)*tau,n,p)
  Pn <- ts[Indexinitial,]+Randomperturbation
  slicenumber <- rep(1:K,times=nv)# The slice where each point is located(fixed)
  if(ifparallel==TRUE){
    if (!require("foreach")) install.packages("foreach")
    if (!require("doParallel")) install.packages("doParallel")
    cl <- makeCluster(detectCores()-1) # Adjust the number of threads according to the actual situation
    registerDoParallel(cl)
    cat(paste0("Parallel computing has been enabled, using ", getDoParWorkers(), " worker threads\n"))
  }
  trts=t(ts)  #The transposition of training samples is helpful for subsequent calculations
  
  #MM algorithm iteration
  Pnnew <- Pn
  for(t in 1:T){
    if(ifparallel==TRUE){
      Pnnew_list <- foreach(i = 1:n, .combine = rbind) %dopar% {
        xdy <- Pn[i,]-trts      # Result are saved by row
        normxy <- sqrt(colSums(xdy*xdy))
        normxy[normxy==0] <- 1  # Handle zero norms
        xdx <- t(Pn[i,]-t(Pn))
        normxx <- sqrt(rowSums(xdx*xdx))
        normxx[normxx==0] <- 1  # Handle zero norms
        xdx=xdx/normxx  # Normalization
        
        # Update point
        update <- 1/mean(1/normxy)*(colMeans(ts/normxy)+(1-lambda)*colMeans(xdx[slicenumber==slicenumber[i],])+lambda*colMeans(xdx))
        as.numeric(update)  
      }
      Pnnew <- as.matrix(Pnnew_list)  
    } else {
      
      for(i in 1:n){
        xdy <- Pn[i,]-trts      # Result are saved by row
        normxy <- sqrt(colSums(xdy*xdy))
        normxy[normxy==0] <- 1  # Handle zero norms
        xdx <- t(Pn[i,]-t(Pn))
        normxx <- sqrt(rowSums(xdx*xdx))
        normxx[normxx==0] <- 1  # Handle zero norms
        xdx=xdx/normxx  # Normalization
        
        # Update point
        Pnnew[i,] <- 1/mean(1/normxy)*(colMeans(ts/normxy)+(1-lambda)*colMeans(xdx[slicenumber==slicenumber[i],])+lambda*colMeans(xdx))             
      }
    }
    Pn <- Pnnew 
    cat(sprintf("Iteration %d/%d complete\n", t, T))
  }
  
  # Close parallel cluster
  if(ifparallel==TRUE){
    stopCluster(cl)
    registerDoSEQ()  # Restore to serial computing
  }
  
  return(list(Pn=Pn,slicenumber=slicenumber))
}

library(randtoolbox)

D=slicedHEDParallel(c(10,20,30),ts=sobol(10000,2,scrambling = 3),lambda =1/2,ifparallel = T)
plot(D$Pn[D$slicenumber==1,],col=2,pch=16,cex=2,xlim=c(0,1),ylim=c(0,1))
plot(D$Pn[D$slicenumber==2,],col=3,pch=16,cex=2,xlim=c(0,1),ylim=c(0,1))
plot(D$Pn[D$slicenumber==3,],col=4,pch=16,cex=2,xlim=c(0,1),ylim=c(0,1))
plot(D$Pn,col=D$slicenumber,pch=D$slicenumber+15,cex=2,xlim=c(0,1),ylim=c(0,1))
