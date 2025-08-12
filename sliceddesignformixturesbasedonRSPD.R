# R packages "LatticeDesign" and "vcd" are needed for Visualization.
print("SSPDMFB: Mixutre Design Based on the Densest Packed Grid Point Set ")
print(" Parameters in SSPDMFB: total sample size, dimension, number of pieces (default equal to dimension), coset decomposition weight (such as c(1,2) when p=3),
      maximum search count ")
SSPDMFB=function(n,p,slice=p,coast,K=100)
{
  # Part I
  #______________________________________________________
  # Orthogonal transformation Q by QR fractorization, B=QR
  
  B=rbind(rep(-1,p-1),diag(rep(1,p-1)))        # A set of basis vectors of hyperplane in which the T_p is
  qrf=qr(B)
  Q=qr.Q(qrf)                                  # Q is a p*(p-1) matrix
  #######################################################
  
  # Generator matrix G (Notice: G is a (p-1)*(p-1) matrix)
  # Scaling parameter l
  # Level parameter m
  G=sqrt(p/(p-1))*diag(rep(1,p-1))-1/((p-1)^(1/2)*(p^(1/2)-1))*matrix(1,p-1,p-1)
  l=(n*p^((p-2)/2)*(p-1)^(-(p-1)/2)/(sqrt(p)/factorial(p-1)))^(1/(p-1))
  m=ceiling((l*sqrt((p-1)/p)+sqrt((p+1)/12))/(p/2/(p-1)))
  #######################################################
  #______________________________________________________
  
  # Part II
  #______________________________________________________
  # Cartesian product function for generating full factorial array
  Cartesianset=function(X) 
  {
    A=X[1,]
    for(i in 2:nrow(X)){
      A=outer(A,X[i,],paste)
    }
    B=strsplit(A,split = " ")
    C=unlist(B)
    D=as.numeric(C)
    return(matrix(D,ncol = nrow(X),byrow = T))
  }
  ########################################################
  
  # Full factorial array F and Large design E
  F=Cartesianset(matrix(rep(seq(-m,m,1),p-1),nrow=p-1,byrow = T))
  
  # Screen out points by a ball of radius "sqrt((p-1)/p)+rho_c/l" for the first time #
  F=F[sqrt(rowSums((F%*%G)^2))<=l*sqrt((p-1)/p)+sqrt((p+1)/12),]
  
  E=1/l*F%*%G
  
  ########################################################
  #_______________________________________________________
  
  # Part III
  #______________________________________________________
  # Transform lattice back into hyperplane {(x_1,...,x_p)|x_1+...x_p=1} and judge which points are in T_p
  
  ETp=t(Q%*%t(E)+rep(1,p)/p)                       # Transformations
  judge=ETp[,1]>=0                                 # Judge 
  for(i in 2:p)
  {
    judge=judge&(ETp[,i]>=0)
  }
  Delta=abs(sum(judge)-n)                           # Whether is the number of points in T_p is n
  #if(Delta==0)
  #{
  #   Design=ETp[judge,]/apply(ETp[judge,],1,sum)   # Normlize to 1 (The problem of calculation accuracy)
  # return(list(Design=Design,Delta=Delta, FULL=F[judge,],Searchnumber=0))
  #}else
  {for(k in 1:K)
    
    # Part IV
    #______________________________________________________
    # Rotation matrix Rstar and perturbation vector delta are used to
    # search design which exactly or nearly has n points  
  {
    # \delta uniformly distributes in (p-1)-dimensional ball of radius rho_c/l=((p-1+2)/12)^(1/2)/l
    library(MASS)
    unorm=mvrnorm(1,rep(0,p-1),diag(rep(1,p-1)))  # Standard Normal r.v. in p-1 dimension
    usphere=unorm/sqrt(sum(unorm^2))              # Uniform r.v. in p-1 sphere
    ur=(runif(1))^(1/(p-1))                       # Inverse of distribution function of radius (Fang and Wang,1994 Section4.2)
    delta=((p-1+2)/12)^(1/2)/l*ur*usphere         
    #######################################################
    
    # Rotation matrix function for m dimension
    rotate=function(m)
    {             
      CPair = matrix(0,m*(m-1)/2,2)
      row = 1
      for(i in 1:(m-1)) for(j in (i+1):m) { CPair[row,1] <- i; CPair[row,2] <- j; row <- row+1; }
      R = diag(m)
      for(a in 1:((m*(m-1)/2)))
      {
        alpha = runif(1,0,2*pi)
        thepair = CPair[a,]
        W = diag(m)
        W[thepair[1],thepair[1]] = cos(alpha); W[thepair[2],thepair[2]] <- cos(alpha); 
        W[thepair[1],thepair[2]] = sin(alpha); W[thepair[2],thepair[1]] <- -sin(alpha); 
        R = R%*%W
      }
      return(R)
    }
    #######################################################
    
    # The lattice after the perturbation and rotation, and judge which points are in T_p
    Rstar=rotate(p-1)
    Enew=Rstar%*%t(E)+delta
    EnewTp=t(Q%*%Enew+rep(1,p)/p)                    # Transformations
    
    judgenew=EnewTp[,1]>=0                                 # Judge 
    for(i in 2:p)
    {
      judgenew=judgenew&(EnewTp[,i]>=0)
    }
    #######################################################
    
    # Whether is the number of points in T_p is exactly n or nearer than before 
    Deltanew=abs(sum(judgenew)-n)                           
    if(Deltanew==0&k>=5&k<=K)
    {
      Design=EnewTp[judgenew,]/apply(EnewTp[judgenew,],1,sum)      # Normlize to 1 (The problem of calculation accuracy)
      Slice=rowSums(t(t(F[judgenew,])*coast))%%slice
      return(list(Design=Design,Slice= Slice,Delta=Deltanew,FULL=F[judgenew,],coast=coast, Searchnumber=k))
    }else
      if(Deltanew<Delta&k>=5&k<=K)
      {
        Delta=Deltanew
        Design=EnewTp[judgenew,]/apply(EnewTp[judgenew,],1,sum)
        Slice=rowSums(t(t(F[judgenew,])*coast))%%slice
      }
  }
    #######################################################   
    return(list(Design=Design,Slice= Slice,Delta=Delta,coast=coast, Searchnumber=k))}
  #__________________________________________________________
}
# library(vcd)
# ddd=SSPDMFB(140,3,7,c(1,4))
# table(ddd$Slice)
# ternaryplot(ddd$Design,col=ddd$Slice+1,pch=ddd$Slice+11,grid_color = "black",cex=0.7,main="",labels_color = "black")
# 
# sl7=apply(t(t(ddd$FULL)*c(1,4)),1,sum)%%7##(1,1) (1,2) (1,6)
# ternaryplot(ddd$Design,col=sl7+1,pch=sl7+11,grid_color = "black",cex=1,main="",labels_color = "black")

