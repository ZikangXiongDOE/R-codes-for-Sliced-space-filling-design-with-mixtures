# Generating Combined Constraints with Gibbs sampling
print("Gibbs_l: Gibbs sampler for generating general mixture experimental
      region with linear constraints ")
Gibbs_l=function(a0,b0,n,A=diag(c(0,0)),lb=0,x0=0)
{
  sum_a=sum(a0)
  sum_b=sum(b0)
  if((sum_b<=1)|(sum_a>=1))
  {
    print("it is necessary that sum(a)<1 and sum(b)>1")
    break
  }
  p=length(a0)
  astar=rep(0,p)
  bstar=rep(1,p)
  for(i in 1:p)
  {
    astar[i]=max(a0[i],1-sum_b+b0[i])
    bstar[i]=min(b0[i],1-sum_a+a0[i])
  }
  L_op=list()
  L_ne=list()
  for(k in 1:(p-1))
  {
    L_op=c(L_op,list(which(A[,k]>0)))   #linear constraint with positive coefficients
    L_ne=c(L_ne,list(which(A[,k]<0)))  #linear constraint with negative coefficients
  }
  udf1=function(k,x)
  {
    lower=max(astar[k],1-sum(x)-bstar[p]+x[k])
    upper=min(bstar[k],1-sum(x)-astar[p]+x[k])
    return(c(lower,upper))
  }
  
  udf2=function(k,x)
  {
    imp=L_op[[k]]
    imn=L_ne[[k]]
    iupper=1/A[imp,k]*(lb[imp]-A[imp,]%*%x+A[imp,k]*x[k]) 
    ilower=1/A[imn,k]*(lb[imn]-A[imn,]%*%x+A[imn,k]*x[k])
    lower=max(astar[k],1-sum(x)-bstar[p]+x[k],ilower)
    upper=min(bstar[k],1-sum(x)-astar[p]+x[k],iupper)
    return(c(lower,upper))
  }
  if(all(A[1,]==rep(0,ncol(A))))
  {udf=udf1}else
  {udf=udf2}
  # Sample set  
  sample=matrix(0,n,p)  
  
  ###initial point
  if(length(x0)==1)   ###Not getting  artificial initial point
  {
    x0=astar[1:(p-1)]   
    if(sum(x0)<(1-bstar[p]))
    {
      dstar=bstar[1:(p-1)]-astar[1:(p-1)]
      x0=x0+(1-bstar[p]-sum(x0)+runif(1)*(bstar[p]-astar[p]))*dstar/sum(dstar)
    }
  }else 
  {x0=x0[1:(p-1)]}
  sample[1,]=c(x0,1-sum(x0))
  x=x0
  for(m in 2:n)
  { array=sample(1:(p-1),p-1)
  for(k in array)
  {
    z=udf(k,x)     ##getting upper and low limits
    sample[m,k]=runif(1,z[1],z[2])  ##sample x_k
    x[k]=sample[m,k]      ##updating present point
  } 
  }
  sample[2:n,p]=1-apply(sample[2:n,],1,sum)
  
  ##Getting sample and initial point 
  return(list(sample=sample,xstar=sample[1,]))
}
