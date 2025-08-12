# 产生有上下界约束的混料均匀设计
print("ITM:条件分布逆变换法产生有上下界约束的混料均匀设计")
ITM=function(u,a,b) 
{
  n=nrow(u)
  q=ncol(u)
  if((sum(b)<=1)|(sum(a)>=1))
  {
    print("it is necessary that sum(a)<1 and sum(b)>1")
    break
  }
  q=q+1             ##mapping (s-1)-dimension to s-dimension 
  y=matrix(0,n,q+1)   ##desin matrix
  delta=matrix(1,n,q)
  d=matrix(0,n,q)        ##new low limits
  phi=matrix(1,n,q)      ##new upper limits
  for(k in 0:(q-2))
  {
    i=q-k
    delta[,i]=1-apply(as.matrix(y[,(i+1):(q+1)]),1,sum)
    d[,i]=apply(as.matrix(cbind(a[i]/delta[,i],1-sum(b[1:i-1])/delta[,i])),1,max)
    phi[,i]=apply(as.matrix(cbind(b[i]/delta[,i],1-sum(a[1:i-1])/delta[,i])),1,min)
    y[,i]=delta[,i]*(1-(u[,i-1]*(1-phi[,i])^(i-1)+(1-u[,i-1])*(1-d[,i])^(i-1))^(1/(i-1)))
  }
  y[,1]=1-apply(as.matrix(y[,2:q]),1,sum)
  
  return(y[,1:q])
}
