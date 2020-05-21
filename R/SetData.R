SetData <-
function( loc1,loc2,n,m,mu,s,p1,p2,p3) 
{
 # initialize
 # n=segment number for one chromosome=5,000
 # m=case+control sample number =1000
 # s=signal segment= 10
 # mu is the mean value if it's in a CNV region
 
 #create y.vec : binary outcome 
 y.vec=rep(0, m)
 y.vec=replace(y.vec,1:(m/2), 1)

 #create x.mat for the whole data
 x.mat=matrix(rnorm(m*n), nrow=m, ncol=n )

 # function to generate a length s=10 interval with mu=2.5
 CNV<-function(mu,s){mu + rnorm(s)}

 # change the fix location different rate for both case 30% and control 10%
 
 for( i in 1:m) 
{
if ( i <= (p1*m/2))
{
x.mat[i,loc1:(loc1+s-1)]=CNV(mu,s)
}
else if (i <= (m/2+p2*m/2) & i > (m/2))
{
x.mat[i,loc1:(loc1+s-1)]=CNV(mu,s)
}
}

 x1.mat=x.mat

 # change the fix location same rate for both case 10% and control 10%

 for( i in 1:m) 
{
if ( i <= (p3*m/2))
{
x.mat[i,loc2:(loc2+s-1)]=CNV(mu,s)
}
else if (i <= (m/2+p3*m/2) & i > (m/2))
{
x.mat[i,loc2:(loc2+s-1)]=CNV(mu,s)
}
}

 x2.mat=x.mat
 return(x2.mat)

## If needed, CNVtest can also detect the CNVs with location shift between carriers
## Set a new loc3 with length s+5 and each CNV carrier has CNV with length s
## there are totally 6 different choices
## CNV regions are overlapped with each other but not exactly the same loci.  

 #loc3=round(runif(1,1,(n-s)))
 # loc3=301
 # loc3interval=c(loc3:(loc3+s-1))
 # chooseloc3=rbind(loc3interval,loc3interval+1,loc3interval+2,loc3interval+3,loc3interval+4,loc3interval+5,loc3interval+6)

 #for( i in 1:m) 
 #{
 #if ( i <= (ploc3*m/2))
 #{
#temp1=ceiling(7*runif(1))
#x.mat[i,chooseloc3[temp1,]]=CNV(mu,s)
#}
#else if (i <= (m/2+qloc3*m/2) & i > (m/2))
#{
#temp2=ceiling(7*runif(1))
#x.mat[i,chooseloc3[temp2,]]=CNV(mu,s)
#}
#}
# x3.mat=x.mat
}

