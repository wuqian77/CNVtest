GLMCNV <-
function(data,L, ncase)
{
  x.mat=data

  T = dim(x.mat)[2]
  N = dim(x.mat)[1]
  
  len = (T-L)*(L+1)+(L+1)*L/2

  x = rep(0, N)

  y=rep(0, N)
  y[1:ncase]=rep(1, ncase)
  
  y_ct = y-mean(y)
  y_sd = sd(y)

###### Find CNV first  ########

start = 1
end = len


 Sabs = rep(0, len)
 Sdup = rep(0, len)
 Sdel = rep(0, len)

 caseCnv_abs=rep(0,end)
 controlCnv_abs=rep(0,end)
 
 caseCnv_dup=rep(0,end)
 controlCnv_dup=rep(0,end)
 
 caseCnv_del=rep(0,end)
 controlCnv_del=rep(0,end)

for (k in  start:end){
t = ceiling(k/(L+1));
j = min(k-(t-1)*L, T);
for (i in 1:N) {
x[i] = sum(x.mat[i,t:j])/sqrt(j-t+1);
}

zabs = abs(x) >  sqrt(2*log(T*(L+1)))
zdup = x >   sqrt(2*log(T*(L+1)))
zdel = x <  -sqrt(2*log(T*(L+1)))


caseCnv_abs[k]=sum(zabs[1:(N/2)])
controlCnv_abs[k]=sum(zabs[(1+N/2):N])

caseCnv_dup[k]=sum(zdup[1:(N/2)])
controlCnv_dup[k]=sum(zdup[(1+N/2):N])

caseCnv_del[k]=sum(zdel[1:(N/2)])
controlCnv_del[k]=sum(zdel[(1+N/2):N])


if (sum(zabs)>0) {
Sabs[k] = sum(zabs*y_ct) / (sqrt(N)*sd(zabs)*y_sd)
}

if (sum(zdup)>0) {
Sdup[k] = sum(zdup*y_ct) / (sqrt(N)*sd(zdup)*y_sd)
}

if (sum(zdel)>0) {
Sdel[k] = sum(zdel*y_ct) / (sqrt(N)*sd(zdel)*y_sd)
}


} # end k loop


###### Get candidates that pass threshold. 


#======================================================================
# Find the non-overlapped CNV regions by peaks
#======================================================================

peaks = function(index, stat, bin){


left = ceiling(index/(bin+1));
right = index-(left-1)*bin;

list = cbind(left, right, stat)

I = rep(0,2);
J = rep(0,2);
xstar = rep(0,2);

t=1;

while (length(list)> 3) {
k.max = which(abs(list[,3])== max(abs(list[,3])));

I[t] = list[min(k.max),1]
J[t] = list[min(k.max),2]
xstar[t] = list[min(k.max),3];

II = which(list[,1]<=J[t] & list[,2]>=I[t]);
list = list[-II,];
t=t+1;
}

if(length(list==3)){
I[t] = list[1]
J[t] = list[2]
xstar[t] = list[3];
}

location = cbind(I, J, xstar);
return(location);
}

###### Get candidates that pass threshold. 



#thres_fix = sqrt(2*log(T*lmax))

# Select CNVs either duplication or deletion:

z2=rbind(caseCnv_abs,controlCnv_abs)
z2_abs=ifelse((caseCnv_abs>0 | controlCnv_abs >0 ),1,0)
thresN2=sum(z2_abs)
rhat=thresN2
thres = sqrt(2*log(thresN2))
ind = which(abs(Sabs)>thres);
can = Sabs[ind];
pe = peaks(ind, can, L);
select_abs = pe[pe[,2]>pe[,1]+1,];

# Only select duplication CNVs

z3=rbind(caseCnv_dup,controlCnv_dup)
z3_dup=ifelse((caseCnv_dup>0 | controlCnv_dup >0 ),1,0)
thresN3=sum(z3_dup)
rhat=thresN3
thres = sqrt(2*log(thresN3))
ind = which(abs(Sdup)>thres);
can = Sdup[ind];
pe = peaks(ind, can, L);
select_dup = pe[pe[,2]>pe[,1]+1,];


# Only select deletion CNVs

z4=rbind(caseCnv_del,controlCnv_del)
z4_del=ifelse((caseCnv_del>0 | controlCnv_del >0 ),1,0)
thresN4=sum(z4_del)
rhat=thresN4
thres = sqrt(2*log(thresN4))
ind = which(abs(Sdel)>thres);
can = Sdel[ind];
pe = peaks(ind, can, L);
select_del = pe[pe[,2]>pe[,1]+1,];

listit = list("Select_all"=select_abs,"Select_duplication_only"=select_dup,"Select_deletion_only"=select_del )
return(listit)
}

