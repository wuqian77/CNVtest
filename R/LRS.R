LRS <-
function(a, b, c, d){

#input a=a sequence of data, 
#b=length of the sequence,
#c=L,
#d=threshold.

x = rep(0, b*(c+1));   for(i in 1:(b-1)){
for (j in i:min(b, (i+c))) { 
x[(i-1)*c+j] = sum(a[i:j])/ sqrt(j-i+1)
}
}

k= which(abs(x)>d);
i = ceiling(k/(c+1)); 
j = k - (i-1)*c;
list = cbind(i,j, x[k]);

start = rep(0,1);
end = rep(0,1);
LR = rep(0,1);
t=1;

while (length(list)> 3) {
ind = which(abs(list[,3]) == max(abs(list[,3])));
start[t] = list[ind,1];
end[t] = list[ind,2];
LR[t] = list[ind,3];
II = which(list[,1]<=end[t] & list[,2]>=start[t]);
list = list[-II,];
t = t+1; 
} 

if(length(list)==3) {
start[t] = list[1];
end[t] = list[2];
LR[t] = list[3];
}

peaks = cbind(start, end, LR); 
return(peaks);
}

