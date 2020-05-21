Simulation_LRS <-
function(n,mu,s,loc)
{

q = length(s);
data = rnorm(n); 

for (i in 1:q){
data[loc[i]:(loc[i]+s[i]-1)] = data[loc[i]:(loc[i]+s[i]-1)] + mu 
}
return(data)
}

