PickCNV_LRS <-
function(data,L)
{

obs = (data-median(data))/mad(data)    # standardization
n = length(obs)
thres = sqrt(2*log(n*L))    # set threshold 
pick = LRS(obs, n, L, thres);  # Call the LRS function in LRS.R

return(pick)
}

