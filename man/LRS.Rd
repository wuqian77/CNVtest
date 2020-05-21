\name{LRS}
\alias{LRS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Likelihood Ratio Selection (LRS)
}
\description{
Likelihood Ratio selection (LRS) function is to identify the sparse segment in whole genome, e.g. CNV regions. 
}
\usage{
LRS(a, b, c, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
input datasets a, which is the sequence of the CNV data for one sample 
}
  \item{b}{
n=b, length of the sequence
}
  \item{c}{
L=c, which is chosen to be longer than the majority of the segments. 
}
  \item{d}{
threshold used in the lrs algorithm: \eqn{d=\sqrt{2log(bc)}}
}
}

\references{
X. Jessie Jeng, T. Tony Cai, and Hongzhe Li (2010). Optimal Sparse Segment Identification with Application in Copy Number Variation Analysis. \eqn{Journal of American Statistical Association}, T&M. 105 (491), 1156-1166. 
}
\seealso{
Simulation_LRS
}

\examples{

data=Simulation_LRS(n=1000,mu=3,s=c(15, 10, 5),loc=c(101, 301, 501))
obs = (data-median(data))/mad(data)    # standardization
n = length(obs)
L=20
thres = sqrt(2*log(n*L))    # set threshold 

pick = LRS(obs, n, L, thres);  # Call the lrs function in LRS.R

k=3
picklist = pick[pick[,1]<(pick[,2]-k+2),]


}

