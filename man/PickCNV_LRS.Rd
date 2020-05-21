\name{PickCNV_LRS}
\alias{PickCNV_LRS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pick the CNV regions
}
\description{
Applicaiton of LRS for CNV identification using SNP data
}
\usage{
PickCNV_LRS(data, L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
DNA sequence of real data for one sample
}
  \item{L}{
it is chosen to be longer than the majority of the segments, usually L < 20. }
}

\references{
X. Jessie Jeng, T. Tony Cai, and Hongzhe Li (2010). Optimal Sparse Segment Identification with Application in Copy Number Variation Analysis. \eqn{Journal of American Statistical Association}, T&M. 105 (491), 1156-1166. 
}
\seealso{
LRS, Simulation_LRS
}
\note{
Note 1: L can be selected as 20 since the CNV regions mostly range below 20 SNPs. 

Note 2: If you want to select CNV regions with more than 2 SNPs, follow the example below. 
}

\examples{

# input
data=Simulation_LRS(n=1000,mu=3,s=c(15, 10, 5),loc=c(101, 301, 501))
L=20    # L is chosen to be longer than the majority of the segments. 
Example_pick=PickCNV_LRS(data,L)

#You will get a list of selected intervals with their likelihood ratio statistics:
# Example_pick
#     start end        LR
#[1,]   101 115 10.510696
#[2,]   301 310  9.845136
#[3,]   501 505  5.688137

# If you want to select CNV regions with more than 2 SNPs as follows:
k=3
picklist = Example_pick[Example_pick[,1]<(Example_pick[,2]-k+2),]
}
