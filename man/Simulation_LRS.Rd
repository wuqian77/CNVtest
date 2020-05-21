\name{Simulation_LRS}
\alias{Simulation_LRS}
\title{
Simulate CNV data for one sample 
}
\description{
Simulate CNV data for one individual. This is used for LRS method
}
\usage{
Simulation_LRS(n, mu, s, loc)
}
\arguments{
  \item{n}{
length of the DNA sequence
}
  \item{mu}{ 
average signal strength for CNV regions (jump size). 
}
  \item{s}{
length of each true CNV region
}
  \item{loc}{
location of each true CNV region
}
}


\references{
X. Jessie Jeng, T. Tony Cai, and Hongzhe Li (2010). Optimal Sparse Segment Identification with Application in Copy Number Variation Analysis. \eqn{Journal of American Statistical Association}, T&M. 105 (491), 1156-1166. 
}
\seealso{
LRS
}
\examples{
data=Simulation_LRS(n=1000,mu=3,s=c(15, 10, 5),loc=c(101, 301, 501))
}
