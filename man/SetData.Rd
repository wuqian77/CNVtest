\name{SetData}
\alias{SetData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate a case-control CNV datasets
}
\description{
Generate a Genome-wide CNV association study via Monte Carlo simulation. Simulate two CNVs: Null CNV at loc2 and Disease-associated CNV at loc1. Null CNV with the same \eqn{p3\%} proportion of CNV carriers in both case and control group. Disease-associated CNV has a frequency of \eqn{p1 \%} in case and \eqn{p2\%} in control. 
}
\usage{
SetData(loc1, loc2, n, m, mu, s, p1, p2, p3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{loc1}{
CNV location for true trait-associated CNV
}
  \item{loc2}{
CNV location for null CNV
}
  \item{n}{
Length of the DNA sequence for each individual
}
  \item{m}{
total sample size m=ncase+ncontrol
}
  \item{mu}{
The shift mean of the CNV
}
  \item{s}{
the length of simulated CNV regions
}
  \item{p1}{
the proportion of CNV carriers for case at loc1
}
  \item{p2}{
the proportion of CNV carriers for control at loc1
}
  \item{p3}{
the same proportion of CNV carriers at loc2 (Null CNV)
}
}
\details{
Each CNV carriers has CNV at the same location. This could be modified by setting a new loc3, which allows CNV carriers has overlapped CNV but does not require fix location. 
}

\references{
X. Jessie Jeng, Qian Wu and Hongzhe Li (2015) A Statistical Method for Identifying Trait-Associated Copy Number Variants. \eqn{Human heredity};79:147-156
}

\seealso{
CNVtest
}

