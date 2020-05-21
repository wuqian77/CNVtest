\name{GLMCNV}
\alias{GLMCNV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify the trait-associated CNVs via generalized linear model
}
\description{
GLMCNV is developed to identify the copy number variants (CNVs) that are associated with disease. Methods for testing the CNV association are limited. Most available methods use a two-step approach, where the CNVs carried by the samples are identified first and then tested for association. But, GLMCNV directly identify the trait-associated CNVs without the need of identifying sample-specific CNVs. 
}
\usage{
GLMCNV(data, L, ncase)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data is the input data matrix with dimension \eqn{N*T}. 

N is the totally sample size. T is the length of the genome. The first ncase row is from case group. 
}
  \item{L}{
L is the maximum length of the sliding window, which is used in the algorithm to detect CNV region. 
}
  \item{ncase}{
ncase is the sample size of case. N - ncase is the sample size of control
}
}

\details{
L need to satisfy some conditions. It won't cover two CNV region at the same time and also can cover every CNV signals. In addition, \eqn{ log L = o(log T)}. Usually L is less than 20.
}

\value{
Output list: 
\item{Select_all}{ Select both CNV duplication and deletion   }
\item{Select_duplication_only}{ Only select CNV duplication}
\item{Select_deletion_only}{ Only select CNV deletion   }
}
\references{
X. Jessie Jeng, Qian Wu and Hongzhe Li (2015) A Statistical Method for Identifying Trait-Associated Copy Number Variants. \eqn{Human heredity};79:147-156
}

\seealso{
\code{\link{SetData}}
}

\examples{
data=SetData(loc1=304,loc2=3058,n=5000,m=1000,mu=2,s=10,p1=0.2,p2=0.1,p3=0.15)
L=15;
ncase=500;
}


