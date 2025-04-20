\name{energy_method_complex}
\alias{energy_method_complex}
\title{Implements the two sample paired or independent energy method}
\description{
  This function takes two samples of complex-valued high-dimensional functional data, implements the energy
  method, and returns a p-value for the global test of equality of mean and a channel-wise p-value
  for each functional coordinate. 
}
\usage{
energy_method_complex(sample_1, sample_2, num_bootstrap_reps, seed, type)
}
\arguments{
  \item{sample_1}{A three dimensional complex array with dimension attribute (K,T,n_1) where K is the number of channels, T is the number of functional recordings, and n_1 is the sample size. }
  \item{sample_2}{A three dimensional complex array with dimension attribute (K,T,n_1) where K is the number of channels, T is the number of functional recordings, and n_2 is the sample size.}
    \item{num_bootstrap_reps}{A number. The number of bootstrap resamples to use when implementing the test}
     \item{seed}{A number. The seed used for randomness in bootstrap procedure }
       \item{type }{A sting. Either "paired" or "independent""}
}
\value{
 A list containg  the p-values of the test for the global hypothesis and channel-wise hypotheses, as
 well as summary information about the samples. 
}
\examples{
K=10
T=100
n_1=10
n_2=20
  sample_1_real=array(rnorm(K*T*n_1), dim=c(K,T,n_1))
  sample_1_complex=array(rnorm(K*T*n_1), dim=c(K,T,n_1))
  
  sample_1<-array(complex(real=sample_1_real, imaginary=sample_1_complex), dim=c(K,T,n_1))
  
  sample_2_real=array(rnorm(K*T*n_2), dim=c(K,T,n_2))
  sample_2_complex=array(rnorm(K*T*n_2), dim=c(K,T,n_2))
  
  sample_2<-array(complex(real=sample_2_real, imaginary=sample_2_complex), dim=c(K,T,n_2))
energy_method_complex(sample_1, sample_2, num_bootstrap_reps=1000, seed=123, type="independent")

}

\author{
  David Colin Decker
}
\references{
  Article on energy method forthcoming
}