\name{QRLMM}
\alias{QRLMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quantile Regression for Linear Mixed-Effects Models
}
\description{
Performs a quantile regression for a LMEM using the Stochastic-Approximation of the EM Algorithm (SAEM) for an unique or a set of quantiles.
}
\usage{
QRLMM(y,x,z,nj,p=0.5,precision=0.0001,MaxIter=300,M=10,cp=0.25,
      beta=NA,sigma=NA,Psi=NA,show.convergence=TRUE,CI=95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
the response vector of dimension \eqn{N} where \eqn{N} is the total of observations.
}
  \item{x}{
design matrix for the fixed effects of dimension \eqn{N x d} where \eqn{d} represents the number of fixed effects including the intercept, if considered. 
}
  \item{z}{
design matrix for the random effects of dimension \eqn{N x q} where \eqn{q} represents the number of random effects. 
}
  \item{nj}{
Vector of dimension \eqn{n} containing the number of observations per subject. Must sum \eqn{N}.
}
  \item{p}{
unique quantile or a set of quantiles related to the quantile regression.
}
  \item{precision}{
the convergence maximum error.
}
  \item{MaxIter}{
the maximum number of iterations of the SAEM algorithm. Default = 300.
}
  \item{M}{
Number of Monte Carlo simulations used by the SAEM Algorithm. Default = 10. For more accuracy we suggest to use \code{M=20}.
}
  \item{cp}{
cut point \eqn{(0 \le cp \le 1)} which determines
the percentage of initial iterations with no memory.
}
  \item{beta}{
fixed effects vector of initial parameters, if desired. 
}
  \item{sigma}{
dispersion initial parameter for the error term, if desired.
}
  \item{Psi}{
Variance-covariance random effects matrix of initial parameters, if desired.
}
  \item{show.convergence}{
if \code{TRUE}, it will show a graphical summary for the convergence of the estimates of all parameters for each quantile in order to assess the convergence.
}
  \item{CI}{
Confidence to be used for the Confidence Interval when a grid of quantiles is provided. Default=95.
}

}
\details{
This algorithm performs the SAEM algorithm proposed by Delyon et al. (1999), a stochastic version of the usual EM Algorithm deriving exact maximum likelihood estimates of the fixed-effects and variance components.

If the initial parameters are not provided, by default, the fixed effects parameter \eqn{\beta} and dispersion parameter \eqn{\sigma} will be the maximum Likelihood Estimates for an Asymmetric Laplace Distribution (obviating the random term). See Yu & Zhang
(2005).

When a grid of quantiles is provided, a graphical summary with point estimates and Confidence Intervals for model parameters is shown and also a graphical summary for the convergence of these estimates (for each quantile), if \code{show.convergence=TRUE}.

If the convergence graphical summary shows that convergence has not be attained, it's suggested to increase \code{M} to 20, to increase the total number of iterations \code{MaxIter} to 500 or both.

About the cut point parameter \code{cp}, a number between 0 and 1 \eqn{(0 \le cp \le 1)} will assure an initial convergence in distribution to a solution neighborhood for the first \code{cp}*\code{MaxIter} iterations and an almost sure convergence for the rest of the iterations. If you do not know how SAEM algorithm works, this parameter SHOULD NOT be changed.

This program uses progress bars that will close when the algorithm ends. They must not be closed before if not the algorithm will stop.
}
\value{
The function returns a list with two objects
\item{conv}{A two elements list with the matrices \code{teta} and \code{se} containing the point estimates and standard error estimate for all parameters along all iterations.}

The second element of the list is \code{res}, a list of 12 elements detailed as

\item{iter}{number of iterations.}
\item{criteria}{attained criteria value.}
\item{beta}{fixed effects estimates.}
\item{sigma}{scale parameter estimate for the error term.}
\item{Psi}{Random effects variance-covariance estimate matrix.}
\item{SE}{Standard Error estimates.}
\item{table}{Table containing the inference for the fixed effects parameters.}
\item{loglik}{Log-likelihood value.}
\item{AIC}{Akaike information criterion.}
\item{BIC}{Bayesian information criterion.}
\item{HQ}{Hannan-Quinn information criterion.}
\item{time}{processing time.}
}
\references{
Delyon, B., Lavielle, M. & Moulines, E. (1999). Convergence of a stochastic approximation
version of the EM algorithm. Annals of Statistics, pages 94-128.

Yu, K., & Zhang, J. (2005). A three-parameter asymmetric Laplace distribution and its extension. Communications in Statistics-Theory and Methods, 34(9-10), 1867-1879.
}
\author{
Christian E. Galarza <\email{cgalarza88@gmail.com}> and 
Victor H. Lachos <\email{hlachos@ime.unicamp.br}>
}
\note{
If a grid of quantiles was provided, the result is a list of the same dimension where each element corresponds to each quantile as detailed above.
}

\seealso{
\code{\link{dALD}}, \code{\link{pALD}}, \code{\link{qALD}}, \code{\link{rALD}}
}
\examples{
#Using the Orthodontic distance growth data

library(nlme)
data(Orthodont)
attach(Orthodont)

y  = distance #response
x  = cbind(1,c(rep(0,64),rep(1,44)),age) #design matrix for fixed effects
z  = cbind(1,age) #design matrix for random effects
nj = rep(4,27)  #balanced

\dontrun{
QRLMM(y,x,z,nj,MaxIter=100) #a quick median regression


#Fit a very quick regression for the three quartiles (Just for having an idea!)
QRLMM(y,x,z,nj,p = c(0.25,0.50,0.75),MaxIter=50,M=10)

#A full profile quantile regression (This might take some time)
QRLMM(y,x,z,nj,p = seq(0.05,0.95,0.05),MaxIter=300,M=10)

#A simple output example
-------------------------------------------------
Quantile Regression for Linear Mixed Model
-------------------------------------------------
Quantile = 0.75
Subjects = 27 ; Observations = 108 ; Balanced = 4
-----------
Estimates
-----------
- Fixed effects
Estimate Std. Error z value Pr(>|z|)
beta 1 17.08405 0.53524 31.91831 0
19
beta 2 2.15393 0.36929 5.83265 0
beta 3 0.61882 0.05807 10.65643 0
sigma = 0.38439
Random effects Varcov matrix
z1 z2
z1 0.16106 -0.00887
z2 -0.00887 0.02839
------------------------
Model selection criteria
------------------------
Loglik AIC BIC HQ
Value -216.454 446.907 465.682 454.52
-------
Details
-------
Convergence reached? = FALSE
Iterations = 300 / 300
Criteria = 0.00381
MC sample = 10
Cut point = 0.25
Processing time = 7.590584 mins
}
}

\keyword{SAEM}
\keyword{quantile}
\keyword{quantile regression}
\keyword{linear mixed models}
\keyword{ALD}