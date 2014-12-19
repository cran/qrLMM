\name{qrLMM-package}
\alias{qrLMM-package}
\docType{package}
\title{
Package for Quantile Regression for Linear Mixed-Effects Model
}
\description{
This package contains a principal function that performs a quantile regression for a Linear Mixed-Effects Model using the Stochastic-Approximation of the EM Algorithm (SAEM) for an unique or a set of quantiles. It also contains the d,p,q and r functions for the Asymmetric Laplace Distribution defined in Koenker and Machado (1999) useful for quantile regression. 

Exploiting the nice hierarchical representation of the ALD, our classical approach follows the Stochastic Approximation of the EM(SAEM) algorithm for deriving exact maximum likelihood estimates of the fixed-effects and variance components.
}
\details{
\tabular{ll}{
Package: \tab qrLMM\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-12-19\cr
License: \tab What license is it under?\cr
}
}
\author{Christian E. Galarza <\email{cgalarza88@gmail.com}> and 
Victor H. Lachos <\email{hlachos@ime.unicamp.br}>

Maintainer: Christian E. Galarza <\email{cgalarza88@gmail.com}>
}
\references{
Delyon, B., Lavielle, M. & Moulines, E. (1999). Convergence of a stochastic approximation
version of the EM algorithm. Annals of Statistics, pages 94-128.

Koenker, R., Machado, J. (1999). Goodness of fit and related inference processes for quantile
regression. J. Amer. Statist. Assoc. 94(3):1296-1309.

Yu, K. & Moyeed, R. (2001). Bayesian quantile regression. Statistics & Probability Letters, 54(4), 437-447.

Yu, K., & Zhang, J. (2005). A three-parameter asymmetric Laplace distribution and its extension. Communications in Statistics-Theory and Methods, 34(9-10), 1867-1879.
}

\keyword{package}
\keyword{quantile}
\keyword{ALD}
\seealso{\code{\link{dALD}}, \code{\link{pALD}}, \code{\link{qALD}}, \code{\link{rALD}}, \code{\link{QRLMM}}
}
\examples{
#See examples for the QRLMM function linked above.
}