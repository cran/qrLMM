\name{ALD}
\alias{dALD}
\alias{pALD}
\alias{qALD}
\alias{rALD}

\title{The Asymmetric Laplace Distribution
}
\description{Density, distribution function, quantile function and random generation for a Three-Parameter Asymmetric Laplace Distribution defined in Koenker and Machado (1999) useful for quantile regression with location parameter equal to \code{mu}, scale parameter \code{sigma} and skewness parameter \code{p}.
}
\usage{
dALD(y, mu = 0, sigma = 1, p = 0.5)
pALD(q, mu = 0, sigma = 1, p = 0.5, lower.tail = TRUE)
qALD(prob, mu = 0, sigma = 1, p = 0.5, lower.tail = TRUE)
rALD(n, mu = 0, sigma = 1, p = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y,q}{
vector of quantiles.
}
  \item{prob}{
vector of probabilities.
}
  \item{n}{
number of observations.
}
  \item{mu}{
location parameter.
}
  \item{sigma}{
scale parameter.
}
  \item{p}{
skewness parameter.
}
  \item{lower.tail}{
logical; if TRUE (default), probabilities are P[X \eqn{\le} x] otherwise, P[X > x].
}
}
\details{
If \code{mu}, \code{sigma} or \code{p} are not specified they assume the default values of 0, 1 and 0.5, respectively, belonging to the Symmetric Standard Laplace Distribution denoted by \eqn{ALD(0,1,0.5)}.

As discussed in Koenker and Machado (1999) and Yu and Moyeed (2001) we say that a random variable
Y is distributed as an ALD with location parameter \eqn{\mu}, scale parameter \eqn{\sigma>0} and skewness parameter \eqn{p} in (0,1), if its probability density function (pdf) is given by

\deqn{f(y|\mu,\sigma,p)=\frac{p(1-p)}{\sigma}\exp
{-\rho_{p}(\frac{y-\mu}{\sigma})}}

where \eqn{\rho_p(.)} is the so called check (or loss) function defined by
\deqn{\rho_p(u)=u(p - I_{u<0})},
with \eqn{I_{.}} denoting the usual indicator function. This  distribution is denoted by \eqn{ALD(\mu,\sigma,p)} and it's \eqn{p}-th quantile is equal to \eqn{\mu}.
}
\value{
\code{dALD} gives the density, \code{pALD} gives the distribution function, \code{qALD} gives the quantile function, and \code{rALD} generates a random sample.

The length of the result is determined by n for \code{rALD}, and is the maximum of the lengths of the numerical arguments for the other functions \code{dALD}, \code{pALD} and \code{qALD}.

The scale parameter \code{sigma} must be positive and non zero.The skew parameter \code{p} must be between zero and one (0<\code{p}<1).
}

\references{
Koenker, R., Machado, J. (1999). Goodness of fit and related inference processes for quantile
regression. J. Amer. Statist. Assoc. 94(3):1296-1309.

Yu, K. & Moyeed, R. (2001). Bayesian quantile regression. Statistics & Probability Letters, 54(4), 437-447.

Yu, K., & Zhang, J. (2005). A three-parameter asymmetric Laplace distribution and its extension. Communications in Statistics-Theory and Methods, 34(9-10), 1867-1879.
}

\author{
Christian E. Galarza <\email{cgalarza88@gmail.com}> and 
Victor H. Lachos <\email{hlachos@ime.unicamp.br}>
}
\note{
The numerical arguments other than \code{n} are recycled to the length of the result.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{QRLMM}}
}

\examples{
## Let's plot an Asymmetric Laplace Distribution!

##Density
sseq = seq(-40,80,0.5)
dens = dALD(y=sseq,mu=50,sigma=3,p=0.75)
plot(sseq,dens,type="l",lwd=2,col="red",xlab="x",ylab="f(x)", main="ALD Density function")

## Distribution Function
df = pALD(q=sseq,mu=50,sigma=3,p=0.75)
plot(sseq,df,type="l",lwd=2,col="blue",xlab="x",ylab="F(x)", main="ALD Distribution function")
abline(h=1,lty=2)

##Inverse Distribution Function
prob = seq(0,1,length.out = 1000)
idf = qALD(prob=prob,mu=50,sigma=3,p=0.75)
plot(prob,idf,type="l",lwd=2,col="gray30",xlab="x",ylab=expression(F^{-1}~(x)))
title(main="ALD Inverse Distribution function")
abline(v=c(0,1),lty=2)

#Random Sample Histogram
sample = rALD(n=10000,mu=50,sigma=3,p=0.75)
hist(sample,breaks = 70,freq = FALSE,ylim=c(0,max(dens)),main="")
title(main="Histogram and True density")
lines(sseq,dens,col="red",lwd=2)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{asymmetric laplace distribution}
\keyword{quantile regression}
\keyword{Laplace}
\keyword{ALD}
