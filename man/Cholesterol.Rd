\name{Cholesterol}
\alias{Cholesterol}
\docType{data}
\title{
Framingham cholesterol study
}
\description{
The Framingham cholesterol study generated a benchmark dataset (Zhang and Davidian, 2001) for longitudinal analysis to examine the role of serum cholesterol as a risk factor for the evolution of cardiovascular disease for 200 randomly selected subjects.
}

\usage{data(Cholesterol)}
\format{
  This data frame contains the following columns:
  \describe{
    \item{\code{newid}}{a numeric vector indicating the subject on which the
      measurement was made. It represents the subject number in the sample.}
      \item{\code{ID}}{a numeric vector indicating the subject on which the
      measurement was made. It represents the subject number in the population.}
    \item{\code{cholst}}{cholesterol level for patient \code{newid}.}
    \item{\code{sex}}{a dichotomous gender (0=female, 1=male).}
    \item{\code{age}}{age of the patient in years.}
    \item{\code{year}}{years elapsed since the start of the study to the current measurement.}
  }
}
\source{
Zhang, D., & Davidian, M. (2001). Linear mixed models with flexible distributions of random effects for longitudinal data. Biometrics, 57(3), 795-802.
}
\references{
https://www.framinghamheartstudy.org/about-fhs/background.php
}
\examples{
\dontrun{
data(Cholesterol)
attach(Cholesterol)

y  = cholst            #response
x  = cbind(1,sex,age)  #design matrix for fixed effects
z  = cbind(1,year)     #design matrix for random effects

#A median regression
median_reg = QRLMM(y,x,newid,nj,MaxIter = 500)
}
}
\keyword{datasets}
