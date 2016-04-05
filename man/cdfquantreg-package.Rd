\name{cdfquantreg-package}
\alias{cdfquantreg-package}
\docType{package}
\title{
Quantile Regression for Random Variables on the Unit Interval
}
\description{Employs a two-parameter family of distributions for modelling random variables on the (0, 1) interval by applying the cumulative distribution function (cdf) of one parent distribution to the quantile function of another. 
}
\details{
\tabular{ll}{
Package: \tab cdfquantreg\cr
Type: \tab Package\cr
Version: \tab 1.0.4\cr
Date: \tab 2015-04-05\cr
License: \tab GPL-3\cr
}
The cdfquantreg package includes 15 members of a two-parameter family of distributions for modelling random variables on the (0, 1) interval (see \link{cdfqrFamily}).  This family has explicit pdfs, cdfs, and quantile functions.  The two parameters include a location parameter and a dispersion parameter. The location parameter models the median and simple functions of the dispersion parameter model other quantiles (see Smithson and Shou, 2015, for details about the distribution family and the models). Separate submodels may be specified for the location and for the dispersion parameters, with different or overlapping sets of predictors in each.\cr\cr
The package offers maximum likelihood (see \link{cdfquantreg}), Bayesian MCMC (see \link{bugsModel}), and bootstrap (see \link{qrBoot}) estimation methods.  All model functions return S3 objects. The bugsModel() function runs OpenBUGS in R. In addition to the usual goodness of fit information, the package provides root-mean-squared errors in both the raw and logit scales, and the gradient. Model diagnostics include raw, Pearson, and deviance residuals (see \link{residuals.cdfqr}), and dfbetas (see \link{influence.cdfqr}). \cr\cr
For each of the 15 distributions, the package provides evaluations of the pdf (\link{dq}), cdf (\link{pq}), and quantile (\link{qq}), as well as random samples from any of them (\link{rq}). Evaluations of skew and kurtosis (\link{qrPwlm}) also are available using probability-weighted L-moments. 
}
\author{
Yiyun Shou (\email{yiyun.shou@@anu.edu.au}) and Michael Smithson (\email{Michael.Smithson@@anu.edu.au})

Maintainer: Yiyun Shou
}
\references{
Smithson, M. and Shou, Y. (2015). CDF-Quantile Distributions for Modeling Random Variables on the Unit Interval. Unpublished Manuscript, The Australian National University, Canberra, Australia.
}

\keyword{ package }
\seealso{
\code{\link{cdfqrFamily}}
}

