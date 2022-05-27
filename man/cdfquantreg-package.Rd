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
Version: \tab 1.3.2\cr
Date: \tab 2022-05-19\cr
License: \tab GPL-3\cr
}

The cdfquantreg package includes 36 members of a two-parameter family of distributions for modelling random variables on the (0, 1) interval (see \link{cdfqrFamily}).  
This family has explicit pdfs, cdfs, and quantile functions.  The two parameters consist of a location parameter and a dispersion parameter. The location parameter models the median and the dispersion parameter models the spread of other quantiles around the median (see Smithson and Shou, 2016, for details about the distribution family and the models). Separate submodels may be specified for the location and for the dispersion parameters, permitting different or overlapping sets of predictors in each.

The package offers maximum likelihood (see \link{cdfquantreg})and bootstrap (see \link{qrBoot}) estimation methods.  All model functions return S3 objects. In addition to the usual goodness of fit information, the package provides root-mean-squared errors in both the raw and logit scales, and the gradient. Model diagnostics include raw, Pearson, and deviance residuals (see \link{residuals.cdfqr}), and dfbetas (see \link{influence.cdfqr}). 

For each distribution, the package provides evaluations of the pdf (\link{dq}), cdf (\link{pq}), and quantile (\link{qq}), as well as random samples from any of them (\link{rq}). Evaluations of skew and kurtosis (\link{qrPwlm}) also are available using probability-weighted L-moments. 
}

\author{
Yiyun Shou and Michael Smithson \cr

Maintainer: Yiyun Shou (\email{yiyun.shou@anu.edu.au})
}
\references{
Shou, Y. and Smithson, M., (2019). {cdfquantreg}: An {R} Package for {CDF}-Quantile Regression. Journal of Statistical Software,88(1), pp.1--30, doi: 10.18637/jss.v088.i01
}

\keyword{ package }

\seealso{
\code{\link{cdfqrFamily}}
}

