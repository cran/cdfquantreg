% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qrGrad.R
\name{qrGrad}
\alias{qrGrad}
\title{Give the Gradient Function for CDF-Quantile Distribution Models}
\usage{
qrGrad(fd, sd)
}
\arguments{
\item{fd}{A string that specifies the parent distribution.}

\item{sd}{A string that specifies the sub-family distribution.}
}
\value{
grad The gradient function of parameter estimates, given a specified cdf-quantile distribution
}
\description{
Give the Gradient Function for CDF-Quantile Distribution models
}
\examples{
qrGrad('t2','t2')
}
