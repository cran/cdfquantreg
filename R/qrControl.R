#' @title Control Optimization Parameters for CDF-Quantile Probability Distributions
#' @description Control Optimization Parameters for CDF-Quantile Probability Distributions.
#' @aliases cdfqr.control
#' @export
#' @param method Characters string specifying the method argument passed to \link[stats]{optim}.
#' @param maxit Integer specifying the maxit argument (maximal number of iterations) passed to \link[stats]{optim}.
#' @param trace Logical or integer controlling whether tracing information on the progress of the optimization should be produced
#' @return A list with the arguments specified.
#' @examples
#' 
#' data(cdfqrExampleData)
#' fit <- cdfquantreg(crc99 ~ vert | confl, 't2', 't2', 
#' data = JurorData,control = cdfqr.control(trace = TRUE))


cdfqr.control<- function(method = "BFGS", maxit = 5000, trace = FALSE){
  
  val <- list(method = method, maxit = maxit,  trace = trace)
  
  val
}

