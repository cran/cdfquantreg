#' @title Censored CDF-Quantile Probability Distributions
#' @aliases cdfquantregC
#' @description \code{cdfquantregC} is the a function to fit a censored cdf quantile regression with a variety of distributions . 
#' 
#' @param formula A formula object, with the dependent variable (DV) on the left of an ~ operator, and predictors on the right. For the part on the right of '~', the specification of the location and dispersion submodels can be separated by '|'. So \code{y ~ X1 | X2} specifies that the DV is \code{y}, \code{X1} is the predictor in the location submodel, and \code{X2} is the predictor in the dispersion submodel.
#' @param data The data in a data.frame format 
#' @param fd A string that specifies the parent distribution.
#' @param sd A string that specifies the child distribution.
#' @param c1 The left censored value, if NULL, the minimum value in the data will be used
#' @param c2 The right censored value, if NULL, the maximum value in the data will be used
#' @param family If `fd` and `sd` are not provided, the name of a member of the family of distributions can be provided (See \code{\link{cdfqrFamily}} for details of family functions)
#' @param censor A string variable to indicate how many censored point is used- only left censored \code{`LC`}, or only right-hand censored \code{`RC`}, or both sides \code{`DB`}.
#' @param start The starting values for model fitting. If not provided, default values will be used.
#' @param control Control optimization parameters  (See \code{\link{cdfqr.control}}))
#' @param ... Currently ignored. 
#'
#' @details The cdfquantreg function fits a quantile regression model with a distributions from the cdf-quantile family selected by the user (Smithson and Shou, 2015). The model is specified in a two-part formula, one part containing the predictors of the location parameter, and the second part containing the predictors of the dispersion parameter.  The models are fitted in two stages, the first of which uses the Nelder-Mead algorithm and the second of which takes the estimates from the first stage and applies the BFGS algorithm to refine the estimates.  
#' 
#' @export
#' @import Formula
#' 
#' @return An object of class \code{cdfquantreg} will be returned. Generic functions such as \link{summary},\link{print} (e.g.,  \link{print.cdfqr}) and \link{coef} can be used to extract output (see \link{summary.cdfqr} for more details about the generic functions that can be used).
#' Class of object is a list with the following output:
#'  \describe{
#'   \item{coefficients}{A named vector of coefficients.}
#'   \item{residuals}{Raw residuals, the difference between the fitted values and the data.}
#'   \item{fitted}{The fitted values, including full model fitted values, fitted values for the mean component, and fitted values for the dispersion component.}
#'   \item{rmse}{The model root mean squared errors}
#'   \item{rmseLogit}{The root mean squared errors between the logit of the fitted values, and the logit of the response values.}
#'   \item{vcov}{The variance-covariance matrix of the coefficient estimates.}
#'   \item{AIC, BIC}{Akaike's Information Criterion and Bayesian Information Criterion.}
#'   \item{deviance}{The deviance for the model.}
#' }
#' 
#' @examples
#' data(cdfqrExampleData)
#' fit <- cdfquantregC(crc99 ~ vert | confl, c1 = 0.001, c2= 0.999, 
#' fd ='t2',sd ='t2', data = JurorData)
#' 
#' summary(fit)

cdfquantregC <- function(formula, fd = NULL, sd = NULL, 
                         data, family = NULL, censor = 'DB',
                         c1 = NULL, c2 = NULL, start = NULL, control = cdfqr.control(...), 
  ...) {
  
  # ***************************** Diagnose the input, extract input values
  input_full <- match.call()
  
  input <- match.call(expand.dots = FALSE)
  input_index <- match(c("formula", "data", "fd", "sd", "family", "start",
                         "c1", "c2", "t1", "t2"), names(input), 
    0L)
  input <- input[c(1L, input_index)]
  
  if (is.null(family) & (is.null(fd) | is.null(sd))) {
    stop("No distribution is specified!")
  }
  
  if (!is.null(family)) {
    fd <- strsplit(family,"-")[[1]][1]
    sd <- strsplit(family,"-")[[1]][2]
  }
  

  
  # check the data if data frame is not provided, try to find the data in the
  # global environment
  if (missing(data)) 
    data <- environment(formula)
  
  for (i in 1:ncol(data)) {
    if (!is.numeric(data[, i])) 
      data[, i] <- factor(data[, i])
  }  #turn character variables to factor 
  
  # check the formula
  formula <- Formula::as.Formula(formula)
  if (length(formula)[2] < 2) {
    # if the user only specifies the mean submodel, add the dispersion submodel intercept in to
    # the formula
    formula <- deparse(formula)
    formula <- paste(formula, "| 1", sep = "")
    formula <- Formula::as.Formula(formula)
  }
  
  # Process unstandard input
  fd <- tolower(fd)
  sd <- tolower(sd)
  
  
  # ***************************** Process data and formula
  data <- model.frame(formula, data)  # get the data with variables in the model only
  ydata <- as.matrix(model.frame(formula, data, rhs = 0))
  
  if (any(ydata <= 0) | any(ydata >= 1)) {
    warning(paste("The values of the dependent variable is rescaled into (0, 1) interval", "\n", "see scaleTR() for details"))
    ydata <- scaleTR(ydata)
  }
   
  
  data_lm <- as.matrix(model.matrix(formula, data, rhs = 1))  # datamatrix for location model
  data_pm <- as.matrix(model.matrix(formula, data, rhs = 2))  # datamatrix for dispersion model
  
  n_lm <- ncol(data_lm)  #number of parameters in mean component
  n_pm <- ncol(data_pm)  #number of parameters in dispersion component
  
  n <- nrow(data)  #number of responses
  k <- n_lm + n_pm  # total number of parameters
  
  # ***************************** Obtain other parameter specifications, get the
  # distribution specific loglikelihood function and grad
  grad <- qrGrad(fd, sd)
  
  # starting values
  if(is.null(start)){
    # if the user does not supply starting values, generate some default values
    start <- rep(0.1, k)
    start_0 <- qrStart(ydata,fd=fd,sd=sd)
    start[1] <- start_0[1] #starting value for mu intercept
    start[n_lm + 1] <- start_0[2] #starting value for sigma intercept
  }
  
  # Get the censored points
  if(is.null(c1)) c1 <- min(ydata)
  if(is.null(c2)) c2 <- max(ydata)
  q1 <- ifelse(ydata == c1, 1, 0)
  q2 <- ifelse(ydata == c2, 1, 0)
  
 
  if(censor == "LC") {c2 = NULL;  q2 = NULL}
  if(censor == "RC") {c1 = NULL;  q1 = NULL}
  
  if(censor == 'DB'){
    cdfqint <- function(h, y, c1, c2, q1, q2, x, z, fd, sd)
    {
      hx = x%*%h[1: length(x[1,])]
      mu = hx
      gz = z%*%h[length(x[1,])+1: length(z[1,])]
      sigma = exp(gz)
      loglik = log(q1*pq(c1, mu, sigma, fd, sd) + 
                     q2*(1 - pq(c2, mu, sigma, fd, sd)) + 
                     pmin(1-q1, 1-q2)*dq(y, mu, sigma, fd, sd))
      -sum(loglik, na.rm = TRUE)
    }
    preopt <- optim(start, fn = cdfqint, hessian = F, x = data_lm, y = ydata, 
                    z = data_pm, c1 = c1, c2 = c2, q1 = q1, q2 = q2, fd = fd, sd = sd,
                    method = "Nelder-Mead")
    
    contr.method <- control$method
    #contr.hessian <- control$hessian
    contr.hessian <- TRUE
    contr.trace <- control$trace
    contr.maxit <- control$maxit
    
    betaopt <- optim(preopt$par, fn = cdfqint, x = data_lm, y = ydata, 
                     z = data_pm, c1 = c1, c2 = c2, q1 = q1, q2 = q2, fd = fd, sd = sd,
                     hessian = contr.hessian, method = contr.method, 
                     control=list(trace = contr.trace , maxit = contr.maxit))
    
    
  }else{
    tt <- c(c1, c2)
    qtc <- c(q1, q2)
    cdfqint1 <- function(h, y, tt, qtc, x, z, fd, sd)
    {
      hx = x%*%h[1: length(x[1,])]
      mu = hx
      gz = z%*%h[length(x[1,])+1: length(z[1,])]
      sigma = exp(gz)
      loglik = log(qtc*pq(tt, mu, sigma, fd, sd) + 
                     pmin(1-qtc)*dq(y, mu, sigma, fd, sd))
      -sum(loglik, na.rm = TRUE)
    }
    preopt <- optim(start, fn = cdfqint1, hessian = F, x = data_lm, y = ydata, 
                    z = data_pm, tt = tt, qtc = qtc, fd = fd, sd = sd,
                    method = "Nelder-Mead")
    
    contr.method <- control$method
    #contr.hessian <- control$hessian
    contr.hessian <- TRUE
    contr.trace <- control$trace
    contr.maxit <- control$maxit
    
    betaopt <- optim(preopt$par, fn = cdfqint1, x = data_lm, y = ydata, 
                     z = data_pm, tt = tt, qtc = qtc, fd = fd, sd = sd,
                     hessian = contr.hessian, method = contr.method, 
                     control=list(trace = contr.trace , maxit = contr.maxit))
    
    # ***************************** Process the output parameter
  }

  

 

  # ***************************** Fit the model
  # Gradient
  gradients <- grad(betaopt$par,ydata,data_lm,data_pm)
  
  # Convergence message from optim
  converge <- betaopt$convergence
  
  # estimation--------------------
  estim <- betaopt$par  #mean estiamte
  serr <- sqrt(sqrt(diag(solve(betaopt$hessian))^2))  # se
  zstat <- estim/serr  # z-test
  prob <- 2 * (1 - pnorm(abs(zstat)))  # p value
  
  est <- cbind(estim <- estim, 
               serr <- serr, 
               zstat <- estim/serr,
               prob <- 2 * (1 - pnorm(abs(zstat))))
  colnames(est) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  
  # Separate output depending on the submodels
  est_lm <- est[1:n_lm, ]  #location model
  if(n_lm==1) est_lm <- matrix(est_lm, ncol=4)
  
  rownames(est_lm) <- colnames(data_lm)
  colnames(est_lm) <- colnames(est)
  
  est_pm <- est[(1 + n_lm):nrow(est), ]  #dispersion model
  if(n_pm==1) est_pm <- matrix(est_pm, ncol=4)
  rownames(est_pm) <- colnames(data_pm)
  colnames(est_pm) <- colnames(est)
  
  coefficients <- list(location = est_lm, dispersion = est_pm)
  
  
  
  # Generall fit-------------------- Loglike
  logLiklihod <- -betaopt$value
  attr(logLiklihod, "nall") <- n
  attr(logLiklihod, "nobs") <- n
  attr(logLiklihod, "df") <- k
  class(logLiklihod) <- "logLik"
  
  # AIC & BIC
  AIC <- AIC(logLiklihod)
  BIC <- BIC(logLiklihod)
  
  # Fitted values (will be quantile estimates corresponding to seq(0,N)/N, where N
  # is the sample size)
  if(n_lm == 1) {fitted_mu <- data_lm * est_lm[, 1]
  }else{
    fitted_mu <- rowSums(t(apply(data_lm, 1, 
                                 function(x) x * est_lm[, 1])))
  }
  if  (fd =="km"){
    #log link for parameter a in Km distribution
    fitted_mu <- exp(fitted_mu)
  }
  

  #dispersion model
  if(n_pm == 1) {fitted_phi <- data_pm * est_pm[, 1]
  }else{
    fitted_phi <- rowSums(t(apply(data_pm, 1, function(x) x * est_pm[, 1])))
  }
 
  fitted_sigma <- exp(fitted_phi)
  

  q_y <- seq(1e-06, 0.999999, length.out = n)
  cdf_y <- q_y[rank(ydata)]
  
  fitted <- qq(cdf_y, fitted_mu, fitted_sigma, fd, sd)
  
 
  # Raw residuals
  residuals <- ydata - fitted
  df.residual <- n - k
  
  # Deviance
  deviance <- 2 * (qrLogLik(ydata, fitted_mu, fitted_sigma, fd, sd) - 
                      qrLogLik(fitted,  fitted_mu, fitted_sigma, fd, sd))
  
  # Parameter variance-covariance matrix
  vcov <- solve(as.matrix(betaopt$hessian))
  rownames(vcov) <- colnames(vcov) <- c(paste('(mu)_', rownames(est_lm), sep = ""),
                                        paste('(sigma)_', rownames(est_pm), sep = ""))
  
  
  # RMSE
  rmse <- sqrt(mean(residuals^2))
  
  rmseLogit <- NULL
  # RMSlogit(E)
  for (i in 1:n){
   if(is.nan(fitted[i])|is.na(fitted[i])){
           rmseLogittemp <- NA
   }else{
         rmseLogittemp <-sqrt(mean((logit(fitted[i]) - logit(ydata[i]))^2))
      }
  
      rmseLogit <- c(rmseLogit, rmseLogittemp)
  }
  
  fitted = switch(fd,
                  list(full = fitted, mu = fitted_mu, sigma = fitted_sigma),
                  km = list(full = fitted, a = fitted_mu, b = fitted_sigma))
  
  # Output
  output <- list(family = list(fd = fd, sd = sd), 
                 coefficients = coefficients, 
                 call = input,  formula = formula, 
                 N = n, k = k, y = ydata, 
                 residuals = residuals, 
                 fitted = fitted, 
                 vcov = vcov, rmse = rmse, rmseLogit = rmseLogit, 
                 logLik = logLiklihod, deviance = deviance, 
                 aic = AIC, bic = BIC, 
                 converged = converge, grad = 0,
                 df.residual = df.residual, optim = betaopt)
  
  class(output) <- "cdfqr"
  
  rm(grad)
  #print.cdfqr(output)
  return(output)
  
} 
