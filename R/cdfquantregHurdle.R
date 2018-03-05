#' @title Zero/One inflated CDF-Quantile Probability Distributions
#' @aliases cdfquantregH
#' @description \code{cdfquantregH} is the a function to fit a Zero/One inflated CDF-Quantile regression with a variety of distributions . 
#' 
#' @param formula A formula object, with the dependent variable (DV) on the left of an ~ operator, and predictors on the right. For the part on the right of '~', the specification of the location and dispersion submodels can be separated by '|'. So \code{y ~ X1 | X2} specifies that the DV is \code{y}, \code{X1} is the predictor in the location submodel, and \code{X2} is the predictor in the dispersion submodel.
#' @param zero.fo A formula object to indicate the predictors for the zero component, only input as \code{~ predictors} 
#' @param one.fo A formula object to indicate the predictors for the one component, only input as \code{~ predictors} 
#' @param type A string variable to indicate whether the model is zero-inflated \code{`ZI`}, or one-inflated \code{`OI`}, or zero-one inflated \code{`ZO`}.
#' @param data The data in a data.frame format 
#' @param fd A string that specifies the parent distribution.
#' @param sd A string that specifies the child distribution.
#' @param family If `fd` and `sd` are not provided, the name of a member of the family of distributions can be provided (See \code{\link{cdfqrFamily}} for details of family functions)
#' @param start The starting values for model fitting. If not provided, default values will be used.
#' @param control Control optimization parameters  (See \code{\link{cdfqr.control}}))
#' @param ... Currently ignored. 
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
#'   \item{vcov}{The variance-covariance matrix of the coefficient estimates.}
#'   \item{AIC, BIC}{Akaike's Information Criterion and Bayesian Information Criterion.}
#' }
#' 
#' @examples
#' data(cdfqrExampleData)
#' # For one-inflated model
#' ipcc_high <- subset(IPCC, mid == 1 & high == 1 & prob!=0)
#' fit <- cdfquantregH(prob ~ valence | valence,one.fo = ~valence,
#'   fd ='t2',sd ='t2', type = "OI", data = ipcc_high)
#' 
#' summary(fit)
#' 
#' # For zero-inflated model
#' ipcc_low <- subset(IPCC, mid == 0 & high == 0 & prob!=1)
#' fit <- cdfquantregH(prob ~ valence | valence, zero.fo = ~valence,
#'   fd ='t2',sd ='t2', type = "ZI", data = ipcc_low)
#'  
#'  
#' # For zero &one-inflated model
#' ipcc_mid <- subset(IPCC, mid == 1 & high == 0)
#' fit <- cdfquantregH(prob ~ valence | valence, zero.fo = ~valence,
#'   one.fo = ~valence,
#'   fd ='t2',sd ='t2', type = "ZO", data = ipcc_mid)
#'  
#'   

cdfquantregH <- function(formula, zero.fo = ~1, one.fo = ~1, fd = NULL, sd = NULL, 
                         data, family = NULL, type = 'ZI',
                         start = NULL, control = cdfqr.control(...), 
                         ...) {
  
  # Diagnose the input, extract input values and check data----
  input_full <- match.call()
  
  input <- match.call(expand.dots = FALSE)
  input_index <- match(c("formula", "data", "fd", "sd", "family", "start",
                         "zero.fo", "one.fo", "type"), names(input), 
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
  
  
  # Process data for cdf model part and formula ----
  datacdf <- model.frame(formula, data)  # get the data with variables in the model only
  ydata <- as.matrix(model.frame(formula, datacdf, rhs = 0))
  
  if (any(ydata < 0) | any(ydata > 1)) {
    warning(paste("The values of the dependent variable is rescaled into (0, 1) interval", "\n", "see scaleTR() for details"))
    ydata <- scaleTR(ydata)
  }
  
  # Process submodels
  data_lm <- as.matrix(model.matrix(formula, datacdf, rhs = 1))  # datamatrix for location model
  data_pm <- as.matrix(model.matrix(formula, datacdf, rhs = 2))  # datamatrix for dispersion model
  
  n_lm <- ncol(data_lm)  #number of parameters in mean component
  n_pm <- ncol(data_pm)  #number of parameters in dispersion component
  
  n <- nrow(datacdf)  #number of responses
  k <- n_lm + n_pm  # total number of parameters
  
  # Process data for zero/one model part and formula----
  qy1 <- ifelse(ydata == 0, 1, 0)
  qy2 <- ifelse(ydata == 1, 1, 0)
  
  zero.fo <- Formula::as.Formula(zero.fo) 
  one.fo <- Formula::as.Formula(one.fo) 
  
  data_zero <- as.matrix(model.matrix(zero.fo, data, rhs = 1))
  n_zm <- ncol(data_zero)  #number of parameters in zero component
  
  data_one <- as.matrix(model.matrix(one.fo, data, rhs = 1))
  n_om <- ncol(data_one)  #number of parameters in one component
  
  if (type =="ZI") k <- n_lm + n_pm + n_zm # total number of parameters
  if (type =="OI") k <- n_lm + n_pm + n_om # total number of parameters
  if (type =="ZO") k <- n_lm + n_pm + n_zm + n_om # total number of parameters

  
  # Get log likelihood----
  quantreg_loglik <- qrLogLikFun(fd, sd)
  # Gradient
  grad <- qrGrad(fd, sd)
  
  # zeroonepart <- function(h, y, qy, x, z, w, fd, sd)
  # {
  #   w<- as.matrix(w)
  #   x<- as.matrix(x)
  #   z<- as.matrix(w)
  #   
  #   hx = w%*%h[1: ncol(w)]
  #   h1 <- h[-c(1:  ncol(w))]
  #   hx <- x %*% h1[1:length(x[1, ])]
  #   mu <- hx
  #   gz <- z %*% h1[length(x[1, ]) + 1:length(z[1, ])]
  #   sigma <- exp(gz)
  #   
  #   loglik <- rep(0, length(y))
  #   for (i in 1:length(y)){
  #     if(qy[i] == 0){
  #       loglik[i] <-  sum(-qy[i]*log(1 + exp(-(hx[i]))) - (1-qy[i])*log(1 + exp(hx[i]))) + 
  #         qrLogLik(y[i], mu[i], sigma[i], fd, sd)
  #     }else{
  #       loglik[i] <- sum(-qy[i]*log(1 + exp(-(hx[i]))) - (1-qy[i])*log(1 + exp(hx[i])))
  #     }
  #   }
  #   -sum(loglik, na.rm = TRUE)
  # }
  
  if(type =="ZI"){
    qy <- qy1;
    n_zo <- n_zm
    ycdf <-  ydata[qy1 == 0]
    xmat = data_lm[qy1==0, ,drop=FALSE]
    zmat = data_pm[qy1==0, ,drop=FALSE]
    
    data_zo <- data.frame(data_zero)
    data_zo$qy <- qy
    zeroone.fo <- paste0(c("qy", as.character(zero.fo)), collapse = "")
    zeroone.fo <- Formula::as.Formula(zeroone.fo) 
  }
  
  if(type =="OI"){
    qy <- qy2;
    n_zo <- n_om
    ycdf <-  ydata[qy2 == 0]
    xmat = data_lm[qy2== 0, ,drop=FALSE]
    zmat = data_pm[qy2== 0, ,drop=FALSE]
    
    data_zo <- data.frame(data_one)
    data_zo$qy <-  qy
    zeroone.fo <- paste0(c("qy", as.character(one.fo)), collapse = "")
    zeroone.fo <- Formula::as.Formula(zeroone.fo) 
  }
  
  if(type =="ZO"){
    ycdf <-  ydata[qy1==0 & qy2== 0]
    xmat = data_lm[qy1==0 & qy2== 0, ,drop=FALSE]
    zmat = data_pm[qy1==0 & qy2== 0, ,drop=FALSE]
    
    data_one <- data.frame(data_one)
    data_one$qy <-  qy2
    one.fo <- paste0(c("qy", as.character(one.fo)), collapse = "")
    one.fo <- Formula::as.Formula(one.fo) 
    
    data_zero <- data.frame(data_zero)
    data_zero$qy <-  qy1
    zero.fo <- paste0(c("qy", as.character(zero.fo)), collapse = "")
    zero.fo <- Formula::as.Formula(zero.fo) 
  }

  # starting values
  if(is.null(start)){
    # if the user does not supply starting values, generate some default values
    start <- rep(0.1, k)
    start_0 <- qrStart(ycdf, fd=fd, sd=sd)
    start[n_zm + 1] <- start_0[1] #starting value for mu intercept
    start[n_zm + n_lm + 1] <- start_0[2] #starting value for sigma intercept
  }
  
  # Esimating First approach-
  # if(type == "ZI" | type == "OI"){
  #   preopt <- optim(start, fn = zeroonepart, hessian = F, x = data_lm, y = ydata, 
  #                   z = data_pm, w = data_pm, qy = qy, fd=fd, sd=sd,
  #                   method = "Nelder-Mead")
  #   
  #   contr.method <- control$method
  #   #contr.hessian <- control$hessian
  #   contr.hessian <- TRUE
  #   contr.trace <- control$trace
  #   contr.maxit <- control$maxit
  #   
  #   betaopt <- optim(preopt$par, fn = zeroonepart,x = data_lm, y = ydata, 
  #                    z = data_pm, w = data_pm, qy = qy, fd=fd, sd=sd,
  #                    hessian = contr.hessian, method = contr.method, 
  #                    control=list(trace = contr.trace , maxit = contr.maxit))
  #   logLiklihod <- -betaopt$value
  # }

  # Estimation - Two stage apporach---------
  if(type == "ZI" | type == "OI"){
    data$qy = qy
    #fit zero/one part--
    ZOpart <- glm(zeroone.fo, data = data, family = binomial)
    
    #extract estimated values
    zstatzo <- coef(ZOpart)/sqrt(diag(vcov(ZOpart)))
    est_zo <-  cbind(estim <-  coef(ZOpart), 
                   serr <- sqrt(diag(vcov(ZOpart))), 
                   zstat <- zstatzo,
                   prob <- 2 * (1 - pnorm(abs(zstatzo))))
    colnames(est_zo) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    
    #extract fitted value and residuals
    fitted_zop <- ZOpart$fitted.values
    residual_zop <- ZOpart$residuals
    deviance_zop <- ZOpart$deviance
    vcov_zop <- vcov(ZOpart)
    if(type == "ZI"){
      est_zm <- est_zo
      est_om <-NULL
      fitted_zp <- fitted_zop
      fitted_op <- NULL
      residual_zp <- residual_zop
      residual_op <- 0
      deviance_zp <- deviance_zop
      deviance_op <- 0
      vcov_zp <- vcov_zop
      vcov_op <- NULL
      mod_zp <- ZOpart
      mod_op <- NULL
    }else{
      est_zm <- NULL
      est_om <- est_zo
      fitted_zp <- NULL
      fitted_op <- fitted_zop
      residual_zp <- 0
      residual_op <- fitted_zop
      deviance_zp <- 0
      deviance_op <- fitted_zop
      vcov_zp <- NULL
      vcov_op <- vcov_zop
      mod_zp <- NULL
      mod_op <- ZOpart
    }
  
  
    
    #fit cdfquantreg part--
    preopt <- optim(start[-c(1:n_zo)], fn = quantreg_loglik, hessian = F, 
                    x = xmat, z = zmat, y = ycdf, method = "Nelder-Mead")
    
    contr.method <- control$method
    #contr.hessian <- control$hessian
    contr.hessian <- TRUE
    contr.trace <- control$trace
    contr.maxit <- control$maxit
    
    betaopt <- optim(preopt$par, fn = quantreg_loglik,x = xmat, z = zmat, y = ycdf,
                     hessian = contr.hessian, method = contr.method, 
                     control=list(trace = contr.trace , maxit = contr.maxit))
    
    logLiklihod <- -betaopt$value + as.numeric(logLik(ZOpart))
    }
  
  if(type == "ZO"){
    data$qy <- qy1
    Zeropart <- glm(zero.fo, data = data, family = binomial)
    
    #extract estimated values
    zstatz <- coef(Zeropart)/sqrt(diag(vcov(Zeropart)))
    est_zm <-  cbind(estim <-  coef(Zeropart), 
                     serr <- sqrt(diag(vcov(Zeropart))), 
                     zstat <- zstatz,
                     prob <- 2 * (1 - pnorm(abs(zstatz))))
    colnames(est_zm) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    
    data$qy = qy2
    Onepart <- glm(one.fo, data = data, family = binomial)
    zstato <- coef(Onepart)/sqrt(diag(vcov(Onepart)))
    est_om <-  cbind(estim <-  coef(Onepart), 
                     serr <- sqrt(diag(vcov(Onepart))), 
                     zstat <- zstato,
                     prob <- 2 * (1 - pnorm(abs(zstato))))
    colnames(est_om) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    
    #extract fitted value and residuals
    fitted_zp <- Zeropart$fitted.values
    residual_zp <- Zeropart$residuals
    deviance_zp <- Zeropart$deviance
    vcov_zp <- vcov(Zeropart)
    mod_zp <- Zeropart
    
    #extract fitted value and residuals
    fitted_op <- Onepart$fitted.values
    residual_op <- Onepart$residuals
    deviance_op <- Onepart$deviance
    vcov_op <- vcov(Onepart)
    mod_op <- Onepart
    
    preopt <- optim(start[-c(1:(n_zm+n_om))], fn = quantreg_loglik, hessian = F, 
                    x = xmat, 
                    y = ycdf, 
                    z = zmat, 
                    method = "Nelder-Mead")
    
    contr.method <- control$method
    #contr.hessian <- control$hessian
    contr.hessian <- TRUE
    contr.trace <- control$trace
    contr.maxit <- control$maxit
    
    betaopt <- optim(preopt$par, fn = quantreg_loglik,
                     x = xmat, 
                     y = ycdf, 
                     z = zmat, 
                     hessian = contr.hessian, method = contr.method, 
                     control=list(trace = contr.trace , maxit = contr.maxit))
    
    # ***************************** 
    logLiklihod <- -betaopt$value + logLik(Zeropart)+ logLik(Onepart)
  }
  
 gradients <- grad(betaopt$par,ycdf,xmat,zmat)
  # Convergence message from optim
  converge <- betaopt$convergence
  
  
  # ***************************** Obtain other parameter specifications, get the
  # distribution specific loglikelihood function and grad
  
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
  

  # Generall fit-------------------- Loglike
  
  attr(logLiklihod, "nall") <- n
  attr(logLiklihod, "nobs") <- n
  attr(logLiklihod, "df") <- k
  class(logLiklihod) <- "logLik"
  
  # AIC & BIC
  AIC <- AIC(logLiklihod)
  BIC <- BIC(logLiklihod)
  
  # Fitted values for cdfquant part-------------
  if(n_lm == 1) {fitted_mu <- xmat * est_lm[, 1]
  }else{
    fitted_mu <- rowSums(t(apply(xmat, 1, 
                                 function(x) x * est_lm[, 1])))
  }
  if  (fd =="km"){
    #log link for parameter a in Km distribution
    fitted_mu <- exp(fitted_mu)
  }
  
  #dispersion model
  if(n_pm == 1) {fitted_phi <- zmat * est_pm[, 1]
  }else{
    fitted_phi <- rowSums(t(apply(zmat, 1, function(x) x * est_pm[, 1])))
  }
  
  fitted_sigma <- exp(fitted_phi)
  
  
  q_y <- seq(1e-06, 0.999999, length.out = n)
  cdf_y <- q_y[rank(ycdf)]
  
  fitted <- qq(cdf_y, fitted_mu, fitted_sigma, fd, sd)
  
  
  # Raw residuals
  residuals <- ycdf - fitted
  df.residual <- n - k
  
  # Deviance
  deviance <- 2 * (qrLogLik(ycdf, fitted_mu, fitted_sigma, fd, sd) - 
                     qrLogLik(fitted,  fitted_mu, fitted_sigma, fd, sd))
  
  

  # Parameter variance-covariance matrix
  vcov <- solve(as.matrix(betaopt$hessian))
  rownames(vcov) <- colnames(vcov) <- c(rownames(est_lm), 
                                        paste('(sigma)_', rownames(est_pm), sep = ""))
  
  # RMSE
  rmse <- sqrt(mean(residuals^2))
  
  rmseLogit <- NULL
  # RMSlogit(E)
  for (i in 1:n){
    if(is.nan(fitted[i])|is.na(fitted[i])){
      rmseLogittemp <- NA
    }else{
      rmseLogittemp <-sqrt(mean((logit(fitted[i]) - logit(ycdf[i]))^2))
    }
    
    rmseLogit <- c(rmseLogit, rmseLogittemp)
  }
  
 
  
  #Combine zero/one and cdf part output
  # Coefficients:
  coefficients <- list(zero = est_zm,
                       one = est_om,
                       location = est_lm, 
                       dispersion = est_pm)
  
  # Fitted values
  fittemp <- rep(NA, n)
  fittemp[qy1 == 0 & qy2 == 0] <- fitted
  fitted <- fittemp
  fittemp[qy1 == 0 & qy2 == 0] <- fitted_mu
  fitted_mu <- fittemp
  fittemp[qy1 == 0 & qy2 == 0] <- fitted_sigma
  fitted_sigma <- fittemp
  
  fittedlist <- switch(fd,
                        list(full = fitted, zerop = fitted_zp, onep = fitted_op, 
                             mu = fitted_mu, sigma = fitted_sigma),
                        km = list(full = fitted, zerop = fitted_zp, onep = fitted_op,
                                  a = fitted_mu, b = fitted_sigma))
  
  #Residuals
  redtemp <-rep(NA, n) 
  redtemp[qy1 == 0 & qy2 == 0] <- residuals
  residuals <- redtemp
  residuals <- list(rcdf = residuals,
                    rzero = residual_zp,
                    rone = residual_op)
  
  #variance-covariance
  vcov <- list(cdf = vcov,
               zero = vcov_zp,
               one = vcov_op)
  
  deviance <- deviance_zp + deviance_op + deviance
  
  # Output
  output <- list(family = list(fd = fd, sd = sd), type = type,call = input, 
                 formula = formula, zero.fo = zero.fo, one.fo = one.fo,
                 N = n, k = k, y = ydata, 
                 coefficients = coefficients, 
                 residuals = residuals, 
                 fitted = fittedlist, 
                 vcov = vcov, 
                 logLik = logLiklihod, 
                 deviance = deviance,
                 aic = AIC, bic = BIC,
                 mod_zp = mod_zp,
                 mod_op = mod_op,
                 df.residual = df.residual)
  
  class(output) <- "cdfqrH"
  
  rm(grad)
  #print.cdfqr(output)
  return(output)
  
} 
