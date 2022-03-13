#' @title CDF-Quantile Finite Tailed Probability Distributions
#' @aliases cdfquantregFT
#' @description \code{cdfquantregFT} is a function to fit a cdf quantile regression with a variety of finite tailed distributions. It can account for data that has boundary values. 
#' 
#' @param formula A formula object, with the dependent variable (DV) on the left of an ~ operator, and predictors on the right. For the part on the right of '~', the specification of the dispersion (sigma; first) and skewness (theta; second) submodels can be separated by '|'. So \code{y ~ X1 | X2} specifies that the DV is \code{y}, \code{X1} is the predictor in the dispersion submodel, and \code{X2} is the predictor in the skewness submodel.
#' @param data The data in a data.frame format 
#' @param fd A string that specifies the parent distribution. At the moment, only "arcsinh", "cauchit" and "t2" can be used. See details. 
#' @param sd A string that specifies the child distribution. At the moment, only "arcsinh", "cauchy" and "t2" can be used. See details. 
#' @param mu.fo A formula object to indicate the predictors for the location submodel if the 3-parameter distribution is used, only input as \code{~ predictors} 
#' @param inner A logic value that indicates if the inner (\code{inner = TRUE}) case or outer (\code{inner = FALSE}) will be used. Currently inner case can only be used for 2-parameter distributions.
#' @param version A string indicates that which version will be used. "V" is the tilt transformation while "W" indicates the Jones Pewsey transformation. 
#' @param family If `fd` and `sd` are not provided, the name of a member of the family of distributions can be provided (see below) for details of family functions)
#' @param start The starting values for model fitting. If not provided, default values will be used.
#' @param ssn The number of searches on optimal starting values to be performed. If model does not converge, can increase this number.
#' @param control Control optimization parameters  (See \code{\link{cdfqr.control}}))
#' @param ... Currently ignored. 
#'
#' @details The cdfquantregFT function fits a quantile regression model with a distributions from the cdf-quantile finite tailed distributions. 
#' Here is the list of currently available distributions.
#' 
#' \bold{Bimodal Shape Distributions}
#'   \tabular{lllcc}{
#'   \bold{Distribution}  \tab \bold{R input} \tab \bold{Alternative Input}  \tab \bold{Available Version}\cr
#'   ArcSinh-ArcSinh	\tab \code{fd = "arcsinh", sd = "arcsinh"}	\tab \code{family = "arcsinh-arcsinh"} \tab	\code{"V", "W"}	\cr
#'   ArcSinh-Cauchy	\tab \code{fd = "arcsinh", sd = "cauchy"}	\tab \code{family = "arcsinh-cauchy"} \tab	\code{"V", "W"}	\cr
#'   Cauchit-ArcSinh	\tab \code{fd = "cauchit", sd = "arcsinh"}	\tab \code{family = "cauchit-arcsinh"} \tab	\code{"V", "W"}	\cr
#'   Cauchit-Cauchy	\tab \code{fd = "cauchit", sd = "cauchy"}	\tab \code{family = "cauchit-cauchy"} \tab	\code{"V", "W"}	\cr
#'   T2-T2	\tab \code{fd = "t2", sd = "t2"}	\tab \code{family = "t2-cauchy"} \tab	\code{"V", "W"}	\cr
#'   } 
#' @export
#' @import Formula
#' 
#' @return An object of class \code{cdfqrFT} will be returned. Generic functions such as \link{summary},\link{print} and \link{coef} can be used to extract output (see \link{summary.cdfqr} for more details about the generic functions that can be used).
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
#' fit <- cdfquantregFT(pnurse ~ Ambulance |Ambulance ,
#'  fd = "arcsinh", sd = "arcsinh",  inner = FALSE, version = "V", data = yoon)
#' summary(fit)

cdfquantregFT <- function(formula,
                          fd = NULL,
                          sd = NULL,
                          mu.fo = NULL,
                          inner = FALSE,
                          version = "V",
                          data,
                          family = NULL,
                          start = NULL, 
                          ssn = 20, 
                          control = cdfqr.control(...),
                          ...) {
  # ***************************** Diagnose the input, extract input values
  input_dist <- mget(names(formals()),sys.frame(sys.nframe()))
  
  input_dist <- input_dist[c(
    "fd",
    "sd",
    "family",
    "inner",
    "version"
  )]
  
  input <- match.call(expand.dots = FALSE)
  input_index <- match(
    c(
      "formula",
      "data",
      "mu.fo",
      "fd",
      "sd",
      "family",
      "inner",
      "version",
      "start"
    ),
    names(input),
    0L
  )
  input <- input[c(1L, input_index)]
  
  if (is.null(family) & (is.null(fd) | is.null(sd))) {
    stop("No distribution is specified!")
  }
  
  if (!is.null(family)) {
    fd <- strsplit(family, "-")[[1]][1]
    sd <- strsplit(family, "-")[[1]][2]
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
  
  if (any(ydata < 0) | any(ydata > 1)) {
    warning(
      paste(
        "The values of the dependent variable is rescaled into [0, 1] interval",
        "\n",
        "see scaleTR() for details"
      )
    )
    ydatatemp = rep(0, length(ydata))
    ydatatemp[ydata == min(ydata)] <- 0
    ydatatemp[ydata == max(ydata)] <- 1
    ydatascale <- scaleTR(ydata)
    ydatatemp[ydata < max(ydata) & ydata > min(ydata)] <-
      ydatascale[ydata < max(ydata) & ydata > min(ydata)]
    ydata <- ydatatemp
    
  }
  
  
  data_tm <- as.matrix(model.matrix(formula, data, rhs = 2))  # datamatrix for theta model
  data_sm <- as.matrix(model.matrix(formula, data, rhs = 1))  # datamatrix for sigma model
  
  n_sm <- ncol(data_sm)  #number of parameters in sigma component
  n_tm <- ncol(data_tm)  #number of parameters in theta component
  
  n <- nrow(data)  #number of responses
  k <- n_sm + n_tm  # total number of parameters
  
  if (!is.null(mu.fo)) {
    mu.fo <- Formula::as.Formula(mu.fo)
    data_mu <- as.matrix(model.matrix(mu.fo, data, rhs = 1))
    n_zm <- ncol(data_mu)  #number of parameters in mu  component
  } else{
    n_zm = 0
    data_mu = NULL
  }
  
  
  k <- n_sm + n_tm + n_zm
  
  
  # ***************************** Obtain other parameter specifications, get the
  # starting values
  
  if (is.null(start)) {
    # if the user does not supply starting values, generate some default values
    start_0 <- qrStart(ydata, fd=fd, sd=sd, skew = TRUE)
    start_s <- start_0[2] #starting value for sigma intercept
    if(n_sm>1) start_s <- c(start_s, rep(0.1, n_sm-1))
    start_t <- start_0[3]
    if(n_tm>1) start_t <- c(start_t, rep(sign(start_t)*(-0.1), n_tm-1))
 
    if(!is.null(mu.fo)) {
      start_m = start_0[1]
      if(n_zm>1) start_m <- c(start_m, rep(0.1, n_zm-1))
      
    }else{
      start_m = NULL
    }
    
    start1 <- c(start_s, start_t, start_m)

    test = tryCatch({
      MLE1 <-  optim(
        start1,
        .qrLogLikFunFT,
        hessian = F,
        sdata = data_sm,
        y = ydata,
        tdata = data_tm,
        fd = fd,
        sd = sd,
        mdata = data_mu,
        inner = inner,
        version = version,
        method = "BFGS"
      )
    }, error=function(e) e)
    
    
    if ("error" %in% class(test)) {

      if(n_zm+n_sm+n_tm > 3){
        start0 <- c(start_s[1], start_t[1], start_m[1])
        
        if(is.null(data_mu)){
          test0 = tryCatch({
            MLE1 <-  optim(
              start0,
              .qrLogLikFunFT,
              hessian = F,
              x = data_sm[, 1],
              y = ydata,
              z = data_tm[, 1],
              fd = fd,
              sd = sd,
              w = data_mu,
              inner = inner,
              version = version,
              method = "BFGS"
            )
          }, error=function(e) e)
        } else {
          test0 = tryCatch({
            MLE1 <-  optim(
              start0,
              .qrLogLikFunFT,
              hessian = F,
              sdata = data_sm[, 1],
              y = ydata,
              tdata = data_tm[, 1],
              fd = fd,
              sd = sd,
              mdata = data_mu[, 1],
              inner = inner,
              version = version,
              method = "BFGS"
            )
          }, error=function(e) e)
        }

        if (!("error" %in% class(test0))){
          for (i in 1:ssn){
          if(n_sm>1) start_s <- c(test$par[2], runif(n_sm-1, -0.1, 0.1))
          if(n_tm>1) start_t <- c(test$par[1], runif(n_tm-1, -0.1, 0.1))
          if(n_zm>1) start_t <- c(test$par[3], runif(n_zm-1, -0.1, 0.1))
          
          start1 <- c(start_s, start_t, start_m)
          test = tryCatch({
            MLE1 <-  optim(
              start1,
              .qrLogLikFunFT,
              hessian = F,
              sdata = data_sm,
              y = ydata,
              tdata = data_tm,
              fd = fd,
              sd = sd,
              mdata = data_mu,
              inner = inner,
              version = version,
              method = "BFGS"
            )
          }, error=function(e) e)
          
          if (!("error" %in% class(test))) {
            break
          }
        }
        
        }
        
      }
      
      
    }
    
    
    if ("error" %in% class(test)) {
     
      for (i in 1:ssn){
        if(!is.null(mu.fo)) {
          start_m <- runif(1, -1, 1)
          if(n_zm > 1) start_m <- c(start_m, runif(n_zm-1, -0.1, 0.1))
          
        }else{
          start_m = NULL
        }
        
        start_s <- runif(1, -1, 1)
        if(n_sm > 1) start_s <- c(start_s, runif(n_sm-1, -0.1, 0.1))
        
        start_t <- runif(1, -1, 1)
        if(n_tm > 1) start_t <- c(start_t, runif(n_tm-1, -0.1, 0.1))
        
        start1 <- c(start_s, start_t, start_m)
        test = tryCatch({
          MLE1 <-  optim(
            start1,
            .qrLogLikFunFT,
            hessian = F,
            sdata = data_sm,
            y = ydata,
            tdata = data_tm,
            fd = fd,
            sd = sd,
            mdata = data_mu,
            inner = inner,
            version = version,
            method = "BFGS"
          )
        }, error=function(e) e)
        if (!("error" %in% class(test))) {
         break
        }
        
        
        
      }
      if("error" %in% class(test)){start1 = start1*0.1}  else {start1 <- test$par}
    } 
    
    
    
    
    
  } else {start1 = start}


  
  contr.method <- control$method
  contr.hessian <- TRUE
  contr.trace <- control$trace
  contr.maxit <- control$maxit
  
  betaopt <-
    optim(
      start1,
      fn = .qrLogLikFunFT,
      sdata = data_sm,
      y = ydata,
      tdata = data_tm,
      fd = fd,
      sd = sd,
      mdata = data_mu,
      inner = inner,
      version = version,
      hessian = contr.hessian,
      method = contr.method,
      control = list(trace = contr.trace , maxit = contr.maxit)
    )
  
  
  # ***************************** Process the output parameter
  # Convergence message from optim
  converge <- betaopt$convergence
  
  # estimation--------------------
  estim <- betaopt$par  #mean estiamte
  serr <- sqrt(sqrt(diag(solve(betaopt$hessian)) ^ 2))  # se
  zstat <- estim / serr  # z-test
  prob <- 2 * (1 - pnorm(abs(zstat)))  # p value
  
  est <- cbind(estim <- estim,
               serr <- serr,
               zstat <- estim / serr,
               prob <- 2 * (1 - pnorm(abs(zstat))))
  colnames(est) <-
    c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  
  # Separate output depending on the submodels
  est_sm <- est[1:n_sm,]  #sigma model
  if (n_sm == 1)
    est_sm <- matrix(est_sm, ncol = 4)
  
  rownames(est_sm) <- colnames(data_sm)
  colnames(est_sm) <- colnames(est)
  
  est_tm <- est[(1 + n_sm):(n_sm + n_tm),]   #theta model
  if (n_tm == 1)
    est_tm <- matrix(est_tm, ncol = 4)
  rownames(est_tm) <- colnames(data_tm)
  colnames(est_tm) <- colnames(est)
  
  if (!is.null(mu.fo)) {
    est_zm <- est[(1 + n_sm + n_tm):nrow(est),]   #mu model
    if (n_zm == 1)
      est_zm <- matrix(est_zm, ncol = 4)
    rownames(est_zm) <- colnames(data_mu)
    colnames(est_zm) <- colnames(est)
    coefficients <- list(dispersion = est_sm,
                         skew = est_tm,
                         location = est_zm)
  } else{
    coefficients <- list(dispersion = est_sm,
                         skew = est_tm,
                         location = NULL)
  }
  
  
  # General fit-------------------- Loglike
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
  if (n_sm == 1) {
    fitted_sigma <- data_sm * est_sm[, 1]
  } else{
    fitted_sigma <- rowSums(t(apply(data_sm, 1,
                                    function(x)
                                      x * est_sm[, 1])))
  }
  fitted_sigma = exp(fitted_sigma)
  #Theta model
  if (n_tm == 1) {
    fitted_theta <- data_tm * est_tm[, 1]
  } else{
    fitted_theta <-
      rowSums(t(apply(data_tm, 1, function(x)
        x * est_tm[, 1])))
  }
  
  #mu model
  if (!is.null(mu.fo)) {
    #Theta model
    if (n_zm == 1) {
      fitted_mu <- data_mu * est_zm[, 1]
    } else{
      fitted_mu <-
        rowSums(t(apply(data_mu, 1, function(x)
          x * est_zm[, 1])))
    }
    
  } else{
    fitted_mu = NULL
  }
  
  
  q_y <- seq(1e-06, 0.999999, length.out = n)
  cdf_y <- q_y[rank(ydata)]
  
  fitted <-
    cdfft(
      cdf_y,
      sigma = fitted_sigma,
      theta = fitted_theta,
      fd = fd,
      sd = sd,
      mu = fitted_mu,
      inner = inner,
      version = version
    )
  
  # Raw residuals
  residuals <- ydata - fitted
  df.residual <- n - k
  
  # Deviance
  loglik1 <-  -sum(log(
    pdfft(
      as.numeric(ydata),
      sigma = as.numeric(fitted_sigma),
      theta = as.numeric(fitted_theta),
      fd = fd,
      sd = sd,
      mu =  as.numeric(fitted_mu),
      inner = inner,
      version = version
    )
  ), na.rm = TRUE)
  loglik2 <-
    -sum(log(
      pdfft(
        as.numeric(fitted) ,
        sigma = as.numeric(fitted_sigma),
        theta = as.numeric(fitted_theta),
        fd = fd,
        sd = sd,
        mu = as.numeric(fitted_mu),
        inner = inner,
        version = version
      )
    ), na.rm = TRUE)
  
  deviance <- 2 * (loglik1 - loglik2)
  
  # Parameter variance-covariance matrix
  vcov <- solve(as.matrix(betaopt$hessian))
  if (!is.null(mu.fo)) {
    rownames(vcov) <-
      colnames(vcov) <- c(
        paste('(sigma)_', rownames(est_sm), sep = ""),
        paste('(theta)_', rownames(est_tm), sep = ""),
        paste('(mu)_', rownames(est_zm), sep = "")
      )
  } else{
    rownames(vcov) <-
      colnames(vcov) <- c(paste('(sigma)_', rownames(est_sm), sep = ""),
                          paste('(theta)_', rownames(est_tm), sep = ""))
  }
  
  # RMSE
  rmse <- sqrt(mean(residuals ^ 2))
  
  rmseLogit <- NULL
  for (i in 1:n) {
    if (is.nan(fitted[i]) | is.na(fitted[i])) {
      rmseLogittemp <- NA
    } else{
      rmseLogittemp <- sqrt(mean((logit(fitted[i]) - logit(ydata[i])) ^ 2))
    }
    
    rmseLogit <- c(rmseLogit, rmseLogittemp)
  }
  
  fitted = list(
    full = fitted,
    sigma = fitted_sigma,
    theta = fitted_theta,
    mu = fitted_mu
  )
  
  # Output
  output <-
    list(
      family = list(
        fd = fd,
        sd = sd,
        inner = inner,
        version = version
      ),
      coefficients = coefficients,
      call = input,
      distinf = input_dist,
      formula = formula,
      N = n,
      k = k,
      y = ydata,
      residuals = residuals,
      fitted = fitted,
      vcov = vcov,
      rmse = rmse,
      rmseLogit = rmseLogit,
      logLik = logLiklihod,
      deviance = deviance,
      aic = AIC,
      bic = BIC,
      converged = converge,
      grad = NULL,
      df.residual = df.residual,
      optim = betaopt
    )
  
  class(output) <- c("cdfqrFT","cdfqr")
  
  return(output)
  
} 
