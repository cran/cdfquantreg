#' @title S3 Methods for getting output from fitted cdfqrr Objects.
#' @description Give S3 Methods for CDF-Quantile Distribution Models
#' @param object The fitted cdfqrH model.
#' @param x The fitted cdfqrH model.
#' @param type The parts of coefficients or variance-covariance matrix to be extracted.Can be "full", "mean",or "sigma". 
#' @param formula. Changes to the formula. See \code{\link[Formula]{update.Formula}} for details.
#' @param evaluate If true evaluate the new updated model else return the call for the new model.
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level the confidence level required.
#' @param digits Number of digits to be retained in printed output.
#' @param zero.fo. Changes to the formula for the zero component, only input as \code{~ predictors}. See \code{\link[Formula]{update.Formula}} for details. 
#' @param one.fo. Changes to the formula for the one component, only input as \code{~ predictors}. See \code{\link[Formula]{update.Formula}} for details. 
#' @param ... Pass onto other functions or currently ignored
#' @examples
#' data(cdfqrHExampleData)
#' ipcc_mid <- subset(IPCC, mid == 1 & high == 0)
#' fit <- cdfquantregH(prob ~ valence | valence, zero.fo = ~valence,
#'   one.fo = ~valence,
#'   fd ='t2',sd ='t2', type = "ZO", data = ipcc_mid)
#' 
#' summary(fit)
#' print(fit)
#' logLik(fit)
#' coef(fit)
#' deviance(fit)
#' vcov(fit)
#' confint(fit)
#' 
#' 
#' @method summary cdfqrH
#' @export
summary.cdfqrH <- function(object, ...) {
  print(object)
  #return(object)
}

#' @method print cdfqrH
#' @export
#' @rdname summary.cdfqrH
print.cdfqrH <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  
  cat("Family: ", paste(x$family$fd, x$family$sd, collapse = "-"), "\n")
  cat("Call: ", deparse(x$call), "\n", fill = TRUE)
  
  if (length(x$coefficients$location)) {
    
    mutitile <- "Mu coefficients (Location submodel)\n"
    if (x$family$fd == 'km'){
      mutitile <- "1st shape parameter coefficients (log link)\n"
    }
    cat(mutitile)
    printCoefmat(x$coefficients$location, digits=digits,signif.legend = FALSE)
    cat("\n")
  } else{
    mutitile2 <- "No coefficients in location submodel\n\n"
    if (x$family$fd == 'km'){
      mutitile2 <- "No coefficient for the 1st shape parameter (log link)\n"
    }
    cat(mutitile2)
  }
  
  if (length(x$coefficients$dispersion)) {
    if (length(x$coefficients$dispersion)) {
      Sigmatitile <- "Sigma coefficients (Dispersion submodel)\n"
      if (x$family$fd == 'km'){
        Sigmatitile <- "2nd shape parameter coefficients (log link)\n"
      }
      cat(Sigmatitile)
      printCoefmat(x$coefficients$dispersion, digits=digits)
      cat("\n")
    } else{
      Sigmatitile2 <- "No coefficients in Dispersion submodel\n"
      if (x$family$fd == 'km'){
        Sigmatitile2 <- "No coefficient for the 2nd shape parameter (log link)\n"
      }
      cat(Sigmatitile2)
    }
  }
  
  if (length(x$coefficients$zero)) {
    
    mutitile <- "Zero component coefficients\n"
    cat(mutitile)
    printCoefmat(x$coefficients$zero, digits=digits,signif.legend = FALSE)
    cat("\n")
  } 
  if (length(x$coefficients$one)) {
    
    mutitile <- "One component coefficients\n"
    cat(mutitile)
    printCoefmat(x$coefficients$one, digits=digits,signif.legend = FALSE)
    cat("\n")
  } 
  
  cat("Log-Likelihood: ", round(x$logLik, digits=digits), "\n", fill = TRUE)
  invisible(x)
  
}

#' @method logLik cdfqrH
#' @export
#' @rdname summary.cdfqrH
logLik.cdfqrH <- function(object, ...) {
  object$logLik
  
}

#' @method nobs cdfqrH
#' @export
#' @rdname summary.cdfqrH
nobs.cdfqrH <- function(object, ...) {
  length(object$residuals)
}

#' @method deviance cdfqrH
#' @export
#' @rdname summary.cdfqrH
deviance.cdfqrH <- function(object, ...) {
  object$deviance
}

#' @method coef cdfqrH
#' @export
#' @rdname summary.cdfqrH
coef.cdfqrH <- function(object, type = c("full","mean","sigma","zero","one"), ...) {
  type <- match.arg(type)
  mean <- object$coefficients$location[,1]
  sigma <- object$coefficients$dispersion[,1]
  if(!is.null(object$coefficients$zero)) {
    zero <- object$coefficients$zero[,1]
    names(zero) <- paste('(zero)_',names(zero),sep = "")
  }else{zero <- NULL}
   
  
  if(!is.null(object$coefficients$one)) {
    one <- object$coefficients$one[,1]
    names(one) <- paste('(one)_',names(one),sep = "")
  }else{one <- NULL}
    
  
  names(sigma) <- paste('(sigma)_',names(sigma),sep = "")
 
  coef <- switch(type, 
               full = {c(mean,sigma, zero, one)}, 
               mean = {mean}, 
               sigma = {sigma}, 
               zero = {zero}, 
               one = {one})
  return(coef)
}

#' @method vcov cdfqrH
#' @export
#' @rdname summary.cdfqrH
vcov.cdfqrH <- function(object, type = c("full","mean","sigma","zero","one"), ...) {
  
  vccdf <- object$vcov$cdf
  zero <- object$vcov$zero
  one <- object$vcov$one
  k <- nrow(object$coefficients$location)
  m <- nrow(object$coefficients$dispersion)
  
  type = match.arg(type)
  
  mean <- vccdf[seq.int(length.out = k), seq.int(length.out = k), drop = FALSE]
  
  dispersion <- vccdf[seq.int(length.out = m) + k, seq.int(length.out = m) + k, drop = FALSE]

  
  vc <- switch(type, full = {object$vcov}, 
               mean = {mean}, 
               sigma = {dispersion}, 
               zero = {zero}, 
               one = {one})
  return(vc)
}

#' @method update cdfqrH
#' @export
#' @rdname summary.cdfqrH
update.cdfqrH <- function(object, formula., zero.fo., 
                          one.fo., ..., evaluate = TRUE) {
  call <- object$call
  if (!missing(formula.)) 
    call$formula <- formula(update(Formula::Formula(formula(call$formul)), formula.))
  
  if (!missing(zero.fo.)) 
    call$zero.fo <- formula(update(Formula::Formula(formula(call$zero.fo)), zero.fo.))
  
  if (!missing(one.fo.)) 
    call$one.fo <- formula(update(Formula::Formula(formula(call$one.fo)), one.fo.))
  
  extras <- match.call(expand.dots = FALSE)$...

  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if (evaluate) 
    eval(call, parent.frame())
  else call
} 

#' @method formula cdfqrH
#' @export
#' @rdname summary.cdfqrH
formula.cdfqrH <- function(x, ...) {
  call <- x$call
  return(formula(call$formula))
} 


#' @method confint cdfqrH
#' @export
#' @rdname summary.cdfqrH
confint.cdfqrH <- function(object, parm, level = 0.95, 
                           type = c("full","mean","sigma","zero","one"), ...) {
 
  type = match.arg(type)
  
  cf_full <- object$coefficients
  cf_mean <- cf_full['location'][[1]]
  rownames(cf_mean) <- paste( rownames(cf_mean), "_m", sep = "")
 
   cf_sigma <- cf_full['dispersion'][[1]]
  rownames(cf_sigma) <- paste( rownames(cf_sigma), "_d", sep = "")
  
  cf_zero <- cf_full['zero'][[1]]
  rownames(cf_zero) <- paste( rownames(cf_zero), "_z", sep = "")
  
  cf_one <- cf_full['one'][[1]]
  rownames(cf_one) <- paste( rownames(cf_one), "_0", sep = "")
  
  #cf <- do.call(rbind, object$coefficients)
  cf <- rbind(cf_mean, cf_sigma, cf_zero, cf_one)
  if(missing(parm)) parm = rownames(cf)
    
  if (type == "mean") {
      if (is.character(parm)){
        if (!all(parm %in% rownames(cf_mean)))
        stop("One or more required parameter(s) is not in the location submodel!")}
      
      if (is.numeric(parm)){
        if (!is.matrix(try(cf_mean[parm, ],T)))
          stop("One or more required parameter(s) is not in the location submodel!")}
      
    cf <- cf_mean<- cf_mean[parm, ,drop=FALSE]
    }
    
  if (type == "sigma") {
    if (is.character(parm)){
      if (!all(parm %in% rownames(cf_sigma)))
        stop("One or more required parameter(s) is not in the dispersion submodel!")}
      
      if (is.numeric(parm)){
        if (!is.matrix(try(cf_sigma[parm, ],T)))
          stop("One or more required parameter(s) is not in the dispersion submodel!")}
      
     cf <- cf_sigma<- cf_sigma[parm, ,drop=FALSE]
    }
  
  if (type == "one") {
    if (is.character(parm)){
      if (!all(parm %in% rownames(cf_one)))
        stop("One or more required parameter(s) is not in the one-component submodel!")}
    
    if (is.numeric(parm)){
      if (!is.matrix(try(cf_one[parm, ],T)))
        stop("One or more required parameter(s) is not in the one-component submodel!")}
    
    cf <- cf_one<- cf_one[parm, ,drop=FALSE]
  }
  
  if (type == "zero") {
    if (is.character(parm)){
      if (!all(parm %in% rownames(cf_zero)))
        stop("One or more required parameter(s) is not in the zero-component submodel!")}
    
    if (is.numeric(parm)){
      if (!is.matrix(try(cf_zero[parm, ],T)))
        stop("One or more required parameter(s) is not in the zero-component submodel!")}
    
    cf <- cf_zero<- cf_zero[parm, ,drop=FALSE]
  }
  if (type == "full") {
    if (is.character(parm)){
      if (!any(parm %in% rownames(cf)))
        stop("One or more required parameter(s) is not in the model!")
      }
      
      if (is.numeric(parm)){
        if (!is.matrix(try(cf[parm, ],T)))
          stop("One or more required parameter(s) is not in the model!")}
      
    cf_mean <- cf_mean[parm[parm%in%rownames(cf_mean)], ,drop=FALSE]
    cf_sigma <- cf_sigma[parm[parm%in%rownames(cf_sigma)], ,drop=FALSE]
    cf_zero <- cf_zero[parm[parm%in%rownames(cf_zero)], ,drop=FALSE] 
    cf_one <- cf_one[parm[parm%in%rownames(cf_one)], ,drop=FALSE] 
    
    cf <- rbind(cf_mean, cf_sigma, cf_zero, cf_one)
    }

  cf_est <- cf[, 1]
  cf_ses <- cf[, 2]
  
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  fac <- qnorm(a)
  ci_name <- paste0(round(a*100, 2),"%")
  ci <- array(NA, dim = c(length(parm), 2L), 
              dimnames = list(parm, ci_name))
  ci[] <- cf_est + cf_ses %o% fac
  
  if(type == "full"){
    location <- ci[rownames(cf_mean), ,drop=FALSE]
    dispersion <- ci[rownames(cf_sigma), ,drop=FALSE]
    zero <- ci[rownames(cf_zero), ,drop=FALSE]
    one <- ci[rownames(cf_one), ,drop=FALSE]
    
    ci <- list(location = location,
               dispersion = dispersion,
               zero = zero,
               one = one)
  }
  
  ci
  
  
} 
