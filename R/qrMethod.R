#' @title S3 Methods for getting output from fitted cdfqr Objects.
#' 
#' @aliases summary.cdfqr
#' 
#' @description Give the S3 Methods for CDF-Quantile Distribution Models
#' 
#' @param x,object The fitted cdfqr model.
#' @param type,submodel The parts of coefficients or variance-covariance matrix to be extracted.Can be "full", "mean",or "sigma". 
#' @param formula. Changes to the formula. See \code{\link[Formula]{update.Formula}} for details.
#' @param zero.fo.,one.fo.,mu.fo., Changes to the formulas for zero/one component for hurdle models, and for location submodel for finite-tailed models.
#' @param evaluate If true evaluate the new updated model else return the call for the new model.
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level the confidence level required.
#' @param digits Number of digits to be retained in printed output.
#' @param ... Pass onto other functions or currently ignored
#' @examples
#' data(cdfqrExampleData)
#' fit <- cdfquantreg(crc99 ~ vert | confl, 't2','t2', data = JurorData)
#' 
#' summary(fit)
#' print(fit)
#' logLik(fit)
#' coef(fit)
#' deviance(fit)
#' vcov(fit)
#' confint(fit)
#' 
#' #Update the model
#' fit2 <- update(fit, crc99 ~ vert*confl | confl)
#' summary(fit2)

#' @export
#' @method summary cdfqr
summary.cdfqr <- function(object, ...) {
  print(object)
}

#' @export
summary.cdfqrH <-function(object, ...) {
  print(object)
}

#' @export
summary.cdfqrFT <- function(object, ...) {
  print(object)
}

#' @export
#' @rdname summary.cdfqr
print.cdfqr <- 
  function(x, digits = max(3, getOption("digits") - 3), ...){
    .print.cdfqr(x, digits = digits)
  }

#' @export
print.cdfqrFT <- 
  function(x, digits = max(3, getOption("digits") - 3), ...){
    .print.cdfqr(x, digits = digits)
  }

#' @export
print.cdfqrH <-
  function(x, digits = max(3, getOption("digits") - 3), ...){
    .print.cdfqr(x, digits = digits)
  }
.print.cdfqr <-
  function(x, digits = max(3, getOption("digits") - 3), ...) {
    object <- x
    cat("Family: ",
        paste(object$family$fd, object$family$sd, collapse = "-"),
        "\n")
    cat("Call: ", deparse(object$call), "\n", fill = TRUE)
    
    if (length(object$coefficients$location)) {
      mutitile <- "Mu coefficients (Location submodel)\n"
      if (object$family$fd == 'km') {
        mutitile <- "1st shape parameter coefficients (log link)\n"
      }
      cat(mutitile)
      printCoefmat(object$coefficients$location,
                   digits = digits,
                   signif.legend = FALSE)
      cat("\n")
    } else{
      mutitile2 <- "No coefficients in location submodel\n\n"
      if (object$family$fd == 'km') {
        mutitile2 <-
          "No coefficient for the 1st shape parameter (log link)\n"
      }
      cat(mutitile2)
    }

    if (length(object$coefficients$dispersion)) {
      if (length(object$coefficients$dispersion)) {
        Sigmatitile <- "Sigma coefficients (Dispersion submodel)\n"
        if (object$family$fd == 'km') {
          Sigmatitile <- "2nd shape parameter coefficients (log link)\n"
        }
        cat(Sigmatitile)
        printCoefmat(object$coefficients$dispersion, digits = digits)
        cat("\n")
      } else{
        Sigmatitile2 <- "No coefficients in Dispersion submodel\n"
        if (object$family$fd == 'km') {
          Sigmatitile2 <-
            "No coefficient for the 2nd shape parameter (log link)\n"
        }
        cat(Sigmatitile2)
      }
    }
    
    if ("cdfqrFT" %in% class(object)) {
      if (length(object$coefficients$skew)) {
        cat("Theta coefficients (Skewness submodel)\n")
        printCoefmat(object$coefficients$skew, digits = digits)
        cat("\n")
      }
      else{
        cat("No coefficients in the skewness submodel\n")
      }
      
    }
    
    if ("cdfqrH" %in% class(object)) {
      if (length(object$coefficients$zero)) {
        mutitile <- "Zero component coefficients\n"
        cat(mutitile)
        printCoefmat(object$coefficients$zero,
                     digits = digits,
                     signif.legend = FALSE)
        cat("\n")
      }
      if (length(object$coefficients$one)) {
        mutitile <- "One component coefficients\n"
        cat(mutitile)
        printCoefmat(object$coefficients$one,
                     digits = digits,
                     signif.legend = FALSE)
        cat("\n")
      }
    }
    
    
    
    
    coverge_msg <- ifelse(object$converged == 0,
                          "successful completion",
                          "has reached the maxit iterations")
    
    cat("Converge: ", coverge_msg, fill = TRUE)
    cat("Log-Likelihood: ",
        round(object$logLik, digits = digits),
        "\n",
        fill = TRUE)
    if (!is.null(object$grad))
      cat("Gradient: ", round(object$grad, digits = digits), "\n", fill = TRUE)
    invisible(object)
    
  }

#' @export
#' @rdname summary.cdfqr
coef.cdfqr <- function(object, type = "full", ...){
  .coef.cdfqr(object, type = type)
}

#' @export
coef.cdfqrFT <-function(object, type = "full", ...){
  .coef.cdfqr(object, type = type)
}

#' @export

coef.cdfqrH <- function(object, type = "full", ...){
  .coef.cdfqr(object, type = type)
}

.coef.cdfqr <-
  function(object,
           type = c("full", "mean", "sigma", "skew", "zero", "one")) {
    type <- match.arg(type)
    
    sigma <- object$coefficients$dispersion[, 1]
    names(sigma) <- paste('(sigma)_', names(sigma), sep = "")
    
    if ("cdfqrFT" %in% class(object)) {
      skew <- object$coefficients$skew[, 1]
      names(skew) <- paste('(theta)_', names(skew), sep = "")
      if (!is.null(object$coefficients$location)) {
        mu <- object$coefficients$location[, 1]
        names(mu) <- paste('(mu)_', names(mu), sep = "")
      } else{
        mu = NULL
      }
      full <-  c(skew, sigma, mu)
    } else{
      skew = NULL
      mu <- object$coefficients$location[, 1]
      names(mu) <- paste('(mu)_', names(mu), sep = "")
      full <-  c(mu, sigma)
    }
    
    
    if ("cdfqrH" %in% class(object)) {
      if (!is.null(object$coefficients$zero)) {
        zero <- object$coefficients$zero[, 1]
        names(zero) <- paste('(zero)_', names(zero), sep = "")
      } else{
        zero <- NULL
      }
      
      
      if (!is.null(object$coefficients$one)) {
        one <- object$coefficients$one[, 1]
        names(one) <- paste('(one)_', names(one), sep = "")
      } else{
        one <- NULL
      }
      full <-  c(full, zero, one)
    } else{
      if (type == "zero" | type == "one") {
        zero = one = "The model is not a hurdle model. No zero/one component is involved"
      }
    }
    
    
    
    
    
    coef <- switch(
      type,
      full = {
        full
      },
      mean = {
        mean
      },
      skew = {
        skew
      },
      sigma = {
        sigma
      },
      zero = {
        zero
      },
      one = {
        one
      }
    )
    return(coef)
  }


#' @export
#' @rdname summary.cdfqr
vcov.cdfqr <- vcov.cdfqrFT <- vcov.cdfqrH <- function(object,
                                                      type = "full", ...){
  .vcov.cdfqr(object, type = type)
}


#' @export
vcov.cdfqrFT <-function(object,
                                                      type = "full", ...){
  .vcov.cdfqr(object, type = type)
}

#' @export
vcov.cdfqrH <- function(object,
                                                      type = "full", ...){
  .vcov.cdfqr(object, type = type)
}

.vcov.cdfqr <- function(object, type = c("full","mean","sigma","skew", "theta", "zero", "one"), ...) {
  

  type = match.arg(type)
  
  if( "cdfqrH" %in% class(object)){
    vc <- object$vcov$cdf
    zero <- object$vcov$zero
    one <- object$vcov$one
  }else{
    vc <- object$vcov
    if(type == "zero"|type == "one"){
      zero = one = "The model is not a hurdle model. No zero/one component is involved"
    }
  }

  
  if( "cdfqrFT" %in% class(object)){
    thetaind <- grep("theta", colnames(vc))
    skenwess <- vc[thetaind, thetaind, drop = FALSE]
    
  }
  
  if(!is.null(object$coefficients$location)){
    muind <- grep("mu", colnames(vc))
    location <- vc[muind, muind, drop = FALSE]
    
  }else{
    location = NULL
  }
    
  sigmaind <- grep("sigma", colnames(vc))
  dispersion <- vc[sigmaind, sigmaind, drop = FALSE]
  
  vc <- switch(type, full = {vc}, 
               mean = {location}, 
               sigma = {dispersion},
               theta = {skenwess}, 
               zero = {zero}, 
               one = {one})
  return(vc)
}

#' @export
#' @rdname summary.cdfqr
update.cdfqr <- function(object, formula.,zero.fo., one.fo., mu.fo., ..., evaluate = TRUE){
  .update.cdfqr(object, formula.,zero.fo., one.fo., mu.fo., ..., evaluate = evaluate)
}

#' @export
update.cdfqrH <- function(object, formula.,zero.fo., one.fo., mu.fo., ..., evaluate = TRUE){
  .update.cdfqr(object, formula.,zero.fo., one.fo., mu.fo., ..., evaluate = evaluate)
}

#' @export
update.cdfqrFT <- function(object, formula.,zero.fo., one.fo., mu.fo., ..., evaluate = TRUE){
  .update.cdfqr(object, formula.,zero.fo., one.fo., mu.fo., ..., evaluate = evaluate)
}

.update.cdfqr <- function(object, formula.,zero.fo., one.fo., mu.fo., ..., evaluate = TRUE) {
  call <- object$call
  if (!missing(formula.))
    call$formula <- formula(update(Formula::Formula(formula(object)), formula.))

  if("cdfqrH" %in% class(object)){
    if (!missing(zero.fo.)) 
      call$zero.fo <- formula(update(Formula::Formula(formula(call$zero.fo)), zero.fo.))
    
    if (!missing(one.fo.)) 
      call$one.fo <- formula(update(Formula::Formula(formula(call$one.fo)), one.fo.))
  }else if("cdfqrFT" %in% class(object)){
    if (!missing(mu.fo.)) 
      call$mu.fo <- formula(update(Formula::Formula(formula(call$zero.fo)), mu.fo.))
  }else{
    if (!missing(mu.fo.)|!missing(zero.fo.)|!missing(one.fo.))
      cat("The supply of zero/one/mu function is ignored as the model is not hurdle/finite tailed model.")
  }
  
  
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

#' @export
#' @rdname summary.cdfqr
confint.cdfqr <- function(object, parm, level = 0.95, submodel = "full", ...){
  .confint.cdfqr(object, parm, level = level, submodel = submodel)
}

.confint.cdfqr <- function(object, parm, level = 0.95, submodel = "full") {
 
  cf_full <- object$coefficients
  cf_mean <- cf_full['location'][[1]]
  cf_sigma <- cf_full['dispersion'][[1]]
  cf <- do.call(rbind, object$coefficients)
  
    
  if (submodel == "location") {
    if(missing(parm)) parm.l = rownames(cf_mean) else parm.l = parm
    
      if (is.character(parm.l)){
        if (!all(parm.l %in% rownames(cf_mean)))
        stop("One or more required parameter(s) is not in the location submodel!")}
      
      if (is.numeric(parm.l)){
        if (!is.matrix(try(cf_mean[parm.l, ],T)))
          stop("One or more required parameter(s) is not in the location submodel!")}
      
    cf <- cf_mean<- cf_mean[parm.l, ,drop=FALSE]
    parm <- parm.l
    }
    
  if (submodel == "sigma") {
    if(missing(parm)) parm.d = rownames(cf_sigma) else parm.d = parm
    
    if (is.character(parm.d)){
      if (!all(parm.d %in% rownames(cf_sigma)))
        stop("One or more required parameter(s) is not in the dispersion submodel!")}
      
      if (is.numeric(parm.d)){
        if (!is.matrix(try(cf_sigma[parm.d, ],T)))
          stop("One or more required parameter(s) is not in the dispersion submodel!")}
      
     cf <- cf_sigma<- cf_sigma[parm.d, ,drop=FALSE]
     parm <- parm.d
    }
    
  if (submodel == "full") {
    if(missing(parm)) {
      parm.l = rownames(cf_mean)
      parm.d = rownames(cf_sigma)}else{
        parm.l = parm[[1]]
        parm.d = parm[[2]]
      }
    
    if (is.character(parm.l)){
      if (!all(parm.l %in% rownames(cf_mean)))
        stop("One or more required parameter(s) is not in the location submodel!")}
    if (is.numeric(parm.l)){
      if (!is.matrix(try(cf_mean[parm.l, ],T)))
        stop("One or more required parameter(s) is not in the location submodel!")}
    if (is.character(parm.d)){
      if (!all(parm.d %in% rownames(cf_sigma)))
        stop("One or more required parameter(s) is not in the dispersion submodel!")}
    
    if (is.numeric(parm.d)){
      if (!is.matrix(try(cf_sigma[parm.d, ],T)))
        stop("One or more required parameter(s) is not in the dispersion submodel!")}
    
      
    cf_mean <- cf_mean[rownames(cf_mean)[parm.l%in%rownames(cf_mean)], ,drop=FALSE]
    cf_sigma <- cf_sigma[rownames(cf_sigma)[parm.d%in%rownames(cf_sigma)], ,drop=FALSE]
    cf <- rbind(cf_mean, cf_sigma)
    parm <- c(parm.l, parm.d)
    }

  cf_est <- cf[, 1]
  cf_ses <- cf[, 2]
  
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  fac <- qnorm(a)
  ci_name <- paste0(round(a*100, 2),"%")
  ci <- array(NA, dim = c(nrow(cf), 2L), 
              dimnames = list(parm, ci_name))
  ci[] <- cf_est + cf_ses %o% fac
  
  if(submodel == "full"){
    location <- ci[1:nrow(cf_mean), ,drop=FALSE]
    dispersion <- ci[(1+nrow(cf_mean)):nrow(ci), ,drop=FALSE]
    ci <- list(location = location,
               dispersion = dispersion)
  }
  
  ci
  
  
} 

#' @export
#' @rdname summary.cdfqr
 formula.cdfqr <- function(x, ...) {
  return(formula(x$call$formula))
} 

#' @export
formula.cdfqrFT <-  function(x, ...) {
  return(formula(x$call$formula))
} 

#' @export
#' @rdname summary.cdfqr
nobs.cdfqr <- function(object, ...) {
  length(object$residuals)
}

#' @export
 nobs.cdfqrH <- function(object, ...) {
  length(object$residuals)
}

#' @export
nobs.cdfqrFT <- function(object, ...) {
  length(object$residuals)
}

#' @export
#' @rdname summary.cdfqr
deviance.cdfqr <- function(object, ...) {
  object$deviance
}
#' @export
 deviance.cdfqrFT <- function(object, ...) {
  object$deviance
}
#' @export
deviance.cdfqrH <- function(object, ...) {
  object$deviance
}


#' @export
#' @rdname summary.cdfqr
logLik.cdfqrH <-  function(object, ...) {
  object$logLik
}

#' @export
logLik.cdfqrFT <- function(object, ...) {
  object$logLik
}

#' @export
logLik.cdfqr <- function(object, ...) {
  object$logLik
}

#' @export
#' @rdname summary.cdfqr
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

#' @export
#' @rdname summary.cdfqr
confint.cdfqrFT <- function(object, parm, level = 0.95, submodel = "full", ...) {

  cf_full <- object$coefficients
  cf_theta <- cf_full['skew'][[1]]
  cf_sigma <- cf_full['dispersion'][[1]]
  if(!is.null(object$coefficients$location)){cf_mu <- cf_full['location'][[1]]}else{
    cf_mu = NULL
  }
  cf <- do.call(rbind, object$coefficients)
  
  
  if (submodel == "theta") {
    if(missing(parm)) parm.t = rownames(cf_theta) else parm.l = parm
    
    if (is.character(parm.t)){
      if (!all(parm.t %in% rownames(cf_theta)))
        stop("One or more required parameter(s) is not in the skewness submodel!")}
    
    if (is.numeric(parm.t)){
      if (!is.matrix(try(cf_theta[parm.t, ],T)))
        stop("One or more required parameter(s) is not in the skewness submodel!")}
    
    cf <- cf_theta<- cf_theta[parm.t, ,drop=FALSE]
    parm <- parm.t
  }
  
  if (submodel == "sigma") {
    if(missing(parm)) parm.d = rownames(cf_sigma) else parm.d = parm
    
    if (is.character(parm.d)){
      if (!all(parm.d %in% rownames(cf_sigma)))
        stop("One or more required parameter(s) is not in the dispersion submodel!")}
    
    if (is.numeric(parm.d)){
      if (!is.matrix(try(cf_sigma[parm.d, ],T)))
        stop("One or more required parameter(s) is not in the dispersion submodel!")}
    
    cf <- cf_sigma<- cf_sigma[parm.d, ,drop=FALSE]
    parm <- parm.d
  }
  
  if (submodel == "sigma"&!is.null(cf_mu)) {
    if(missing(parm)) parm.l = rownames(cf_mu) else parm.l = parm
    
    if (is.character(parm.l)){
      if (!all(parm.d %in% rownames(cf_mu)))
        stop("One or more required parameter(s) is not in the dispersion submodel!")}
    
    if (is.numeric(parm.l)){
      if (!is.matrix(try(cf_mu[parm.l, ],T)))
        stop("One or more required parameter(s) is not in the dispersion submodel!")}
    
    cf <- cf_mu<- cf_mu[parm.l, ,drop=FALSE]
    parm <- parm.l
  }else if (submodel == "sigma"&is.null(object$coefficients$mu)) {
    stop("The model does not have the location component!")
  }
  
  if (submodel == "full") {
    if(missing(parm)) {
      parm.t = rownames(cf_theta)
      parm.d = rownames(cf_sigma)
      if(!is.null(cf_mu)){
        parm.l <-rownames(cf_mu)
      }else{
        parm.l <- NULL
      }
    }else{
      parm.t = parm[[1]]
      parm.d = parm[[2]]
      if(is.null(cf_mu)) parm.l = NULL else parm.l = parm[[3]] 
    }
    
    if (is.character(parm.t)){
      if (!all(parm.t %in% rownames(cf_theta)))
        stop("One or more required parameter(s) is not in the skewness submodel!")}
    if (is.numeric(parm.t)){
      if (!is.matrix(try(cf_theta[parm.t, ],T)))
        stop("One or more required parameter(s) is not in the skewness submodel!")}
    if (is.character(parm.d)){
      if (!all(parm.d %in% rownames(cf_sigma)))
        stop("One or more required parameter(s) is not in the dispersion submodel!")}
    
    if (is.numeric(parm.d)){
      if (!is.matrix(try(cf_sigma[parm.d, ],T)))
        stop("One or more required parameter(s) is not in the dispersion submodel!")}
    
    
    
    cf_theta <- cf_theta[rownames(cf_theta)[parm.t%in%rownames(cf_theta)], ,drop=FALSE]
    cf_sigma <- cf_sigma[rownames(cf_sigma)[parm.d%in%rownames(cf_sigma)], ,drop=FALSE]
    if(!is.null(cf_mu)){
      cf_mu <- cf_mu[rownames(cf_mu)[parm.l%in%rownames(cf_mu)], ,drop=FALSE]
    }
    
    cf <- rbind(cf_sigma, cf_theta, cf_mu)
    parm <- c(parm.d, parm.t, parm.l)
  }
  
  cf_est <- cf[, 1]
  cf_ses <- cf[, 2]
  
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  fac <- qnorm(a)
  ci_name <- paste0(round(a*100, 2),"%")
  ci <- array(NA, dim = c(nrow(cf), 2L), 
              dimnames = list(parm, ci_name))
  ci[] <- cf_est + cf_ses %o% fac
  
  if(submodel == "full"){
    sigma<- ci[1:nrow(cf_sigma), ,drop=FALSE]
    theta <- ci[(1+nrow(cf_sigma)):(nrow(cf_sigma) +nrow(cf_theta)), ,drop=FALSE]
    if(!is.null(cf_mu)){
      mu <- ci[(1+ nrow(cf_sigma) +nrow(cf_theta)):nrow(ci), ,drop=FALSE]
    }else{
      mu = NULL
    }
    
    ci <- list(
      sigma = sigma,
      theta = theta,
      mean = mu)
  }
  
  ci
  
  
} 




