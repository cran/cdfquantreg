#' @title Methods for Cdfqr Objects
#' @aliases predict.cdfqr fitted.cdfqr
#' @description Methods for obtaining the fitted/predicted values for a fitted cdfqr object.
#' @param object A cdfqr model fit object
#' @param type A character that indicates whether the full model prediction/fitted values are needed, or values for the `mu` and `sigma` submodel only. 
#' @param newdata Optional. A data frame in which to look for variables with which to predict. If not provided, the fitted values are returned 
#' @param quant A number or a numeric vector (must be in (0, 1)) to specify the quantile(s) of the predicted value (when `newdata` is provided, and predicted values for responses are required). The default is to use median to predict response values.   
#' @param plot if a plot is needed. 
#' @param ... currently ignored
#' @examples
#' data(cdfqrExampleData)
#' fit <- cdfquantreg(crc99 ~ vert | confl, 't2','t2', data = JurorData)
#' 
#' plot(predict(fit))
#' plot(predict(fit))
#' @method predict cdfqr 
#' @export 
#' 
predict.cdfqr <- function(object, newdata = NULL, type = c("full","mu","sigma", "theta", "one", "zero"), quant = 0.5,...) {
  type <- match.arg(type)
  
 if (!is.null(newdata)){
   tt <- terms(object$formula)
   Terms <- delete.response(tt)
   fd <- object$family$fd
   sd <- object$family$sd
   m <- model.frame(Terms, newdata)
   X <- model.matrix(Terms, m)
   

   if(!is.null(object$coefficients$location)){
     betas.lm <- object$coefficients$location[,"Estimate"]
     term.lm <- rownames(object$coefficients$location)
     
     if(length(term.lm) < 2){
       predictlm <- X[, term.lm] * betas.lm
     }else{
       predictlm <- X[, term.lm] %*% betas.lm
     }
   }else{
     predictlm = NULL
   }

   
   betas.pm <- object$coefficients$dispersion[,"Estimate"]
   term.pm <- rownames(object$coefficients$dispersion)
   
   if(length(term.pm) < 2){
     predictpm <- X[, term.pm] * betas.pm
   }else{
     predictpm <- X[, term.pm] %*% betas.pm
   }
  
   predictsigma <- exp(predictpm)
   
   
   if("cdfqrFT"%in%class(object)){
     
     betas.tm <- object$coefficients$skew[,"Estimate"]
     term.tm <- rownames(object$coefficients$skew)
     
     if(length(term.tm) < 2){
       predicttm <- X[, term.tm] * betas.tm
     }else{
       predicttm <- X[, term.tm] %*% betas.tm
     }
     
     predicttheta <- predicttm
     
   }
   
   n <- length(predictsigma)
   if (length(quant) == 1) {
     q_y <- rep(quant, n)
   }else if (length(quant) == n){
     q_y <- quant
   }
   
   if("cdfqrFT"%in%class(object)){
     modelspecific <- object$distinf
     pred <- qqft(q_y, predictsigma, predicttheta, fd, sd, 
                  mu = predictlm, inner = object$distinf$inner, 
                  version = object$distinf$version)
     
   }else{
     pred <- qq(q_y, predictlm, predictsigma, fd, sd)
   }

   
   
   if("cdfqrH"%in%class(object)){
     if (!is.null(object$mod_zp)){
       pred_zero <- predict(object$mod_zp, newdata = newdata,
                            type = "response")
     }
     if (!is.null(object$mod_op)){
       pred_one <- predict(object$mod_op, newdata = newdata,
                           type = "response")
     }
     
   }

   
   
   
   fitted <- switch(type, full = {
     pred
   }, mu = {
     predictlm
   }, sigma = {
     predictsigma
   }, theta = {
     predicttheta
   }, zero = {
     pred_zero
   }, one = {
     pred_one
   })
 }else{
   fitted <- switch(type, full = {
     object$fitted$full
   }, mu = {
     object$fitted$mu
   }, sigma = {
     object$fitted$sigma
   }, theta = {
     object$fitted$theta
   }, zero = {
     object$fitted$zerop
   }, one = {
     object$fitted$onep
   })
 }
  
  return(fitted)
}

#' @export 
#' @rdname predict.cdfqr
fitted.cdfqr <- function(object, type = c("full","mu","sigma", "theta", "one", "zero"), plot = FALSE, ...) {
  type <- match.arg(type)
  
  fitted <- switch(type, full = {
    object$fitted$full
  }, mu = {
    object$fitted$mu
  }, sigma = {
    object$fitted$sigma
  }, theta = {
    object$fitted$theta
  }, zero = {
    object$fitted$zerop
  }, one = {
    object$fitted$onep
  })
  
  if (plot) {
    if("cdfqrH"%in%class(object)){
      cat("Plot is not available for hurdle model and finite tailed model currently")
    }else{
      plot1 <- plot(object)
      return(list(fitted = fitted,
                  plot = plot1))
    }

  }else{
    return(fitted)
    }

}
