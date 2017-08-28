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
predict.cdfqr <- function(object, newdata = NULL, type = c("full","mu","sigma"), quant = 0.5,...) {
  type <- match.arg(type)
  
 if (!is.null(newdata)){
   tt <- terms(object$formula)
   Terms <- delete.response(tt)
   fd <- object$family$fd
   sd <- object$family$sd
   m <- model.frame(Terms, newdata)
   X <- model.matrix(Terms, m)
   
   
   betas.lm <- object$coefficients$location[,"Estimate"]
   betas.pm <- object$coefficients$dispersion[,"Estimate"]
   term.lm <- rownames(object$coefficients$location)
   term.pm <- rownames(object$coefficients$dispersion)
   
   if(length(term.lm) < 2){
     predictlm <- X[, term.lm] * betas.lm
   }else{
     predictlm <- X[, term.lm] %*% betas.lm
   }
   
   if(length(term.pm) < 2){
     predictpm <- X[, term.pm] * betas.pm
   }else{
     predictpm <- X[, term.pm] %*% betas.pm
   }
  
   predictsigma <- exp(predictpm)
   
   n <- length(predictsigma)
   if (length(quant) == 1) {
     q_y <- rep(quant, n)
   }else if (length(quant) == n){
     q_y <- quant
   }
  
   pred <- qq(q_y, predictlm, predictsigma, fd, sd)
   fitted <- switch(type, full = {
     pred
   }, mu = {
     predictlm
   }, sigma = {
     predictsigma
   })
 }else{
   fitted <- switch(type, full = {
     object$fitted$full
   }, mu = {
     object$fitted$full_mu
   }, sigma = {
     object$fitted$full_sigma
   })
 }
  
  return(fitted)
}

#' @method fitted cdfqr 
#' @export 
#' @rdname predict.cdfqr
fitted.cdfqr <- function(object, type = c("full","mu","sigma"), plot = FALSE, ...) {
  type <- match.arg(type)
  
  fitted <- switch(type, full = {
    object$fitted$full
  }, mu = {
    object$fitted$mu
  }, sigma = {
    object$fitted$sigma
  })
  if (plot) {
    plot1 <- plot(object)
    return(list(fitted = fitted,
           plot = plot1))
  }else{
    return(fitted)
    }

}