#' @title Register method for cdfqr object functions
#' @description Register method for cdfqr object functions.
#' @aliases residuals.cdfqr
#' @param object The cdfqr model project  
#' @param type The type of residuals to be extracted: \code{'raw'}, \code{'pearson'},\code{'std.pearson'}, or \code{'deviance'},
#' @param ... currently ignored 
#' @return residuals of a specified type.
#' @method residuals cdfqr
#' 
#' @examples
#' data(cdfqrExampleData)
#' fit <- cdfquantreg(crc99 ~ vert | confl, 't2','t2', data = JurorData)
#' 
#' residuals(fit, "pearson")

#' @export
residuals.cdfqr <- function(object, type = c("raw","pearson", "deviance"), ...) {
  
  type <- match.arg(type)
  
  residuals <- object$residuals  # Obtain the raw residuals
  pearson <- object$residuals/as.numeric(sqrt(var(object$fitted$full)))

  type <- match.arg(type)
  
  ydata <- object$y  # observed data
  n <- length(ydata)  # number of observations
  fitted <- fitted(object,"full")  # model fitted values

  sigma <- fitted(object, type = "sigma")  # fitted sigma values values

  dist <- object$family
  fd <- dist$fd
  sd <- dist$sd
  if(!is.null(object$fitted$mu))
  {
    mu <- fitted(object, type = "mu")  # fitted mu values
  }else{
    mu <- NULL
  }
  
  if("cdfqrFT" %in% class(object)){
    theta <- fitted(object, type = "theta")  
    deviance_r <- sign(residuals) * sqrt(2 * abs(
      -log(pdfft( ydata, sigma = sigma, theta = theta, fd = fd,sd = sd, mu = mu, inner = object$distinf$inner, version = object$distinf$version)) + 
        log(pdfft(fitted, sigma = sigma, theta = theta,fd = fd,sd = sd, mu = mu, inner = object$distinf$inner, version = object$distinf$version))))
  }else{
    deviance_r <- sign(residuals) * sqrt(2 * abs(-qrLogLik(ydata, mu, sigma, fd, sd, total = FALSE) + qrLogLik(fitted, 
                                                                                                               mu, sigma, fd, sd, total = FALSE)))
  }

  
  if("cdfqrH" %in% class(object)){
    if(!is.na(object$mod_zp)){
      pearson_zero <- residuals(object$mod_zp, type = "pearson")
      deviance_zero <- residuals(object$mod_zp, type = "deviance")
    }else{
      pearson_zero = NULL
      deviance_zero = NULL
    }
    if(!is.na(object$mod_op)){
      pearson_one <- residuals(object$mod_op, type = "pearson")
      deviance_one <- residuals(object$mod_op, type = "deviance")
    }else{
      pearson_one = NULL
      deviance_one = NULL
    }
    pearson = cbind(cdfqr =pearson,
                   zero_part = pearson_zero,
                   one_part = pearson_one
                     )
    deviance_r = cbind(cdfqr = deviance_r,
                    zero_part = deviance_zero,
                    one_part = deviance_one
    )
  }


res <- switch(type, raw = {
    residuals
  }, pearson = {
    pearson
  },  deviance = {
    deviance_r
  })
  
  return(res)
} 
