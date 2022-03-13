#' @title Plot Fitted Values/Residuals of A Cdfqr Object or Distribution
#' @aliases plot.cdfqr
#' @description Plot Fitted Values/Residuals of A cdfqr Object or Distribution
#' @param x If the plot is based on the fitted values, provide a fitted cdfqr object, alternatively, mu and sigma, and the distribution can be specified.
#' @param n The number of random variates to be generated for user specified plot. 
#' @param type Currently only fitted values are available for generating plots.
#' @param fd A string that specifies the parent distribution.
#' @param sd A string that specifies the sub-family distribution.
#' @param mu Location parameter value
#' @param sigma Sigma parameter value
#' @param theta Skew parameter value
#' @param inner If finite-tailed distribution is used: a logic value that indicates if the inner (\code{inner = TRUE}) case or outer (\code{inner = FALSE}) will be used. Currently inner case can only be used for 2-parameter distributions.
#' @param version If finite-tailed distribution is used: A string indicates that which version will be used. "V" is the tilt parameter function while "W" indicates the Jones Pewsey transformation. 
#' @param ... other plot parameters pass onto \code{\link[graphics]{plot}}.
#'
#' @examples
#' data(cdfqrExampleData)
#' fit <- cdfquantreg(crc99 ~ vert | confl, 't2','t2', data = JurorData)
#' plot(fit)
#' 
#' 
#' @method plot cdfqr
#' @export
#' @import graphics
#' @importFrom MASS truehist

plot.cdfqr <- function(x, mu = NULL, sigma = NULL, theta = NULL,
                       fd = NULL, sd = NULL, n = 10000, 
                       inner = TRUE, version = "V",
                       type = c("fitted"), ...) {
  if("cdfqrH"%in%class(x)){
    stop("Plot is currently not available for hurdle models")
  }
  # If plot based on the fitted model
  if(is.null(sigma)){
    ydata <- x$y
    fd <- x$family$fd
    sd <- x$family$sd
    
    if(!is.null(x$fitted$mu)){
      mu <- mean(x$fitted$mu)
    }
 
    sigma <- exp(mean(log(x$fitted$sigma)))
    
    
    if("cdfqrFT"%in%class(x)){
      theta <- mean(x$fitted$theta)
      fitted <- x$fitted$full
      fit_d <- pdfft(fitted, sigma, theta, fd, sd, mu = mu, 
                     inner = x$distinf$inner, 
                     version = x$distinf$version)
    }else{
      fitted <- scaleTR(x$fitted$full) # smooth the boundary values
      fit_d <- dq(fitted, mu, sigma, fd, sd)
    }
  
    xtem <- data.frame(x = as.numeric(fitted), y = as.numeric(fit_d))
    xtem <- xtem[order(xtem$x),]
    
    par(mfrow=c(1,2))
    MASS::truehist(ydata, col = "white", ymax = max(xtem$y) + 0.1,
                   main = "Data (histogram) \n fitted by model (line)",
                   xlab = 'Observations', ylab = 'Density',
                   ...)
    graphics::lines(xtem$x, xtem$y, lty = 1, lwd = 2)

    
    # plot data against fitted values
    plot(ydata,x$fitted$full, xlab = 'Observations', ylab = 'Fitted',
         main = "Observations vs. Fitted",type ='p', pch=20)
    graphics::abline(0, 1)
    
  }else{
    #based on input mu and sigma, and given distribution
    if("cdfqrFT"%in%class(x)){
      ydata <- as.numeric(rqft(n, sigma, theta, fd, sd, mu = mu, 
                                     inner = inner, 
                                     version = version))
      fit_d <- as.numeric(pdfft(x, sigma, theta, fd, sd, mu = mu, 
                                inner = x$distinf$inner, 
                                version = x$distinf$version))
    } else{
      ydata <- as.numeric(rq(n, mu, sigma, fd, sd))
      fit_d <- as.numeric(dq(ydata, mu, sigma, fd, sd))
    }
    
   xtem <- data.frame(x = ydata, y =fit_d)
   xtem <- xtem[order(xtem$x),]
   xlim <- c(0,1)
   MASS::truehist(ydata, col = "white", ymax = max(xtem$y) + 0.1, ...)
   graphics::lines(xtem$x, xtem$y, lty = 1, lwd = 2)
  }
  
  invisible()
}

