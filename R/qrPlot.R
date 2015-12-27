#' @title Plot Fitted Values/Residuals of A Cdfqr Object or Distribution
#' @aliases plot.cdfqr
#' @description Plot Fitted Values/Residuals of A Cdfqr Object or Distribution
#' @param x If the plot is based on the fitted values, provide a fitted cdfqr object.
#'
#' @param mu,sigma, fd, sd alternatively, mu and sigma, and the distribution can be specified
#' @param n The number of random variates to be generated for user specified plot. 
#' @param type Currently only fitted values are available for generating plots.
#' @param fd A string that specifies the parent distribution.
#' @param sd A string that specifies the sub-family distribution.
#' @param ... other plot parameters pass onto \code{\link[graphics]{plot}}.
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

plot.cdfqr <- function(x, mu = NULL, sigma = NULL, fd = NULL, sd = NULL, n = 10000, 
                       type = c("fitted"),...) {
  
  # If plot based on the fitted model
  if(is.null(mu)){
    ydata <- x$y
    fit_d <- density(x$fitted$full)
 
  }else{
    #based on input mu and sigma, and given distribution
   ydata <- rq(n, mu, sigma, fd, sd)
   fit_d <- density(rq(n, mu, sigma, fd, sd))
  }
  
  MASS::truehist(ydata, col = "white", ...)
  graphics::lines(fit_d, lty = 1, lwd = 2)
  invisible()
}

