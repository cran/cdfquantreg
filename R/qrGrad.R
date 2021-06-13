#' @title Give the Gradient Function for CDF-Quantile Distribution Modles
#' @description Give the Gradient Function for CDF-Quantile Distribution Modles. 
#' @aliases qrGrad
#' @param fd A string that specifies the parent distribution.
#' @param sd A string that specifies the sub-family distribution.
#' 
#' @return grad The gradient function of parameter estimates, given a specified cdf-quantile distribution
#' 
#' @export
#' 
#' @examples
#' qrGrad('t2','t2')
qrGrad <- function(fd, sd) {
  
  # arcsinh-XX-----
  if (fd == "arcsinh") {
    
    # arcsinh-arcsinh
    if (sd == "arcsinh") {
       a_arcsinh_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <-
          cbind(x * rep(
            -(1 / (exp(gz) * sqrt(
              1 + (1 / (-1 + y) + 1 / (2 * y - 2 * y ^ 2) + hx) ^ 2 / exp(2 * gz)
            ))) + (2 * exp(asinh(((1 - 2 * y) / (2 * (-1 + y) * y) - hx
            ) / exp(gz)) - gz)) / ((1 +
                            exp(asinh(((1 - 2 * y) /(2 * (-1 + y) * y) - hx) / exp(gz)
                            ))) * sqrt(1 + (1 / (-1 + y) + 1 / (2 * y - 2 * y ^ 2) + hx
                            ) ^ 2 / exp(2 * gz))) - 
              (2 * (-1 + y) * y * (-1 + 2 * y ^ 2 * hx + y * (2 - 2 * hx))) /
              (1 + 4 * y * (-1 + hx) + 4 * y ^ 4 * (exp(2 * gz) + hx ^ 2) + 4 * y ^ 2 *
                  (1 + exp(2 * gz) - 3 * hx + hx ^ 2) - 8 * y ^ 3 * (exp(2 * gz) + hx * (-1 + hx))),
            length(x[1,])), z * rep(-(((( 1 - 2 * y ) / (2 * (-1 + y) * y) - hx)) / (exp(gz) * sqrt(
              1 + (1 / (-1 + y) + 1 / (2 * y - 2 * y ^ 2) + hx) ^ 2 / exp(2 * gz)
            ))) + (2 * exp(asinh(((1 - 2 * y) / (2 * (-1 + y) * y) - hx) / exp( gz
            )) - gz) * ((1 - 2 * y) / (2 * (-1 + y) * y) - hx)) / 
              ((1 + exp(asinh(((1 - 2 * y) / (2 * (-1 + y) * y) - hx) / exp(gz)
            ))) * sqrt(1 + (1 / (-1 + y) + 1 / (2 * y - 2 *y ^2) + hx
            ) ^ 2 / exp(2 * gz))) - (4 * exp(2 * gz) * (-1 + y) ^ 2 * y ^ 2) / (
              1 + 4 * y * (-1 + hx) + 4 * y ^ 4 * (exp(2 * gz) + hx ^ 2) + 
                4 * y ^ 2 * (1 + exp(2 * gz) - 3 * hx + hx ^ 2) -
                8 * y ^ 3 * (exp(2 * gz) + hx * (-1 + hx))),length(z[1,])))
        colSums(gd, na.rm = TRUE)
       }
       grad <- a_arcsinh_grad
    }

    # arcsinh-burr7
    if (sd == "burr7") {
      a_burr7_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep(((-((atanh(1 - 2 * y) + hx)/(1 + (atanh(1 - 2 * y) + hx)^2/exp(2 * 
          gz))) + exp(gz)/sqrt(1 + (atanh(1 - 2 * y) + hx)^2/exp(2 * gz)) - 
          (2 * exp(asinh((atanh(1 - 2 * y) + hx)/exp(gz)) + gz))/
          ((1 + exp(asinh((atanh(1 - 2 * y) + hx)/exp(gz)))) * sqrt(1 + 
          (atanh(1 - 2 * y) + hx)^2/exp(2 * gz)))))/exp(2 * gz), length(x[1, ])), z * rep((-1 + 
          (atanh(1 - 2 * y) + hx)^2/(exp(2 * gz) * (1 + (atanh(1 - 2 * y) + hx)^2/exp(2 * gz))) - 
          (atanh(1 - 2 * y) + hx)/(exp(gz) * sqrt(1 + (atanh(1 - 2 * y) + hx)^2/exp(2 * gz))) + 
          (2 * exp(asinh((atanh(1 - 2 * y) + hx)/exp(gz)) - gz) * (atanh(1 - 2 * y) + hx))/((1 + 
          exp(asinh((atanh(1 - 2 * y) + hx)/exp(gz)))) * sqrt(1 + (atanh(1 - 2 * y) + hx)^2/
          exp(2 * gz)))), length(z[1, ])))
        
        colSums(gd, na.rm = TRUE)
      }
      grad <- a_burr7_grad
    }
   
    # arcsinh-burr8
    if (sd == "burr8") {
      a_burr8_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep(((-((-log(tan((pi * y)/2)) + hx)/(1 + (log(tan((pi * 
                                                                              y)/2)) - hx)^2/exp(2 * gz))) - exp(gz)/sqrt(1 + (log(tan((pi * 
                                                                                                                                          y)/2)) - hx)^2/exp(2 * gz)) + (2 * exp(asinh((log(tan((pi * y)/2)) - 
                                                                                                                                                                                          hx)/exp(gz)) + gz))/((1 + exp(asinh((log(tan((pi * y)/2)) - hx)/exp(gz)))) * 
                                                                                                                                                                                                                 sqrt(1 + (log(tan((pi * y)/2)) - hx)^2/exp(2 * gz)))))/exp(2 * 
                                                                                                                                                                                                                                                                              gz), length(x[1, ])), z * rep((-1 + (log(tan((pi * y)/2)) - hx)^2/(exp(2 * 
                                                                                                                                                                                                                                                                                                                                                       gz) * (1 + (log(tan((pi * y)/2)) - hx)^2/exp(2 * gz))) + (2 * exp(asinh((log(tan((pi * 
                                                                                                                                                                                                                                                                                                                                                                                                                                           y)/2)) - hx)/exp(gz)) - gz) * (log(tan((pi * y)/2)) - hx))/((1 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          exp(asinh((log(tan((pi * y)/2)) - hx)/exp(gz)))) * sqrt(1 + (log(tan((pi * 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  y)/2)) - hx)^2/exp(2 * gz))) + (-log(tan((pi * y)/2)) + hx)/(exp(gz) * 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 sqrt(1 + (log(tan((pi * y)/2)) - hx)^2/exp(2 * gz)))), length(z[1, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ])))
        
        colSums(gd, na.rm = TRUE)
      }
      grad <- a_burr8_grad
    }
  
    # arcsinh-Cauchy
    if (sd == "cauchy") {
      a_cauchy_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        
        gd <- cbind(x * rep((-((cot(pi * y) + hx)/(exp(2 * gz) + cot(pi * 
                                                                       y)^2 + 2 * cot(pi * y) * hx + hx^2)) + 1/(exp(gz) * sqrt(1 + (cot(pi * 
                                                                                                                                           y) + hx)^2/exp(2 * gz))) - (2 * exp(asinh((cot(pi * y) + hx)/exp(gz)) - 
                                                                                                                                                                                 gz))/((1 + exp(asinh((cot(pi * y) + hx)/exp(gz)))) * sqrt(1 + (cot(pi * 
                                                                                                                                                                                                                                                      y) + hx)^2/exp(2 * gz)))), length(x[1, ])), z * rep(((-(exp(3 * 
                                                                                                                                                                                                                                                                                                                    gz)/(exp(2 * gz) + cot(pi * y)^2 + 2 * cot(pi * y) * hx + hx^2)) - 
                                                                                                                                                                                                                                                                                                              (cot(pi * y) + hx)/sqrt(1 + (cot(pi * y) + hx)^2/exp(2 * gz)) + 
                                                                                                                                                                                                                                                                                                              (2 * exp(asinh((cot(pi * y) + hx)/exp(gz))) * (cot(pi * y) + hx))/((1 + 
                                                                                                                                                                                                                                                                                                                                                                                    exp(asinh((cot(pi * y) + hx)/exp(gz)))) * sqrt(1 + (cot(pi * 
                                                                                                                                                                                                                                                                                                                                                                                                                                              y) + hx)^2/exp(2 * gz)))))/exp(gz), length(z[1, ])))
        colSums(gd, na.rm = TRUE)
      }
      grad <- a_cauchy_grad
    }
    
    # arcsinh-logistic
    if (sd == "logistic") {
      a_logistic_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <-  hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep(((-((log(1 - y) - log(y) + hx)/(1 + (log(1 - y) - log(y) + hx)^2/exp(2 * 
        gz))) - exp(gz)/sqrt(1 + (log(1 - y) - log(y) + hx)^2/exp(2 * gz)) + (2 * exp(asinh((-log(1 - 
        y) + log(y) - hx)/exp(gz)) + gz))/((1 + exp(asinh((-log(1 - y) + log(y) - hx)/exp(gz)))) * 
        sqrt(1 + (log(1 - y) - log(y) + hx)^2/exp(2 * gz)))))/exp(2 * gz), length(x[1, ])), 
        z * rep((-1 + (log(1 - y) - log(y) + hx)^2/(exp(2 * gz) * (1 + (log(1 - y) - log(y) + 
            hx)^2/exp(2 * gz))) + (log(1 - y) - log(y) + hx)/(exp(gz) * sqrt(1 + (log(1 - y) - 
            log(y) + hx)^2/exp(2 * gz))) - (2 * exp(asinh((-log(1 - y) + log(y) - hx)/exp(gz)) - 
            gz) * (log(1 - y) - log(y) + hx))/((1 + exp(asinh((-log(1 - y) + log(y) - hx)/exp(gz)))) * 
            sqrt(1 + (log(1 - y) - log(y) + hx)^2/exp(2 * gz)))), length(z[1, ])))
          
        colSums(gd, na.rm = TRUE)
      }
      grad <- a_logistic_grad
    }
    
    # arcsinh-t2
    if (sd == "t2") {
      a_t2_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        
        a1 <- a2 <- y
        for (i in 1:length(y)) {
          c1 <- ((y[i] >= 0) & (y[i] < 1/2))  #1st situation
          c2 <- (y[i] >= 1/2 & y[i] <= 1)  #2ed situation
          
          if (c1) {
            a1[i] <- -sqrt((1 - 2 * y[i])^2/(2 * (1 - y[i]) * y[i]))
          } else if (c2) {
            a1[i] <- sqrt((1 - 2 * y[i])^2/(2 * (1 - y[i]) * y[i]))
          } else {
            a1[i] <- 0
          }
        }
        
        
        gd <- cbind(x * rep((-(((-a1) * (1 + exp(asinh((a1 - hx)/exp(gz)))) + 
                                  (1 + exp(asinh((a1 - hx)/exp(gz)))) * hx - exp(gz) * (-1 + exp(asinh((a1 - 
                                                                                                          hx)/exp(gz)))) * sqrt((a1^2 + exp(2 * gz) - 2 * a1 * hx + hx^2)/exp(2 * 
                                                                                                                                                                                gz))))/((1 + exp(asinh((a1 - hx)/exp(gz)))) * (a1^2 + exp(2 * gz) - 
                                                                                                                                                                                                                                 2 * a1 * hx + hx^2))), length(x[1, ])), z * rep(-((exp(gz) * (exp(gz) + 
                                                                                                                                                                                                                                                                                                 exp(asinh((a1 - hx)/exp(gz)) + gz) + (a1 - hx) * sqrt((a1^2 + exp(2 * 
                                                                                                                                                                                                                                                                                                                                                                     gz) - 2 * a1 * hx + hx^2)/exp(2 * gz)) + exp(asinh((a1 - hx)/exp(gz))) * 
                                                                                                                                                                                                                                                                                                 (-a1 + hx) * sqrt((a1^2 + exp(2 * gz) - 2 * a1 * hx + hx^2)/exp(2 * 
                                                                                                                                                                                                                                                                                                                                                                   gz))))/((1 + exp(asinh((a1 - hx)/exp(gz)))) * (a1^2 + exp(2 * gz) - 
                                                                                                                                                                                                                                                                                                                                                                                                                    2 * a1 * hx + hx^2))), length(z[1, ])))
        
        
        colSums(gd, na.rm = TRUE)
      }
      grad <- a_t2_grad
    }
    
  }
  
    # burr7-XX-----
  if (fd =="burr7"){
    
    if (sd == 'arcsinh'){
      b7_arcsinh_grad <- function(h, y, x, z) {
    hx <- x %*% h[1:length(x[1, ])]
    gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
    gd <- cbind(x * rep((2 * tanh(((1 - 2 * y)/(2 * (-1 + y) * y) - hx)/exp(gz)))/exp(gz), 
        length(x[1, ])), z * rep((-1 + (2 * tanh(((1 - 2 * y)/(2 * (-1 + y) * y) - hx)/exp(gz)) * 
        ((1 - 2 * y)/(2 * (-1 + y) * y) - hx))/exp(gz)), length(z[1, ])))
      colSums(gd, na.rm = TRUE)
      }
      grad <- b7_arcsinh_grad
    }
    
    if (sd == 'burr7'){
    b7_burr7_grad <- function(h, y, x, z) {
    hx <- x %*% h[1:length(x[1, ])]
    gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
    gd <- cbind(x * rep((-2 * tanh((atanh(1 - 2 * y) + hx)/exp(gz)))/exp(gz), length(x[1, ])), 
        z * rep((-exp(-gz)) * (exp(gz) - 2 * atanh(1 - 2 * y) * tanh((atanh(1 - 2 * y) + hx)/exp(gz)) - 
            2 * hx * tanh((atanh(1 - 2 * y) + hx)/exp(gz))), length(z[1, ])))
    colSums(gd, na.rm = TRUE)
      }
      grad <- b7_burr7_grad
    }
  
    if (sd == 'burr8'){
      b7_burr8_grad <- function(h, y, x, z) {
    hx <- x %*% h[1:length(x[1, ])]
    gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
    # gd <- cbind(x * rep((2 * (1 - (2 * exp((2 * hx)/exp(gz)))/(exp((2 * hx)/exp(gz)) + (-(0 + 
    #     1))^(2/exp(gz)) * ((-1 + exp((0 + 1) * pi * y))/(1 + exp((0 + 1) * pi * y)))^(2/exp(gz)))))/exp(gz), 
    #     length(x[1, ])), z * rep(-(((exp(gz) * (-1 + exp((0 + 1) * pi * y))^(2/exp(gz)) + exp(((0 + 
    #     1) * pi + 2 * hx + exp(gz) * gz)/exp(gz)) * (1 + exp((0 + 1) * pi * y))^(2/exp(gz)) + 
    #     (0 + 1) * (-1 + exp((0 + 1) * pi * y))^(2/exp(gz)) * pi - (0 + 1) * exp(((0 + 1) * 
    #     pi + 2 * hx)/exp(gz)) * (1 + exp((0 + 1) * pi * y))^(2/exp(gz)) * pi + 2 * ((-1 + exp((0 + 
    #     1) * pi * y))^(2/exp(gz)) + exp(((0 + 1) * pi + 2 * hx)/exp(gz)) * (1 + exp((0 + 1) * 
    #     pi * y))^(2/exp(gz))) * log(-1 + exp((0 + 1) * pi * y)) - 2 * ((-1 + exp((0 + 1) * 
    #     pi * y))^(2/exp(gz)) + exp(((0 + 1) * pi + 2 * hx)/exp(gz)) * (1 + exp((0 + 1) * pi * 
    #     y))^(2/exp(gz))) * log(1 + exp((0 + 1) * pi * y)) - 4 * (-1 + exp((0 + 1) * pi * y))^(2/exp(gz)) * 
    #     log(-1 + exp((0 + 1) * pi * y)) + 4 * (-1 + exp((0 + 1) * pi * y))^(2/exp(gz)) * log(1 + 
    #     exp((0 + 1) * pi * y)) + 2 * (-1 + exp((0 + 1) * pi * y))^(2/exp(gz)) * hx - 2 * exp(((0 + 
    #     1) * pi + 2 * hx)/exp(gz)) * (1 + exp((0 + 1) * pi * y))^(2/exp(gz)) * hx))/exp(gz)/((-1 + 
    #     exp((0 + 1) * pi * y))^(2/exp(gz)) + exp(((0 + 1) * pi + 2 * hx)/exp(gz)) * (1 + exp((0 + 
    #     1) * pi * y))^(2/exp(gz)))), length(z[1, ])))
    
    gd <- cbind(x * rep( -2*tanh((hx -log(tan((pi*y)/2)))/exp(gz))/exp(gz), length(x[1, ])),
                z * rep((-exp(gz) + 2*(hx-log(tan((pi*y)/2)))*tanh((hx-log(tan((pi*y)/2)))/exp(gz)))/exp(gz)^2, length(z[1, 
                                                                                                                         ])))
    
    
    colSums(gd, na.rm = TRUE)
      }
      grad <- b7_burr8_grad
    }
    
    if (sd == 'cauchy'){
      b7_cauchy_grad <- function(h, y, x, z) {
    hx <- x %*% h[1:length(x[1, ])]
    gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
    gd <- cbind(x * rep((-2 * tanh((cot(pi * y) + hx)/exp(gz)))/exp(gz), length(x[1, ])), z * 
        rep((-1 + (2 * tanh((cot(pi * y) + hx)/exp(gz)) * (cot(pi * y) + hx))/exp(gz)), length(z[1, 
            ])))
    colSums(gd, na.rm = TRUE)
      }
      grad <- b7_cauchy_grad
    }
  
    
    if (sd == 'logistic'){
      b7_logistic_grad <- function(h, y, x, z) {
    hx <- x %*% h[1:length(x[1, ])]
    gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
    gd <- cbind(x * rep((2 * tanh((-log(1 - y) + log(y) - hx)/exp(gz)))/exp(gz), length(x[1, 
        ])), z * rep((-1 - (2 * tanh((-log(1 - y) + log(y) - hx)/exp(gz)) * (log(1 - y) - log(y) + 
        hx))/exp(gz)), length(z[1, ])))
    colSums(gd, na.rm = TRUE)
      }
      grad <- b7_logistic_grad
    }
    
    if (sd == 't2'){
      b7_t2_grad <- function(h, y, x, z) {
    hx <- x %*% h[1:length(x[1, ])]
    gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
    gd <- cbind(x * rep((2 * tanh(((-1 + 2 * y)/(sqrt(2) * sqrt((-(-1 + y)) * y)) - hx)/exp(gz)))/exp(gz), 
        length(x[1, ])), z * rep(-(((exp(gz) * sqrt((-(-1 + y)) * y) + tanh(((-1 + 2 * y)/(sqrt(2) * 
        sqrt((-(-1 + y)) * y)) - hx)/exp(gz)) * (sqrt(2) - 2 * sqrt(2) * y + 2 * sqrt((-(-1 + 
        y)) * y) * hx)))/(exp(gz) * sqrt((-(-1 + y)) * y))), length(z[1, ])))
    colSums(gd, na.rm = TRUE)
      }
      grad <- b7_t2_grad
    }
    
    
    
  }
  
  
  # burr8-XX-----
  if (fd == "burr8") {
    # burr8 arcsinh
    if (sd == "arcsinh") {
      b8_arcsinh_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep(((exp((1 - 2 * y)/(exp(gz) * ((-1 + y) * y))) - exp((2 * hx)/exp(gz)))/(exp(gz) * 
        (exp((1 - 2 * y)/(exp(gz) * ((-1 + y) * y))) + exp((2 * hx)/exp(gz))))), length(x[1, 
        ])), z * rep(-(((2 * exp(1/(exp(gz) * ((-1 + y) * y)) + gz) * (-1 + y) * y + 2 * exp((2 * 
        (1/(-1 + y) + hx))/exp(gz) + gz) * (-1 + y) * y + exp(1/(exp(gz) * ((-1 + y) * y))) * 
        (-1 + 2 * y^2 * hx + y * (2 - 2 * hx)) + exp((2 + 2 * (-1 + y) * hx)/(exp(gz) * (-1 + 
        y))) * (1 - 2 * y^2 * hx + 2 * y * (-1 + hx))))/exp(gz)/(2 * (exp(1/(exp(gz) * ((-1 + 
        y) * y))) + exp((2 + 2 * (-1 + y) * hx)/(exp(gz) * (-1 + y)))) * (-1 + y) * y)), length(z[1, 
        ])))
        colSums(gd, na.rm = TRUE)
      }
      grad <- b8_arcsinh_grad
    } 
    
    # burr8 burr7
    if (sd == "burr7") {
      b8_burr7_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep((-exp(-gz)) * tanh((atanh(1 - 2 * y) + hx)/exp(gz)), 
                            length(x[1, ])), z * rep((-exp(-gz)) * (exp(gz) - atanh(1 - 2 * 
                                                                                      y) * tanh((atanh(1 - 2 * y) + hx)/exp(gz)) - hx * tanh((atanh(1 - 
                                                                                                                                                      2 * y) + hx)/exp(gz))), length(z[1, ])))
        colSums(gd, na.rm = TRUE)
      }
      grad <- b8_burr7_grad
    }
    
    # burr8 burr8
    if (sd == "burr8") {
      b8_burr8_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep((-1) * (((exp((2 * hx)/exp(gz)) - tan((pi * y)/2)^(2/exp(gz))))/(exp(gz) * 
                                                                                               (exp((2 * hx)/exp(gz)) + tan((pi * y)/2)^(2/exp(gz))))), length(x[1, 
                                                                                                                                                                 ])), z * rep(-(((exp((2 * hx)/exp(gz)) * (exp(gz) + log(tan(pi * 
                                                                                                                                                                                                                               y/2)) - hx) + tan((pi * y)/2)^(2/exp(gz)) * (exp(gz) - 2 * log(tan((pi * 
                                                                                                                                                                                                                                                                                                     y)/2)) + log(tan(pi * y/2)) + hx)))/(exp(gz) * (exp((2 * hx)/exp(gz)) + 
                                                                                                                                                                                                                                                                                                                                                       tan((pi * y)/2)^(2/exp(gz))))), length(z[1, ])))
        colSums(gd, na.rm = TRUE)
      }
      grad <- b8_burr8_grad
    } 
    
    # burr8 Cauchy
    if (sd == "cauchy") {
      b8_cauchy_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep(((-1) * ((-1 + exp((2 * (cot(pi * y) + hx))/exp(gz))))/(exp(gz) * 
                                                                                      (1 + exp((2 * (cot(pi * y) + hx))/exp(gz))))), length(x[1, ])), 
                    z * rep(-(((exp(gz) + exp((2 * cot(pi * y) + 2 * hx + exp(gz) * 
                                                 gz)/exp(gz)) - (-1 + exp((2 * (cot(pi * y) + hx))/exp(gz))) * 
                                  cot(pi * y) + hx - exp((2 * (cot(pi * y) + hx))/exp(gz)) * hx))/(exp(gz) * 
                                                                                                     (1 + exp((2 * (cot(pi * y) + hx))/exp(gz))))), length(z[1, ])))
        
        colSums(gd, na.rm = TRUE)
      }
      grad <- b8_cauchy_grad
    }
    
    
    # burr8 logistic
    if (sd == "logistic") {
      b8_logistic_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep((-exp(-gz)) * (-1 + (2 * exp((2 * hx)/exp(gz)))/(exp((2 * hx)/exp(gz)) + 
        y^(2/exp(gz))/(1 - y)^(2/exp(gz)))), length(x[1, ])), z * rep(-(((exp((2 * hx)/exp(gz) + 
        gz) * (1 - y)^(2/exp(gz)) + exp(gz) * y^(2/exp(gz)) + ((-exp((2 * hx)/exp(gz))) * (1 - 
        y)^(2/exp(gz)) + y^(2/exp(gz))) * log(1 - y) + (exp((2 * hx)/exp(gz)) * (1 - y)^(2/exp(gz)) - 
        y^(2/exp(gz))) * log(y) - exp((2 * hx)/exp(gz)) * (1 - y)^(2/exp(gz)) * hx + y^(2/exp(gz)) * 
        hx))/exp(gz)/(exp((2 * hx)/exp(gz)) * (1 - y)^(2/exp(gz)) + y^(2/exp(gz)))), length(z[1, 
        ])))
        colSums(gd, na.rm = TRUE)
      }
      grad <- b8_logistic_grad
    } 
    
    # burr8 t2
    if (sd == "t2") {
      b8_t2_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        
        for (i in 1:length(y)) {
          c1 <- ((y[i] <= 0) | (y[i] >= 1))  #1st situation
          if (c1) {
            gd <- NA
          }
        }
        gd <- cbind(x * rep((-exp(-gz)) * tanh(((1 - 2 * y)/(sqrt(2) * sqrt((-(-1 + 
                                                                                 y)) * y)) + hx)/exp(gz)), length(x[1, ])), z * rep(-(((2 * exp(gz) * 
                                                                                                                                          sqrt((-(-1 + y)) * y) + tanh(((1 - 2 * y)/(sqrt(2) * sqrt((-(-1 + 
                                                                                                                                                                                                         y)) * y)) + hx)/exp(gz)) * (-sqrt(2) + 2 * sqrt(2) * y - 2 * sqrt((-(-1 + 
                                                                                                                                                                                                                                                                                y)) * y) * hx)))/(exp(gz) * (2 * sqrt((-(-1 + y)) * y)))), length(z[1, 
                                                                                                                                                                                                                                                                                                                                                    ])))
        
        
        colSums(gd, na.rm = TRUE)
      }
      grad <- b8_t2_grad
    }
    
  }
  
  
    #Cauchit-XX-----
  if (fd == "cauchit") {
    # cauchit-arcsinh
    if (sd == "arcsinh") {
       c_arcsinh_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
    
      gd <- cbind(x * rep((-(4 * (-1 + y) * y * (-1 + 2 * y + 2 * (-1 + y) * y * hx))/(4 * exp(2 * 
        gz) * (-1 + y)^2 * y^2 + (-1 + 2 * y + 2 * (-1 + y) * y * hx)^2)), length(x[1, ])), 
        z * rep(1 - (8 * exp(2 * gz) * (-1 + y)^2 * y^2)/(4 * exp(2 * gz) * (-1 + y)^2 * y^2 + 
            (-1 + 2 * y + 2 * (-1 + y) * y * hx)^2), length(z[1, ])))
        colSums(gd, na.rm = TRUE)
       }
       grad <- c_arcsinh_grad
    }  
   
    # cauchit-burr7
    if (sd == "burr7") {
       c_burr7_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
    
      gd <- cbind(x * rep(((-2 * (atanh(1 - 2 * y) + hx))/(exp(2 * gz) + atanh(1 - 2 * y)^2 + 
        2 * atanh(1 - 2 * y) * hx + hx^2)), length(x[1, ])), z * rep((1 - (2 * exp(2 * gz))/(exp(2 * 
        gz) + atanh(1 - 2 * y)^2 + 2 * atanh(1 - 2 * y) * hx + hx^2)), length(z[1, ])))
        colSums(gd, na.rm = TRUE)
       }
       grad <- c_burr7_grad
    }  
    
    # cauchit-burr8
    if (sd == "burr8") {
       c_burr8_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
    
      gd <- cbind(x * rep((2 * (log(tan((pi * y)/2)) - hx))/(exp(2 * gz) + log(tan((pi * y)/2))^2 - 
        2 * hx * log(tan((pi * y)/2)) + hx^2), length(x[1, ])), z * rep((1 - (2 * exp(2 * gz))/(exp(2 * 
        gz) + log(tan((pi * y)/2))^2 - 2 * hx * log(tan((pi * y)/2)) + hx^2)), length(z[1, 
        ])))
        colSums(gd, na.rm = TRUE)
       }
       grad <- c_burr8_grad
    }  
    
    
     # cauchit-Cauchy
    if (sd == "cauchy") {
     c_cauchy_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
    
    gd <- cbind(x * rep(((-2 * (cot(pi * y) + hx))/(exp(2 * gz) + cot(pi * y)^2 + 2 * cot(pi * 
        y) * hx + hx^2)), length(x[1, ])), z * rep((1 - (2 * exp(2 * gz))/(exp(2 * gz) + cot(pi * 
        y)^2 + 2 * cot(pi * y) * hx + hx^2)), length(z[1, ])))
        colSums(gd, na.rm = TRUE)
     } 
     grad <- c_cauchy_grad
    }
  
     # cauchit-Logistic 
    if (sd == "logistic") {
      c_logistic_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
       
        # gd <- cbind(x * rep(((-2 * (log(1 - y) - log(y) + hx))/(exp(2 * gz) + (log(1 - y) - log(y))^2 + 
        # 2 * hx * (log(1 - y) - log(y)) + hx^2)), length(x[1, ])), z * rep((1 - (2 * exp(2 * 
        # gz))/(exp(2 * gz) + (log(1 - y) - log(y))^2 + 2 * hx * (log(1 - y) - log(y)) + hx^2)), 
        # length(z[1, ])))
        
        # gz <- cbind(
        #   x*rep(-log(-(exp(gz)/(pi*(-1 + y)*y*(hx^2 + exp(gz)^2 - 2*hx*log(y/(1 - y)) + log(y/(1 - y))^2)))), length(x[1, ])), 
        #   z*rep(-log(-(exp(gz)/(pi*(-1 + y)*y*(hx^2 + exp(gz)^2 - 2*hx*log(y/(1 - y)) + log(y/(1 - y))^2)))), length(x[1, ]))
        # ) 
      lldq <- function(u) {
          loglik <- -log(dq(u[1], u[2], exp(u[3]), fd = "cauchit", sd = "logistic"))
        } 
      n <- length(y)
      pdfgrad <- matrix(c(rep(0,3*n)), ncol = 3)
      for (i in 1:n) {
      pdfgrad[i,] <- pracma::grad(lldq, x0 <- c(y[i], hx[i], gz[i]))}
      
      pdfgrad <- cbind(
        x*rep(pdfgrad[, 2], length(x[1, ])),
        z*rep(pdfgrad[, 3], length(z[1, ]))
      )
      colSums(pdfgrad, na.rm = TRUE)
      }
      
      grad <- c_logistic_grad
    }
    
     # cauchit-t2 
    if (sd == "t2") {
      c_t2_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        
        a1 <- a2 <- y
        for (i in 1:length(y)) {
        c1 <- ((y[i] >= 0) & (y[i] < 1/2))  #1st situation
        c2 <- (y[i] >= 1/2 & y[i] <= 1)  #2ed situation
        
        if (c1) {
            a1[i] <- -sqrt((1 - 2 * y[i])^2/(2 * (1 - y[i]) * y[i]))
        } else if (c2) {
            a1[i] <- sqrt((1 - 2 * y[i])^2/(2 * (1 - y[i]) * y[i]))
        } else {
            a1[i] <- 0
        }
        }
      gd <- cbind(x * rep((-(2 * (-a1 + hx))/(a1^2 + exp(2 * gz) - 2 * a1 * hx + hx^2)), length(x[1, 
        ])), z * rep(-(((-a1^2 + exp(2 * gz) + 2 * a1 * hx - hx^2))/(a1^2 + exp(2 * gz) - 2 * 
        a1 * hx + hx^2)), length(z[1, ])))
        colSums(gd, na.rm = TRUE)
      }
      grad <- c_t2_grad
    } 
    
  }
  
  
  
  # logit-XX-----
  if (fd == "logit") {
    # logit-arcsinh
    if (sd == "arcsinh") {
      l_arcsinh_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep(((exp(1/(exp(gz) * (1 - y))) - exp((1/(2 * y - 2 * y^2) + hx)/exp(gz))))/(exp(gz) * 
        (exp(1/(exp(gz) * (1 - y))) + exp((1/(2 * y - 2 * y^2) + hx)/exp(gz)))), length(x[1, 
        ])), z * rep(-(((2 * exp(1/(exp(gz) * (1 - y)) + gz) * (-1 + y) * y + 2 * exp((1/(2 * 
        y - 2 * y^2) + hx)/exp(gz) + gz) * (-1 + y) * y + exp(1/(exp(gz) * (1 - y))) * (-1 + 
        2 * y^2 * hx + y * (2 - 2 * hx)) + exp((1/(2 * y - 2 * y^2) + hx)/exp(gz)) * (1 - 2 * 
        y^2 * hx + 2 * y * (-1 + hx))))/exp(gz)/(2 * (exp(1/(exp(gz) * (1 - y))) + exp((1/(2 * 
        y - 2 * y^2) + hx)/exp(gz))) * (-1 + y) * y)), length(z[1, ])))
        colSums(gd, na.rm = TRUE)
      }
      grad <- l_arcsinh_grad
    }
    
    # logit-burr7
    if (sd == "burr7") {
      l_burr7_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep((-exp(-gz)) * tanh(((1/2) * (atanh(1 - 2 * y) + hx))/exp(gz)), length(x[1, 
        ])), z * rep((-exp(-gz)) * (exp(gz) - atanh(1 - 2 * y) * tanh(((1/2) * (atanh(1 - 2 * 
        y) + hx))/exp(gz)) - hx * tanh(((1/2) * (atanh(1 - 2 * y) + hx))/exp(gz))), length(z[1, 
        ])))
        
        colSums(gd, na.rm = TRUE)
      }
      grad <- l_burr7_grad
    }
    
    # logit-burr8
    if (sd == "burr8") {
      l_burr8_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep(((1 - (2 * exp((hx)/exp(gz)))/(exp((hx)/exp(gz)) + 
          tan((pi * y)/2)^exp(-gz))))/exp(gz), length(x[1, ])), z * rep(-(((exp((hx)/exp(gz)) * 
          (exp(gz) + log(tan((pi * y)/2)) - hx) + tan((pi * y)/2)^exp(-gz) * 
          (exp(gz) - log(tan((pi * y)/2)) + hx)))/(exp(gz) * (exp((hx)/exp(gz)) + 
          tan((pi * y)/2)^exp(-gz)))), length(z[1, ])))
        
        colSums(gd, na.rm = TRUE)
      }
      grad <- l_burr8_grad
    } 

    # logit-Cauchy
    if (sd == "cauchy") {
      l_cauchy_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep(((-(-1 + exp((cot(pi * y) + hx)/exp(gz))))/(exp(gz) * 
          (1 + exp((cot(pi * y) + hx)/exp(gz))))), length(x[1, ])), z * rep(-(((exp(gz) + 
          exp((cot(pi * y) + hx + exp(gz) * gz)/exp(gz)) - (-1 + exp((cot(pi * 
          y) + hx)/exp(gz))) * cot(pi * y) + hx - exp((cot(pi * y) + hx)/exp(gz)) * 
          hx))/(exp(gz) * (1 + exp((cot(pi * y) + hx)/exp(gz))))), length(z[1, 
          ])))
        
        colSums(gd, na.rm = TRUE)
      }
      grad <- l_cauchy_grad
    } 

    # logit-logistic
    if (sd == "logistic") {
      l_logistic_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep(((1 - exp((hx)/exp(gz)) * (-1 + 1/y)^exp(-gz)))/(exp(gz) * 
          (1 + exp((hx)/exp(gz)) * (-1 + 1/y)^exp(-gz))), length(x[1, ])), 
          z * rep(((-(1/(1 + exp((hx)/exp(gz)) * (-1 + 1/y)^exp(-gz)))) * 
          (exp(gz) + exp((hx)/exp(gz) + gz) * (-1 + 1/y)^exp(-gz) + (1 - 
            exp((hx)/exp(gz)) * (-1 + 1/y)^exp(-gz)) * log(-1 + 1/y) + 
            hx - exp((hx)/exp(gz)) * (-1 + 1/y)^exp(-gz) * hx))/exp(gz), 
          length(z[1, ])))
        
        colSums(gd, na.rm = TRUE)
      }
      grad <- l_logistic_grad
    }
    
    # logit-t2
    if (sd == "t2") {
      l_t2_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        
        a1 <- a2 <- y
        for (i in 1:length(y)) {
          c1 <- ((y[i] >= 0) & (y[i] < 1/2))  #1st situation
          c2 <- (y[i] >= 1/2 & y[i] <= 1)  #2ed situation
          
          if (c1) {
          a1[i] <- -sqrt((1 - 2 * y[i])^2/(2 * (1 - y[i]) * y[i]))
          } else if (c2) {
          a1[i] <- sqrt((1 - 2 * y[i])^2/(2 * (1 - y[i]) * y[i]))
          } else {
          a1[i] <- 0
          }
          
        }
        
        
        gd <- cbind(x * rep(((exp(a1/exp(gz)) - exp((hx)/exp(gz))))/(exp(gz) * 
          (exp(a1/exp(gz)) + exp((hx)/exp(gz)))), length(x[1, ])), z * rep(-(((exp(a1/exp(gz) + 
          gz) + exp((hx)/exp(gz) + gz) + a1 * (-exp(a1/exp(gz)) + exp((hx)/exp(gz))) + 
          exp(a1/exp(gz)) * hx - exp((hx)/exp(gz)) * hx))/(exp(gz) * (exp(a1/exp(gz)) + 
          exp((hx)/exp(gz))))), length(z[1, ])))
        
        colSums(gd, na.rm = TRUE)
      }
      grad <- l_t2_grad
    }
    

  }
  
 
  # t2-XX-----
  if (fd == "t2") {
    # t2-ArcSinh 
    if (sd == "arcsinh") {
        t2_arcsinh_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep((-(6 * (-1 + y) * y * (-1 + 2 * y^2 * hx + y * (2 - 2 * hx)))/(1 + 
        4 * y * (-1 + hx) + 4 * y^4 * (2 * exp(2 * gz) + hx^2) + 4 * y^2 * (1 + 2 * exp(2 * 
        gz) - 3 * hx + hx^2) - 8 * y^3 * (2 * exp(2 * gz) + hx * (-1 + hx)))), length(x[1, 
        ])), z * rep(2 * (1 - (12 * exp(2 * gz) * (-1 + y)^2 * y^2)/(1 + 4 * y * (-1 + hx) + 
        4 * y^4 * (2 * exp(2 * gz) + hx^2) + 4 * y^2 * (1 + 2 * exp(2 * gz) - 3 * hx + hx^2) - 
        8 * y^3 * (2 * exp(2 * gz) + hx * (-1 + hx)))), length(z[1, ])))
        
        colSums(gd, na.rm = TRUE)
        }
        grad <- t2_arcsinh_grad
    }
    
    # t2-burr7
    if (sd == "burr7") {
      t2_burr7_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep((-(3 * (atanh(1 - 2 * y) + hx))/(2 * exp(2 * 
          gz) + atanh(1 - 2 * y)^2 + 2 * atanh(1 - 2 * y) * hx + hx^2)), 
          length(x[1, ])), z * rep(2 * (1 - (3 * exp(2 * gz))/(2 * exp(2 * 
          gz) + atanh(1 - 2 * y)^2 + 2 * atanh(1 - 2 * y) * hx + hx^2)), 
          length(z[1, ])))
        
        colSums(gd, na.rm = TRUE)
      }
      grad <- t2_burr7_grad
    }
    
    # t2-burr8
    if (sd == "burr8") {
      t2_burr8_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep((3 * (log(tan((pi * y)/2)) - hx))/(2 * exp(2 * 
          gz) + log(tan((pi * y)/2))^2 - 2 * hx * log(tan((pi * y)/2)) + 
          hx^2), length(x[1, ])), z * rep(2 * (1 - (3 * exp(2 * gz))/(2 * 
          exp(2 * gz) + log(tan((pi * y)/2))^2 - 2 * hx * log(tan((pi * y)/2)) + 
          hx^2)), length(z[1, ])))
        
        colSums(gd, na.rm = TRUE)
      }
      grad <- t2_burr8_grad
      
    } 
  
    # t2-Cauchy
    if (sd == "cauchy") {
      t2_cauchy_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep((-(3 * sqrt(y) * (-sqrt(2 - 2 * y) + 2 * sqrt(2 - 
          2 * y) * y - 2 * sqrt(y) * hx + 2 * y^(3/2) * hx))/(-1 + 4 * sqrt(2 - 
          2 * y) * y^(3/2) * hx - 2 * sqrt(2) * sqrt((-(-1 + y)) * y) * hx - 
          2 * y * (-2 + 2 * exp(2 * gz) + hx^2) + 2 * y^2 * (-2 + 2 * exp(2 * 
          gz) + hx^2))), length(x[1, ])), z * rep(2 - (12 * exp(2 * gz) * 
          (-1 + y) * y)/(-1 + 4 * sqrt(2 - 2 * y) * y^(3/2) * hx - 2 * sqrt(2) * 
          sqrt((-(-1 + y)) * y) * hx - 2 * y * (-2 + 2 * exp(2 * gz) + hx^2) + 
          2 * y^2 * (-2 + 2 * exp(2 * gz) + hx^2)), length(z[1, ])))
        
        colSums(gd, na.rm = TRUE)
      }
      grad <- t2_cauchy_grad
    }
    
    # t2-t2
    if (sd == "t2") {
      t2_t2_grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep((-(3 * sqrt(y) * (-sqrt(2 - 2 * y) + 2 * sqrt(2 - 
          2 * y) * y - 2 * sqrt(y) * hx + 2 * y^(3/2) * hx))/(-1 + 4 * sqrt(2 - 
          2 * y) * y^(3/2) * hx - 2 * sqrt(2) * sqrt((-(-1 + y)) * y) * hx - 
          2 * y * (-2 + 2 * exp(2 * gz) + hx^2) + 2 * y^2 * (-2 + 2 * exp(2 * 
          gz) + hx^2))), length(x[1, ])), z * rep(2 - (12 * exp(2 * gz) * 
          (-1 + y) * y)/(-1 + 4 * sqrt(2 - 2 * y) * y^(3/2) * hx - 2 * sqrt(2) * 
          sqrt((-(-1 + y)) * y) * hx - 2 * y * (-2 + 2 * exp(2 * gz) + hx^2) + 
          2 * y^2 * (-2 + 2 * exp(2 * gz) + hx^2)), length(z[1, ])))
        
        colSums(gd, na.rm = TRUE)
      }
      grad <- t2_t2_grad
    }
    
    # t2-Cauchy
    if (sd == "logistic") {
      t2_logistic_grad <- function(h, y, x, z) {
   hx = x%*%h[1: length(x[1,])]
  gz = z%*%h[length(x[1,])+1:length(z[1,])]
  gd <- cbind(x*rep((-(3*(log(1 - y) - log(y) + hx))/(2*exp(2*gz) + 
      log(1 - y)^2 + log(y)^2 - 2*hx*log(y) + hx^2 -2*log(1 - y)*(log(y) - hx))), 
    length(x[1,])), z*rep(-((2*(exp(2*gz) - log(1 - y)^2 - 
        log(y)^2 + 2*hx*log(y) - hx^2 + 2*log(1 - y)*(log(y) - hx)))/(2*exp(2*gz) + 
            log(1 - y)^2 + log(y)^2 - 2*hx*log(y) + hx^2 - 2*log(1 - y)*(log(y) - hx))), length(z[1,])))
  

        colSums(gd, na.rm = TRUE)
      }
      grad <- t2_logistic_grad
    }  

  }
  
 #KM KM-----
  if (fd == "km" | sd == "km") {
    km_grad <- function(h, y, x, z) {
      hx <- x %*% h[1:length(x[1, ])]
      gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
      gd <- cbind(x * rep(((-1 + y^exp(hx) + exp(hx) * (-1 + exp(gz) * y^exp(hx)) * 
        log(y)))/(-1 + y^exp(hx)), length(x[1, ])), z * rep((1 + exp(gz) * 
        log(1 - y^exp(hx))), length(z[1, ])))
      colSums(gd, na.rm = TRUE)
    }
    grad <- km_grad
  }
  
  
  grad
} 
