#' @title Function to Give the Log Likelihood Function
#' @aliases qrLogLikFun
#' @description Function to compute the (negative) log likelihood for fitting cdfquantile models.
#' 
#' @param fd A string that specifies the parent distribution.
#' @param sd A string that specifies the sub-family distribution.
#'
#' @return The log-likelihood calculation function given a specified cdfquantile distribution.
#' 
#' @export
#' @examples
#' qrLogLikFun('t2','t2')
#' 
qrLogLikFun <- function(fd, sd) {
  
  # arcsinh-XX-----
  if (fd == "arcsinh") {
    
    # arcsinh-arcsinh
    if (sd == "arcsinh") {
      a_arcsinh_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- asinh(((1 - 2 * y)/(2 * (-1 + y) * y) - mu)/sigma) - 
          2 * log(1 + exp(asinh(((1 - 2 * y)/(2 * (-1 + y) * y) - mu)/sigma))) - 
          log(1 - y) - log(y) + log(1 - 2 * y + 2 * y^2) - 
          (1/2) * log(1 + 4 * y * (-1 + mu) + 4 * y^4 * (mu^2 + sigma^2) + 
                        4 * y^2 * (1 - 3 * mu + mu^2 + sigma^2) - 
                        8 * y^3 * (-mu + mu^2 + sigma^2))
        
        -sum(loglik, na.rm = TRUE)
      }
      quantreg_loglik <- a_arcsinh_loglik
    }
    
    # arcsinh-burr7
    if (sd == "burr7") {
      a_b7_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- -log(2) - log(sigma) + asinh((mu + atanh(1 - 2 * y))/sigma) - 
      2 * log(1 + exp(asinh((mu + atanh(1 - 2 * y))/sigma))) - 
      (1/2) * log(1 + (mu + atanh(1 - 2 * y))^2/sigma^2) - log(1 - y) - log(y)
        
        -sum(loglik, na.rm = TRUE)
      }
       quantreg_loglik <- a_b7_loglik
    }
    
    # arcsinh-burr8
    if (sd == "burr8") {
      a_b8_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- log(pi) - log(sigma) + asinh((-mu + log(tan((pi * y)/2)))/sigma) - 
          2 * log(1 + exp(asinh((-mu + log(tan((pi * y)/2)))/sigma))) + 
          log(csc(pi * y)) - (1/2) * log(1 +(mu - log(tan((pi * y)/2)))^2/sigma^2)
        
        -sum(loglik, na.rm = TRUE)
      }
       quantreg_loglik <- a_b8_loglik
    }  
    
    # arcsinh-Cauchy
    if (sd == "cauchy") {
      a_cauchy_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- log(pi) + asinh((mu + cot(pi * y))/sigma) - 
          2 * log(1 + exp(asinh((mu + cot(pi *   y))/sigma))) - 
          (1/2) * log(mu^2 + sigma^2 + 2 * mu * cot(pi * y) + cot(pi * y)^2) + 
          2 * log(csc(pi * y))
        
        -sum(loglik, na.rm = TRUE)
      }
      quantreg_loglik <- a_cauchy_loglik
    }
    
    # arcsinh-logisitc
    if (sd == "logistic") {
      a_logistic_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <-  -log(sigma) + asinh((-mu - log(1 - y) + log(y))/sigma) -
          2 * log(1 + exp(asinh((-mu - log(1 - y) + log(y))/sigma))) - 
          (1/2) * log(1 + (mu + log(1 - y) - log(y))^2/sigma^2) - log(1 - y) - log(y)
        
        -sum(loglik, na.rm = TRUE)
      }
      quantreg_loglik <- a_logistic_loglik
    }

    # arcsinh-t2
    if (sd == "t2") {
      a_t2_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
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
          
          
          c3 <- (y[i] > 0 & y[i] < 1/2)  #3rd situation
          c4 <- (y[i] > 1/2 & y[i] < 1)  #4th situation
          c5 <- (y[i] == 1 | y[i] == 0)
          if (y[i] == 0.5) {
            a2[i] <- 2 * sqrt(2)
          } else if (c3) {
            a2[i] <- (1 - 2 * y[i])/sqrt(8 * (1 - 2 * y[i])^2 * ((1 - y[i]) * y[i])^3)
          } else if (c4) {
            a2[i] <- (2 * y[i] - 1)/sqrt(8 * (1 - 2 * y[i])^2 * ((1 - y[i]) * y[i])^3)
          } else if (c5) {
            a2[i] <- NA  # NO IDEA    
          } else {
            a2[i] <- 0
          }
        }
        
        lik <- (exp(asinh((-mu + a1)/sigma)) * a2)/((1 + exp(asinh((-mu + a1)/sigma)))^2 * sigma * 
                                                      sqrt(1 + (mu - a1)^2/sigma^2))
        loglik <- log(lik)
        -sum(loglik, na.rm = TRUE)
      }
      quantreg_loglik <- a_t2_loglik
    }
  }
  
  # burr7-XX-----
  
  if (fd =="burr7"){
    if (sd == 'arcsinh'){
       b7_arcsinh_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
       loglik <- -log(4) - log(sigma) + 2 * log(sech((-mu - (1 - 2 * y)/(2 * (1 - y) * y))/sigma)) - 
        2 * log(1 - y) - 2 * log(y) + log(1 - 2 * (1 - y) * y)
        
        -sum(loglik, na.rm = TRUE)
       }
       quantreg_loglik <- b7_arcsinh_loglik
    } 

    if (sd == 'burr7'){
       b7_burr7_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- 2 * log(sech((mu + atanh(1 - 2 * y))/sigma)) - log(4 * sigma * y - 4 * sigma * 
        y^2)
        
        -sum(loglik, na.rm = TRUE)
       }
       quantreg_loglik <- b7_burr7_loglik
    }
    
    if (sd == 'burr8'){
       b7_burr8_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
       loglik <- Re((1/(2 * sigma)) * (-2 * (0 + 1i) * pi + 4 * mu + (0 + 1i) * pi * sigma + sigma * 
        log(16) + 2 * sigma * log(pi) - 2 * sigma * log(sigma) + 4 * log(-1 + exp((0 + 1i) * 
        pi * y)) - 4 * log(1 + exp((0 + 1i) * pi * y)) - 2 * sigma * log(-1 + exp(2 * (0 + 1i) * 
        pi * y)) - 4 * sigma * log(exp((2 * mu)/sigma) + ((-(0 + 1i))^(2/sigma) * (-1 + exp((0 + 
        1) * pi * y))^(2/sigma))/(1 + exp((0 + 1i) * pi * y))^(2/sigma)) + 2 * (0 + 1i) * pi * 
        sigma * y))
        
        -sum(loglik, na.rm = TRUE)
       }
       quantreg_loglik <- b7_burr8_loglik
    }
    
    if (sd == 'cauchy'){
       b7_cauchy_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- -log((2 * sigma)/pi) + 2 * log(csc(pi * y)) + 2 * log(sech((mu + cot(pi * y))/sigma))
        
        -sum(loglik, na.rm = TRUE)
       }
       quantreg_loglik <- b7_cauchy_loglik
    }
    
    if (sd == 'logistic'){
       b7_logistic_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
       loglik <- -log(2) - log(sigma) + 2 * log(sech((-mu - log(1 - y) + log(y))/sigma)) - log(1 - 
        y) - log(y)
        
        -sum(loglik, na.rm = TRUE)
       }
       quantreg_loglik <- b7_logistic_loglik
    }
    
     if (sd == 't2'){
       b7_t2_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        
        a1 <- y
    for (i in 1:length(y)) {
        c1 <- ((y[i] >= 0) & (y[i] < 1))  #1st situation
        c2 <- (y[i] == 0 & y[i] == 1)  #2ed situation
        
        if (c1) {
            a1[i] <- sech(((-1 + 2 * y[i])/(sqrt(2) * sqrt((-(-1 + y[i])) * y[i])) - mu[i])/sigma[i])^2/(4 * 
                sqrt(2) * ((-(-1 + y[i])) * y[i])^(3/2) * sigma[i])
            # sech((-mu[i] + (-1 + 2*y[i])/(sqrt(2)*sqrt(1 - y[i])*sqrt(y[i])))/sigma[i])^2/
            # (4*sqrt(2)*sigma[i]*(1 - y[i])^(3/2)*y[i]^(3/2))
            
        } else if (c2) {
            a1[i] <- NA
        } else {
            a1[i] <- 0
        }
    }
    
    loglik <- log(a1)
        
        -sum(loglik, na.rm = TRUE)
       }
       quantreg_loglik <- b7_t2_loglik
    }
  }
  
  # burr8-XX-----
  if (fd == "burr8") {
    
    if (sd == 'arcsinh'){
       b8_arcsinh_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
       loglik <- (-(1/(2 * sigma))) * (2 * sigma * log(pi) + 2 * sigma * log(sigma) + 2 * sigma * 
        log(exp((2 * (mu - 1/(1 - y)))/sigma) + exp(1/(sigma * (-1 + y) * y))) + 4 * sigma * 
        log(1 - y) + 4 * sigma * log(y) - 2 * sigma * log(1 - 2 * (1 - y) * y) + 2 * (1/(1 - 
        y)) - 2 * mu * (1/(1 - y)) + 1/((1 - y) * y) + 2 * mu * (y/(1 - y)))
        
        -sum(loglik, na.rm = TRUE)
       }
       quantreg_loglik <- b8_arcsinh_loglik
    } 
    
    # burr8 burr7
    if (sd == "burr7") {
      b8_burr7_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- -log(2 * pi * y * sigma - 2 * pi * y^2 * sigma) + 
          log(sech((mu + atanh(1 - 2 *  y))/sigma))
        -sum(loglik, na.rm = TRUE)
      }
      quantreg_loglik <- b8_burr7_loglik
    }
    
    # burr8 burr8
    if (sd == "burr8") {
      b8_burr8_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- (mu + sigma * log(2) + sigma * log(csc(pi * y)) + log(tan(pi * y/2)) - sigma * 
                     log(exp((2 * mu)/sigma) * sigma + sigma * tan((pi * y)/2)^(2/sigma)))/sigma
        -sum(loglik, na.rm = TRUE)
      }
      quantreg_loglik <- b8_burr8_loglik
    }
    
    # burr8 Cauchy
    if (sd == "cauchy") {
      b8_cauchy_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- (mu + sigma * log(2) - sigma * log(sigma) + cot(pi * y) - 
                     sigma * log(1 + exp((2 * (mu + cot(pi * y)))/sigma)) + 
                     2 * sigma * log(csc(pi * y)))/sigma
        
        -sum(loglik, na.rm = TRUE)
      }
       quantreg_loglik <- b8_cauchy_loglik
    }
    
    # burr8 logistic
    if (sd == 'logistic'){
       b8_logistic_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- -((-mu - sigma * log(2) + sigma * log(pi) + sigma * log(sigma) + (1 + sigma) * 
        log(1 - y) + (-1 + sigma) * log(y) + sigma * log(exp((2 * mu)/sigma) + y^(2/sigma)/(1 - 
        y)^(2/sigma)))/sigma)
        
        -sum(loglik, na.rm = TRUE)
       }
       quantreg_loglik <- b8_logistic_loglik
    }
    
    # burr8 t2
    if (sd == "t2") {
      b8_t2_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        a1 <- a2 <- y
        
        for (i in 1:length(y)) {
          c1 <- ((y[i] >= 0) & (y[i] < 1))  #1st situation
          c2 <- (y[i] == 0 & y[i] == 1)  #2ed situation
          
          if (c1) {
            a1[i] <- sech(((1 - 2 * y[i])/(sqrt(2) * sqrt((-(-1 + y[i])) * y[i])) + mu[i])/sigma[i])/(2 * sqrt(2) * pi * ((-(-1 + y[i])) * y[i])^(3/2) * sigma[i])
            # sech((mu[i] + (1 - 2*y[i])/(sqrt(2)*sqrt(1 - y[i])*sqrt(y[i])))/sigma[i])/
            # (2*sqrt(2)*pi*sigma[i]*(1 - y[i])^(3/2)*y[i]^(3/2))
          } else if (c2) {
            a1[i] <- NA
          } else {
            a1[i] <- 0
          }
        }
        
        loglik <- log(a1)
        -sum(loglik, na.rm = TRUE)
      }
      quantreg_loglik <- b8_t2_loglik
    }
  }
  
  #Cauchit-XX-----
  if (fd == "cauchit") {
    # cauchit-arcsinh
    if (sd == "arcsinh") {
      c_arcsinh_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- log((2 * sigma)/pi) + log(1 + 2 * (-1 + y) * y) - log(4 * sigma^2 * (-1 + y)^2 * 
        y^2 + (-1 + 2 * y + 2 * mu * (-1 + y) * y)^2)
        -sum(loglik, na.rm = TRUE)
      }
      quantreg_loglik <- c_arcsinh_loglik
    } 
    
    # cauchit-BURR7
    if (sd == "burr7") {
      c_burr7_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- -log(2) - log(pi) + log(sigma) - log(mu^2 + sigma^2 + 2 * mu * atanh(1 - 2 * 
        y) + atanh(1 - 2 * y)^2) - log(1 - y) - log(y)
        
        -sum(loglik, na.rm = TRUE)
      }
      quantreg_loglik <- c_burr7_loglik
    } 
    
    # cauchit-BURR8
    if (sd == "burr8") {
      c_burr8_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- log(sigma) + log(csc(pi * y)) - log(mu^2 + sigma^2 - 2 * mu * log(tan((pi * y)/2)) + 
        log(tan((pi * y)/2))^2)
        
        -sum(loglik, na.rm = TRUE)
      }
      quantreg_loglik <- c_burr8_loglik
    } 
     
    # cauchit-Cauchy
    if (sd == "cauchy") {
      c_cauchy_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- log(sigma) - log(mu^2 + sigma^2 + 2 * mu * cot(pi * y) + cot(pi * y)^2) + 2 * 
        log(csc(pi * y))
        -sum(loglik, na.rm = TRUE)
      }
      quantreg_loglik <- c_cauchy_loglik
    }
  
    # cauchit-Logistic 
    if (sd == "logistic") {
      c_logistic_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- log(sigma) - log(mu^2 + sigma^2 - 2 * mu * (-log(1 - y) + log(y)) + (-log(1 - 
        y) + log(y))^2) - log(1 - y) - log(y)
        -sum(loglik, na.rm = TRUE)
      }
      quantreg_loglik <- c_logistic_loglik
    }
    
     # cauchit-t2 
    if (sd == "t2") {
      c_t2_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        
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
        
        
        c3 <- (y[i] > 0 & y[i] < 1/2)  #3rd situation
        c4 <- (y[i] > 1/2 & y[i] < 1)  #4th situation
        c5 <- (y[i] == 1 | y[i] == 0)
        if (y[i] == 0.5) {
            a2[i] <- 2 * sqrt(2)
        } else if (c3) {
            a2[i] <- (1 - 2 * y[i])/sqrt(8 * (1 - 2 * y[i])^2 * ((1 - y[i]) * y[i])^3)
        } else if (c4) {
            a2[i] <- (2 * y[i] - 1)/sqrt(8 * (1 - 2 * y[i])^2 * ((1 - y[i]) * y[i])^3)
        } else if (c5) {
            a2[i] <- NA  # NO IDEA    
        } else {
            a2[i] <- 0
        }
    }
    
    
    lik <- a2/(pi * sigma * (1 + (mu - a1)^2/sigma^2))
    loglik <- log(lik)
    -sum(loglik, na.rm = TRUE)
      }
      quantreg_loglik <- c_t2_loglik
    }
    
  }
  
  #logit-XX-----
  if (fd == "logit") {
    
  # logit-arcsinh
    if (sd == "arcsinh") {
      l_arcsinh_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
     sigma <- exp(gz)
      loglik <- -log(2) - log(sigma) - 2 * log(exp(1/(sigma - sigma * y)) + exp(mu/sigma + 1/(2 * 
        sigma * y - 2 * sigma * y^2))) - 2 * log(1 - y) - 2 * log(y) + log(1 - 2 * y + 2 * 
        y^2) + 1/(2 * sigma * y - 2 * sigma * y^2) + 2 * (y/(2 * sigma * y - 2 * sigma * y^2)) + 
        2 * mu * (y/(2 * sigma * y - 2 * sigma * y^2)) - 2 * mu * (y^2/(2 * sigma * y - 2 * 
        sigma * y^2))
        
        -sum(loglik, na.rm = TRUE)
      }
      quantreg_loglik <- l_arcsinh_loglik
    }    
    
    # logit-burr7
    if (sd == "burr7") {
      l_burr7_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- 2 * log(sech((mu + atanh(1 - 2 * y))/(2 * sigma))) - 
          log(8 * sigma * y - 8 * sigma * y^2)
        -sum(loglik, na.rm = TRUE)
      }
      quantreg_loglik <- l_burr7_loglik
    }
    
    # logit-burr8
    if (sd == "burr8") {
      l_burr8_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <-log(csc(pi * y)) + ((mu + sigma * log(pi) - sigma * log(sigma)) 
                  + log(tan((pi *   y)/2)) - 2 * sigma * log(exp(mu/sigma) + 
                  tan((pi * y)/2)^(1/sigma)))/sigma
        -sum(loglik, na.rm = TRUE)
      }
      quantreg_loglik <- l_burr8_loglik
    }

    # logit-Cauchy
    if (sd == "cauchy") {
      l_cauchy_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- (mu + sigma * log(pi) - sigma * log(sigma) + cot(pi * y) - 
                     2 * sigma * log(1 + exp((mu + cot(pi * y))/sigma)) + 
                     2 * sigma * log(csc(pi * y)))/sigma
        -sum(loglik, na.rm = TRUE)
      }
      quantreg_loglik <- l_cauchy_loglik
    }
    
    # logit-logistic
    if (sd == "logistic") {
      l_logistic_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- ((mu - sigma * log(sigma)) - (-1 + sigma) * log(-1 + 1/y) - 
          2 * sigma * log(y + exp(mu/sigma) * (-1 + 1/y)^(1/sigma) * y))/sigma
        -sum(loglik, na.rm = TRUE)
      }
      quantreg_loglik <- l_logistic_loglik
    }
  
    # logit-t2
    if (sd == "t2") {
      l_t2_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
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
          
          
          c3 <- (y[i] > 0 & y[i] < 1/2)  #3rd situation
          c4 <- (y[i] > 1/2 & y[i] < 1)  #4th situation
          c5 <- (y[i] == 1 | y[i] == 0)
          if (y[i] == 0.5) {
            a2[i] <- 2 * sqrt(2)
          } else if (c3) {
            a2[i] <- (1 - 2 * y[i])/sqrt(8 * (1 - 2 * y[i])^2 * ((1 - y[i]) * y[i])^3)
          } else if (c4) {
            a2[i] <- (2 * y[i] - 1)/sqrt(8 * (1 - 2 * y[i])^2 * ((1 - y[i]) * y[i])^3)
          } else if (c5) {
            a2[i] <- NA  # NO IDEA    
          } else {
            a2[i] <- 0
          }
        }
        
        lik <- (exp((mu + a1)/sigma) * a2)/((exp(mu/sigma) + exp(a1/sigma))^2 * sigma)
        loglik <- log(lik)
        -sum(loglik, na.rm = TRUE)
      }
      quantreg_loglik <- l_t2_loglik
    }
    

  }
  
  
  # t2-XX-----
  if (fd == "t2") {
     # t2-ArcSinh 
    if (sd == "arcsinh") {
      t2_arcsinh_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
       loglik <- log(4) + 2 * log(sigma) + 3 * log(1 - y) + log(y) - (3/2) * log(1 + 4 * (-1 + 
        mu) * y + 4 * (1 - 3 * mu + mu^2 + 2 * sigma^2) * y^2 - 8 * (-mu + mu^2 + 2 * sigma^2) * 
        y^3 + 4 * (mu^2 + 2 * sigma^2) * y^4)
        -sum(loglik, na.rm = TRUE)
      }
      quantreg_loglik <- t2_arcsinh_loglik
    }

    # t2-burr7
    if (sd == "burr7") {
      t2_burr7_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- -log(2) + 2 * log(sigma) - (3/2) * log(mu^2 + 2 * sigma^2 + 
                   2 * mu * atanh(1 - 2 * y) + atanh(1 - 2 * y)^2) - log(1 - y) - log(y)
        
        -sum(loglik, na.rm = TRUE)
      }
      quantreg_loglik <- t2_burr7_loglik
    }
    
    # t2-burr8
    if (sd == "burr8") {
      t2_burr8_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- (log(pi) + 2 * log(sigma)) + log(csc(pi * y)) - 
          (3/2) * log(mu^2 + 2 * sigma^2 - 2 * mu * log(tan((pi * y)/2)) + log(tan((pi * y)/2))^2)
        -sum(loglik, na.rm = TRUE)
      }
      quantreg_loglik <- t2_burr8_loglik
      
    } 
    
    # t2-Cauchy
    if (sd == "cauchy") {
      t2_cauchy_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- log(pi) + log(sigma) - log(mu^2 + 2 * sigma^2 + 2 * mu * cot(pi * y)
                   + cot(pi * y)^2) - (1/2) * log(2 + (mu + cot(pi * y))^2/sigma^2) + 2 * log(csc(pi * y))
        -sum(loglik, na.rm = TRUE)
      }
      quantreg_loglik <- t2_cauchy_loglik
    }
    
    # t2-t2
    if (sd == "t2") {
      t2_t2_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- 2 * log(sigma) - (3/2) * log(1 + 2 * sqrt(2) * mu * sqrt(1 - y) * sqrt(y) + 2 * 
                   (-2 + mu^2 + 2 * sigma^2) * y - 
                  4 * mu * sqrt(2 - 2 * y) * y^(3/2) - 2 * (-2 + mu^2 + 2 * sigma^2) * y^2)
        -sum(loglik, na.rm = TRUE)
      }
      quantreg_loglik <- t2_t2_loglik
    }
    
    # t2-logistic
    if (sd == "logistic") {
      t2_logistic_loglik <- function(h, y, x, z) {
  hx = x%*%h[1: length(x[1,])]
  mu = hx
  gz = z%*%h[length(x[1,])+1: length(z[1,])]
  sigma = exp(gz)
  loglik =log(sigma) -1/2 * log(2 + (mu + log(1 - y) - log(y))^2/sigma^2) -log(mu^2 + 2* sigma^2 -2* mu *(-log(1 - y) + log(y)) + (-log(1 - y) + log(y))^2) -log(1 - y) - log(y)
  
        -sum(loglik, na.rm = TRUE)
      }
      quantreg_loglik <- t2_logistic_loglik
    }  

  }
  
 
  # km-----
  if (fd == "km" | sd == "km") {
    km_km_loglik <- function(h, y, x, z) {
      hx <- x %*% h[1:length(x[1, ])]
      a <- exp(hx)
      gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
      b <- exp(gz)
      loglik <- log(a) + log(b) + (a - 1) * log(y) + (b - 1) * log(1 - y^a)
      -sum(loglik, na.rm = TRUE)
    }
    quantreg_loglik <- km_km_loglik
  }
  
  
  quantreg_loglik
} 
