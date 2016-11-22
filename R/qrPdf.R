#' @title The Famility of Distributions 
#' @aliases dq rq qq pq
#' 
#' @description Density function, distribution function, quantile function, and random generation of variates for a specified cdf-quantile distribution with \code{mean} equal to mean and standard deviation equal to \code{sd}.
#' 
#' @param x vector of quantiles. 
#' @param p vector of probabilities.
#' @param n Number of random samples.
#' @param q vector of quantiles. 
#' @param mu vector of means.
#' @param sigma vector of standard deviations.
#' @param fd A string that specifies the parent distribution.
#' @param sd A string that specifies the sub-family distribution.
#' @export
#' @import pracma
#' 
#' @return \code{dq} gives the density, \code{rq} generates random variates, \code{qq} gives the quantile function, and \code{pq} gives the cumulative densitity of specified distribution.
#' @examples
#' x <- rq(5, mu = 0.5, sigma = 1, 't2','t2'); x
#'
#' dq(x, mu = 0.5, sigma = 1, 't2','t2')
#' 
#' qtil <- pq(x, mu = 0.5, sigma = 1, 't2','t2');qtil
#' 
#' qq(qtil , mu = 0.5, sigma = 1, 't2','t2')
#' 

dq <- function(x, mu, sigma, fd, sd) {
  
  fd <- tolower(fd)
  sd <- tolower(sd)
  
  if (any(x <= 0) | any(x >= 1)) 
    stop(paste("x must be between 0 and 1", "\n", ""))
 
  # arcsinh-XX-----
  if (fd == "arcsinh") {
    if (sd == "arcsinh") {
      dent <- function(x, mu, sigma) 
        ((1 + 2*(-1 + x)*x)*(1 + x*(-2 - 2*(-1 + x)*mu + (-1 + x)*sqrt(4 + (-1 + 2*x + 2*(-1 + x)*x*mu)^2/((-1 + x)^2*x^2*sigma^2))*sigma)))/
        ((-1 + x)*x*sqrt(1 + ((1/2)*(1/(-1 + x) + 1/x) + mu)^2/sigma^2)*
           (1 + x*(-2 - 2*(-1 + x)*mu + (-1 + x)*
                     (2 + sqrt(4 + (-1 + 2*x + 2*(-1 + x)*x*mu)^2/
                                 ((-1 + x)^2*x^2*sigma^2)))*sigma))^2)
        
        }

    if (sd == "burr7") {
      dent <- function(x, mu, sigma) 
          -(exp(asinh((mu + atanh(1 - 2*x))/sigma))/(2*(1 + exp(asinh((mu + 
              atanh(1 - 2*x))/sigma)))^2*(-1 + x)* x*sigma* sqrt(1 + 
                  (mu + atanh(1 - 2*x))^2/sigma^2)))

    }  
    
    if (sd == "burr8") {
      dent <- function(x, mu, sigma)
        (exp(asinh((-mu + log(tan((pi * x)/2)))/sigma)) * pi * 
            csc(pi * x))/((1 + exp(asinh((-mu + log(tan((pi * x)/2)))/sigma)))^2 * 
            sigma * sqrt(1 + (mu - log(tan((pi * x)/2)))^2/sigma^2))
    }  
    
    if (sd == "cauchy") {
      dent <- function(x, mu, sigma) 
        (exp(asinh((mu + cot(pi * x))/sigma)) * pi * 
        csc(pi * x)^2)/((1 + exp(asinh((mu + cot(pi * x))/sigma)))^2 * 
        sigma * sqrt((mu^2 + sigma^2 + 2 * mu * cot(pi * x) + cot(pi * x)^2)/sigma^2))
    }
  
    if (sd == "logistic") {
      dent <- function(x, mu, sigma) 
        -(exp(asinh((-mu + log(x/(1 - x)))/sigma))/((1 + 
            exp(asinh((-mu + log(x/(1 - x)))/sigma)))^2*(-1 + x)*x*sigma* 
            sqrt(1 + (mu - log(x/(1 - x)))^2/sigma^2)))
    }
    
    if (sd == "t2") {
      dent <- function(x, mu, sigma) {
        c1 <- ((x >= 0) & (x < 1/2))  #1st situation
        c2 <- (x >= 1/2 & x <= 1)  #2ed situation
        
        a1 <- (c1 + c2 * -1) * (-sqrt((1 - 2 * x)^2/(2 * (1 - x) * x)))
        
        a2 <- x
        a2[x != 0.5] <- (c1[x != 0.5] + c2[x != 0.5] * -1) * (1 - 2 * x[x != 
                                                                          0.5])/sqrt(8 * (1 - 2 * x[x != 0.5])^2 * ((1 - x[x != 0.5]) * x[x != 
                                                                                                                                            0.5])^3)
        
        a2[x == 0.5] <- 2 * sqrt(2)
        
        (exp(asinh((-mu + a1)/sigma)) * a2)/((1 + exp(asinh((-mu + a1)/sigma)))^2 * 
                                               sigma * sqrt(1 + (mu - a1)^2/sigma^2))
      }
      
      
    }
    
  }
  
   # burr7-XX-----
    if (fd =="burr7"){
    if (sd == "arcsinh") {
      dent <- function(x, mu, sigma) 
        ((1 + 2*(-1 + x)*x)*sech(((1 - 2*x)/(2*(-1 + x)*x) - mu)/sigma)^2)/(4*(-1 + x)^2*x^2*sigma)
    }
    if (sd == "burr7") {
      dent <- function(x, mu, sigma) sech((mu + atanh(1 - 2*x))/sigma)^2/(4*x*sigma - 4*x^2*sigma)
    } 
      
    if (sd == "burr8") {
      dent <- function(x, mu, sigma) (pi*csc((pi*x)/2)*sec((pi*x)/2)*sech((mu - log(tan((pi*x)/2)))/sigma)^2)/(4*sigma)
    }
      
    if (sd == "cauchy") {
      dent <- function(x, mu, sigma) (pi*csc(pi*x)^2*sech((mu + cot(pi*x))/sigma)^2)/(2*sigma)
    }   
    
    if (sd == "logistic") {
      dent <- function(x, mu, sigma) -(sech((-mu + log(1/(1 - x)) + log(x))/sigma)^2/(2*(-1 + x)* x*sigma))
    }
      
    if (sd == "t2") {
      dent <- function(x, mu, sigma) sech(((-1 + 2*x)/(sqrt(2)* sqrt((-(-1 + x))*x)) - mu)/sigma)^2/(4* sqrt(2)*((-(-1 + x))*x)^(3/2)*sigma)
      }
  }
  
  
  # burr8-XX-----
  if (fd == "burr8") {
    if (sd == "arcsinh") {
      dent <- function(x, mu, sigma) 
        (exp((1 + 2*x + 2*(-1 + x)*x*mu)/(2*(-1 + x)*x*sigma))*(1 + 2*(-1 + x)*x))/((exp(1/((-1 + x)*x*sigma)) + 
            exp((2*(1/(-1 + x) + mu))/sigma))*pi* (-1 + x)^2*x^2*sigma)
    } 
    
    if (sd == "burr7") {
      dent <- function(x, mu, sigma) sech((mu + atanh(1 - 2 * x))/sigma)/(2 * 
                                                                            pi * x * sigma - 2 * pi * x^2 * sigma)
    }
    
    if (sd == "burr8") {
      dent <- function(x, mu, sigma) (2 * exp(mu/sigma) * csc(pi * x) * tan((pi * 
                                                                               x)/2)^(1/sigma))/(exp((2 * mu)/sigma) * sigma + sigma * tan((pi * 
                                                                                                                                              x)/2)^(2/sigma))
    }
    
    if (sd == "cauchy") {
      dent <- function(x, mu, sigma) (2 * exp((mu + cot(pi * x))/sigma) * csc(pi * 
                                                                                x)^2)/((1 + exp((2 * (mu + cot(pi * x)))/sigma)) * sigma)
    }
    
    if (sd == "logistic") {
      dent <- function(x, mu, sigma) 
        (2*exp(mu/sigma)*(-(x/(-1 + x)))^(-1 + 1/sigma))/(pi*(-1 + x)^2*(exp((2*mu)/sigma) + 
            (-(x/(-1 + x)))^(2/sigma))*sigma)        
        } 
    
    if (sd == "t2") {
      
      dent <- function(x, mu, sigma) sech(((1 - 2 * x)/(sqrt(2) * sqrt((1 - x) * 
                                                                         x)) + mu)/sigma)/(2 * sqrt(2) * pi * ((1 - x) * x)^(3/2) * sigma)
    }
    
  }
  
  
      #Cauchit-XX-----
  if (fd == "cauchit") {
    # cauchit-arcsinh
    if (sd == "arcsinh") {
   dent <- function(x, mu, sigma) (2*(1 + 2*(-1 + x)*x)*sigma)/(pi*((-1 + 2*x + 2*(-1 + x)*x*mu)^2 + 4*(-1 + x)^2*x^2*sigma^2))
    }
    
    # cauchit-Cauchy
    if (sd == "cauchy") {
      dent <- function(x, mu, sigma) (sigma*csc(pi*x)^2)/(mu^2 + sigma^2 + 2*mu*cot(pi*x) +cot(pi*x)^2)
    }
    
    # cauchit-burr7
    if (sd == "burr7") {
   dent <- function(x, mu, sigma) -(sigma/(2*pi*(-1 + x)*x*(mu^2 + sigma^2 + 2*mu*atanh(1 - 2*x) + atanh(1 - 2*x)^2)))
     
     }
    
    # cauchit-burr8
    if (sd == "burr8") {
   dent <- function(x, mu, sigma) (sigma*csc(pi*x))/(mu^2 + sigma^2 - 2*mu*log(tan((pi*x)/2)) + log(tan((pi*x)/2))^2)
     
     }
    
     # cauchit-Logistic 
    if (sd == "logistic") {
    dent <- function(x, mu, sigma) -(sigma/(pi*(-1 + x)*x*(mu^2 + sigma^2 - 2*mu*log(x/(1 - x)) + log(x/(1 - x))^2)))
    }
    
     # cauchit-t2 
    if (sd == "t2") {
    dent <- function(x, mu, sigma) {
       c1 <- ((x >= 0) & (x < 1/2))  #1st situation
        c2 <- (x >= 1/2 & x <= 1)  #2ed situation
        
        a1 <- (c1 + c2 * -1) * (-sqrt((1 - 2 * x)^2/(2 * (1 - x) * x)))
        
        a2 <- x
        a2[x != 0.5] <- (c1[x != 0.5] + c2[x != 0.5] * -1) * (1 - 2 * x[x != 
          0.5])/sqrt(8 * (1 - 2 * x[x != 0.5])^2 * ((1 - x[x != 0.5]) * x[x != 
          0.5])^3)
        
        a2[x == 0.5] <- 2 * sqrt(2)
        
        a2/(pi*sigma*((mu - a1)^2/sigma^2 + 1))
    }
      
      
      }
    
  }
  
  
  # logit-XX-----
  if (fd == "logit") {
    
    if (sd == "arcsinh") {
      dent <- function(x, mu, sigma) 
        (exp((1 + 2*x + 2*x*mu - 2*x^2*mu)/(2*x*sigma - 2*x^2*sigma))*(1 - 2*x + 2*x^2))/ (2*(exp(1/(sigma - x*sigma)) + exp(mu/sigma + 1/(2*x*sigma - 2*x^2*sigma)))^2*(-1 + x)^2*x^2*sigma)
    }
    
    if (sd == "burr7") {
      dent <- function(x, mu, sigma) sech((mu + atanh(1 - 2*x))/(2*sigma))^2/(8*x*sigma -8*x^2*sigma)
    }
    
    if (sd == "burr8") {
      dent <- function(x, mu, sigma) (exp(mu/sigma) * pi * csc(pi * x) * tan((pi * 
        x)/2)^(1/sigma))/(sigma * (exp(mu/sigma) + tan((pi * x)/2)^(1/sigma))^2)
    }
    
    if (sd == "cauchy") {
      dent <- function(x, mu, sigma) (exp((mu + cot(pi * x))/sigma) * pi * csc(pi * 
        x)^2)/((1 + exp((mu + cot(pi * x))/sigma))^2 * sigma)
    }
    
    if (sd == "logistic") {
      dent <- function(x, mu, sigma) (exp(mu/sigma) * (-1 + 1/x)^(-1 + 1/sigma))/(((x + 
        exp(mu/sigma) * (-1 + 1/x)^(1/sigma) * x)^2) * sigma)
    }
    

    
    if (sd == "t2") {
      dent <- function(x, mu, sigma) {
        c1 <- ((x >= 0) & (x < 1/2))  #1st situation
        c2 <- (x >= 1/2 & x <= 1)  #2ed situation
        
        a1 <- (c1 + c2 * -1) * (-sqrt((1 - 2 * x)^2/(2 * (1 - x) * x)))
        
        a2 <- x
        a2[x != 0.5] <- (c1[x != 0.5] + c2[x != 0.5] * -1) * (1 - 2 * x[x != 
          0.5])/sqrt(8 * (1 - 2 * x[x != 0.5])^2 * ((1 - x[x != 0.5]) * x[x != 
          0.5])^3)
        
        a2[x == 0.5] <- 2 * sqrt(2)
        (exp((mu + a1)/sigma) * a2)/((exp(mu/sigma) + exp(a1/sigma))^2 * 
          sigma)
      }
    }
    

    
    
  }
  
 
  # t2-XX-----
  if (fd == "t2") {
            # t2-ArcSinh 
    if (sd == "arcsinh") {
      dent <- function(x, mu, sigma) (4*(1 - 2*x + 2*x^2))/((-1 + x)^2* x^2*sigma*((1 + 4*x*(-1 + mu) + 
          4*x^4*(mu^2 + 2*sigma^2) + 4*x^2*(1 - 3*mu + mu^2 + 2*sigma^2) - 
          8*x^3*(-mu + mu^2 + 2*sigma^2))/((-1 + x)^2*x^2*sigma^2))^(3/2))
    }
    
    
    if (sd == "cauchy") {
      dent <- function(x, mu, sigma) (pi * sigma * csc(pi * x)^2)/((mu^2 + 2 * 
        sigma^2 + 2 * mu * cot(pi * x) + cot(pi * x)^2) * sqrt(2 + (mu + 
        cot(pi * x))^2/sigma^2))
    }
    
    
    if (sd == "t2") {
      dent <- function(x, mu, sigma) 1/(sigma * ((1 - 4 * sqrt(2 - 2 * x) * x^(3/2) * 
        mu + 2 * sqrt(2) * sqrt((1 - x) * x) * mu + 2 * x * (-2 + mu^2 + 
        2 * sigma^2) - 2 * x^2 * (-2 + mu^2 + 2 * sigma^2))/sigma^2)^(3/2))
    }
    
    if (sd == "burr7") {
      dent <- function(x, mu, sigma) -(1/(2 * (-1 + x) * x * sigma * ((mu^2 + 
        2 * sigma^2 + 2 * mu * atanh(1 - 2 * x) + atanh(1 - 2 * x)^2)/sigma^2)^(3/2)))
    }
    
    if (sd == "burr8") {
      dent <- function(x, mu, sigma) (pi * csc(pi * x))/(sigma * ((mu^2 + 2 * 
        sigma^2 - 2 * mu * log(tan((pi * x)/2)) + log(tan((pi * x)/2))^2)/sigma^2)^(3/2))
    }
    
    if (sd == "logistic") {
      dent <- function(x, mu, sigma) -(sigma/((-1 + x)*x*sqrt(2 + (mu - log(x/(1 - x)))^2/sigma^2)*(mu^2 + 2*sigma^2 - 2*mu*log(x/(1 - x)) + log(x/(1 - x))^2)))
    }
    
  }
  
 
  
  if (fd == "km" | sd == "km") {
    dent <- function(x, mu, sigma) 1 - (1 - x^mu)^sigma
  }
  
  # Return output
  dent(x, mu, sigma)
}



