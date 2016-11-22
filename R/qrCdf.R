#' @rdname dq
#' @export
pq <- function(q, mu, sigma, fd, sd) {
  
  fd <- tolower(fd)
  sd <- tolower(sd)
  
  if (any(q <= 0) | any(q >= 1)) 
    stop(paste("q must be between 0 and 1", "\n", ""))
  
  # arcsinh-XX----
  if (fd == "arcsinh") {
    if (sd == "arcsinh") {
      prt <- function(q, mu, sigma) 1/(1 + exp(-asinh(((1 - 2*q)/(2*(-1 + q)*q) - mu)/sigma)))
    }
    
    if (sd == "burr7") {
      prt <- function(q, mu, sigma) 1/(1 + exp(-asinh((-atanh(1 - 2*q) - mu)/sigma)))
    }
    
    if (sd == "burr8") {
      prt <- function(q, mu, sigma) 1/(1 + exp(asinh((mu - log(tan(pi * q/2)))/sigma)))
    } 
    
    if (sd == "cauchy") {
      prt <- function(q, mu, sigma) 1/(1 + exp(-asinh((-cot(pi * q) - mu)/sigma)))
    }
    
    if (sd == "t2") {
      prt <- function(q, mu, sigma) {
        w <- function(q) {
          x <- q
          for (i in 1:length(q)) {
            c1 <- ((q[i] >= 0) & (q[i] < 1/2))  #1st situation
            c2 <- (q[i] >= 1/2 & q[i] <= 1)  #2ed situation
            
            if (c1) {
              x[i] <- -sqrt((1 - 2 * q[i])^2/(2 * (1 - q[i]) * q[i]))
            } else if (c2) {
              x[i] <- sqrt((1 - 2 * q[i])^2/(2 * (1 - q[i]) * q[i]))
            }
          }
          x
        }
        1/(1 + exp(asinh((mu - w(q))/sigma)))
      }
    }
    

    
    if (sd == "logistic") {
      prt <- function(q, mu, sigma) 1/(1 + exp(-asinh((log(q/(1 - q)) - mu)/sigma)))
    }
  }

  # burr7-XX----
  if (fd =="burr7"){
      
    if (sd == "arcsinh") {
      prt <- function(q, mu, sigma)  (tanh(((1 - 2*q)/(2*(-1 + q)*q) - mu)/sigma) + 1)/2
    }
    
    if (sd == "burr7") {
      prt <- function(q, mu, sigma)  (tanh((-atanh(1 - 2*q) - mu)/sigma) + 1)/2
    }    
    if (sd == "burr8") {
      prt <- function(q, mu, sigma)  (tanh((log(tan((pi*q)/2)) - mu)/sigma) + 1)/2
    }
    
    if (sd == "cauchy") {
      prt <- function(q, mu, sigma)  (tanh((tan((2*pi*q - pi)/2) - mu)/sigma) + 1)/2
    } 
    
    if (sd == "logistic") {
      prt <- function(q, mu, sigma)  (tanh((log(q/(1 - q)) - mu)/sigma) + 1)/2
    }
    
    if (sd == "t2") {
      prt <- function(q, mu, sigma) {
        w <- function(q) {
          x <- q
          for (i in 1:length(q)) {
            c1 <- ((q[i] >= 0) & (q[i] < 1/2))  #1st situation
            c2 <- (q[i] >= 1/2 & q[i] <= 1)  #2ed situation
            
            if (c1) {
              x[i] <- -sqrt((1 - 2 * q[i])^2/(2 * (1 - q[i]) * q[i]))
            } else if (c2) {
              x[i] <- sqrt((1 - 2 * q[i])^2/(2 * (1 - q[i]) * q[i]))
            }
          }
          x
        }
        
        (tanh((w(q) - mu)/sigma) + 1)/2
      }
    } 
    
      }
  
  # burr8-XX----
  if (fd == "burr8") {
    if (sd == "arcsinh") {
      prt <- function(q, mu, sigma) 2*(atan(exp(((1 - 2*q)/(2*(-1 + q)*q) - mu)/sigma))/pi)
    } 
    
    if (sd == "burr7") {
      prt <- function(q, mu, sigma) 2 * (atan(exp((-atanh(1 - 2 * q) - mu)/sigma))/pi)
    }
    
    if (sd == "burr8") {
      prt <- function(q, mu, sigma) 2 * (atan(exp((log(tan((pi * q)/2)) - mu)/sigma))/pi)
    }
    
    if (sd == "cauchy") {
      prt <- function(q, mu, sigma) 2 * (atan(exp((tan((2 * pi * q - pi)/2) - 
                                                     mu)/sigma))/pi)
    }
    
    if (sd == "logistic") {
      prt <- function(q, mu, sigma) 2*(atan(exp((log(-(q/(-1 + q))) - mu)/sigma))/pi)
    }
    
    
    if (sd == "t2") {
      
      prt <- function(q, mu, sigma) {
        w <- function(q) {
          x <- q
          for (i in 1:length(q)) {
            c1 <- ((q[i] >= 0) & (q[i] < 1/2))  #1st situation
            c2 <- (q[i] >= 1/2 & q[i] <= 1)  #2ed situation
            
            if (c1) {
              x[i] <- -sqrt((1 - 2 * q[i])^2/(2 * (1 - q[i]) * q[i]))
            } else if (c2) {
              x[i] <- sqrt((1 - 2 * q[i])^2/(2 * (1 - q[i]) * q[i]))
            }
          }
          x
        }
        
        2 * (atan(exp((w(q) - mu)/sigma))/pi)
      }
    }
    
  }


  #Cauchit-XX-----
  if (fd == "cauchit") {
    # cauchit-arcsinh
    if (sd == "arcsinh") {
   prt <- function(q, mu, sigma) 1/2 + atan(((1 - 2*q)/(2*(-1 + q)*q) - mu)/sigma)/pi
    }
    
    # cauchit-Cauchy
    if (sd == "cauchy") {
     prt <- function(q, mu, sigma) 1/2 + atan((tan((2*pi*q - pi)/2) - mu)/sigma)/pi
    }
    
    # cauchit-burr7
    if (sd == "burr7") {
   prt <- function(q, mu, sigma) 1/2 + atan((-atanh(1 - 2*q) - mu)/sigma)/pi
   }
  
    # cauchit-burr8
    if (sd == "burr8") {
   prt <- function(q, mu, sigma) 1/2 + atan((log(tan((pi*q)/2)) - mu)/sigma)/pi
   }
      
     # cauchit-Logistic 
    if (sd == "logistic") {
    prt <- function(q, mu, sigma) 1/2 + atan((log(q/(1 - q)) - mu)/sigma)/pi
    }
    
     # cauchit-t2 
    if (sd == "t2") {
      prt <- function(q, mu, sigma) {
        w <- function(q) {
          x <- q
          for (i in 1:length(q)) {
          c1 <- ((q[i] >= 0) & (q[i] < 1/2))  #1st situation
          c2 <- (q[i] >= 1/2 & q[i] <= 1)  #2ed situation
          
          if (c1) {
            x[i] <- -sqrt((1 - 2 * q[i])^2/(2 * (1 - q[i]) * q[i]))
          } else if (c2) {
            x[i] <- sqrt((1 - 2 * q[i])^2/(2 * (1 - q[i]) * q[i]))
          }
          }
          x
        }
        1/2 + atan((w(q) - mu)/sigma)/pi
      }
      
    }
    
  }
  
  
  
  # logit-XX----
  if (fd == "logit") {
    
    if (sd == "arcsinh") {
      prt <- function(q, mu, sigma) 1/(1 + exp(-((1 - 2*q)/(2*(-1 + q)*q) - mu)/sigma))
    }  
    
    if (sd == "burr7") {
      prt <- function(q, mu, sigma) 1/(1 + exp(-(-atanh(1 - 2*q) - mu)/sigma))
    }
    
    if (sd == "burr8") {
      prt <- function(q, mu, sigma) 1/(1 + exp((mu - log(tan((pi * q)/2)))/sigma))
    }
    
    if (sd == "cauchy") {
      prt <- function(q, mu, sigma) 1/(1 + exp((mu + cot(pi * q))/sigma))
    } 
    
    
    if (sd == "logistic") {
      prt <- function(q, mu, sigma) 1/(1 + exp((mu + log(-1 + 1/q))/sigma))
    }
    
    
    if (sd == "t2") {
      prt <- function(q, mu, sigma) {
        w <- function(q) {
          x <- q
          for (i in 1:length(q)) {
          c1 <- ((q[i] >= 0) & (q[i] < 1/2))  #1st situation
          c2 <- (q[i] >= 1/2 & q[i] <= 1)  #2ed situation
          
          if (c1) {
            x[i] <- -sqrt((1 - 2 * q[i])^2/(2 * (1 - q[i]) * q[i]))
          } else if (c2) {
            x[i] <- sqrt((1 - 2 * q[i])^2/(2 * (1 - q[i]) * q[i]))
          }
          }
          x
        }
        1/(1 + exp((mu - w(q))/sigma))
      }
      
    }
    
  }
  
  
  # t2-XX----
  if (fd == "t2") {
     if (sd == "arcsinh") {
      prt <- function(q, mu, sigma) {
       1/2 + ((1 - 2*q)/(2*(-1 + q)*q) - mu)/(sigma*(2*sqrt(2 + (((1 - 2*q)/(2*(-1 + q)*q) - mu)/sigma)^2)))
        }
    }
    
    if (sd == "cauchy") {
      prt <- function(q, mu, sigma) (1/2) * (1 - (mu + cot(pi * q))/(sigma * 
        sqrt(2 + (mu + cot(pi * q))^2/sigma^2)))
    }
    
    
    if (sd == "t2") {
      prt <- function(q, mu, sigma) {
        w <- function(q) {
          x <- q
          for (i in 1:length(q)) {
          c1 <- ((q[i] >= 0) & (q[i] < 1/2))  #1st situation
          c2 <- (q[i] >= 1/2 & q[i] <= 1)  #2ed situation
          
          if (c1) {
            x[i] <- -sqrt((1 - 2 * q[i])^2/(2 * (1 - q[i]) * q[i]))
          } else if (c2) {
            x[i] <- sqrt((1 - 2 * q[i])^2/(2 * (1 - q[i]) * q[i]))
          }
          }
          x
        }
        
        (1/2) * (1 + (-mu + w(q))/(sigma * sqrt(2 + (mu - w(q))^2/sigma^2)))
      }
    }
    
    if (sd == "burr7") {
      prt <- function(q, mu, sigma) {
        w <- function(q) {
          x <- q
          for (i in 1:length(q)) {
          c1 <- ((q[i] >= 0) & (q[i] < 1/2))  #1st situation
          c2 <- (q[i] >= 1/2 & q[i] <= 1)  #2ed situation
          
          if (c1) {
            x[i] <- -sqrt((1 - 2 * q[i])^2/(2 * (1 - q[i]) * q[i]))
          } else if (c2) {
            x[i] <- sqrt((1 - 2 * q[i])^2/(2 * (1 - q[i]) * q[i]))
          }
          }
          x
        }
        
        (1/2) * (1 - (mu + atanh(1 - 2 * q))/(sigma * sqrt(2 + (mu + atanh(1 - 
          2 * q))^2/sigma^2)))
      }
      
    }
    
    if (sd == "burr8") {
      prt <- function(q, mu, sigma) (1/2) * (1 + (-mu + log(tan((pi * q)/2)))/(sigma * 
        sqrt(2 + (mu - log(tan((pi * q)/2)))^2/sigma^2)))
      
    }
    
    if (sd == "logistic") {
      prt <- function(q, mu, sigma) 1/2 + (log(q/(1 - q)) - mu)/(sigma*(2*sqrt(2 + ((log(q/(1 - q)) - mu)/sigma)^2)))
    }
    
  }
  

  
  if (fd == "km" | sd == "km") {
    prt <- function(q, mu, sigma) 1 - (1 - q^mu)^sigma
  }
  
  prt(q, mu, sigma)
}

 