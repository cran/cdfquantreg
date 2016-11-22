#' @rdname dq
#' @export
rq <- function(n, mu, sigma, fd, sd) {
  nq <- runif(n)
  samp <- qq(nq, mu, sigma, fd, sd)
  samp
}

#' @rdname dq
#' @export
qq <- function(p, mu, sigma, fd, sd) {
  
  fd <- tolower(fd)
  sd <- tolower(sd)
  
  if (any(p <= 0) | any(p >= 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  # arcsinh-XX----
  if (fd == "arcsinh") {
    
    if (sd == "arcsinh") {
      qant <- function(p, mu, sigma) 
        1/(1 + exp(-asinh(((1 - 2*p)/(2*(-1 + p)*p) + mu/sigma)/(1/sigma))))
    }
    
    if (sd == "burr7") {
      qant <- function(p, mu, sigma) 
         1/2* (1 + tanh(mu + (sigma - 2 * p *sigma)/( 2 *(-1 + p)* p)))
    }
    
    if (sd == "burr8") {
      qant <- function(p, mu, sigma) 
        2 * (atan(exp(mu + (sigma - 2 * p * sigma)/(2 * (-1 + p) * p)))/(pi))
    }
    
    if (sd == "cauchy") {
      qant <- function(p, mu, sigma) 
        1/2 + atan(mu + (sigma - 2 * p * sigma)/(2 * (-1 + p) * p))/(pi)
    }
    
    if (sd == "logistic") {
     qant <- function(p, mu, sigma)  
       #exp(mu + sigma*Re((0+1i)*sin((-(0+1i))*log(-(p/(p - 1))))))/(1 + exp(mu + sigma*Re((0+1i)*sin((-(0+1i))*log(-(p/(p - 1)))))))
        exp(mu)/(exp(mu) + exp((sigma - 2*p*sigma)/(2*p - 2*(p^2))))
     }
   
    if (sd == "t2") {
      qant <- function(p, mu, sigma) 1/2 + (((1 - 2 * p)/(2 * (-1 + p) * p) + 
                                               mu/sigma) * sigma)/(2 * sqrt(2 + (((1 - 2 * p)/(2 * (-1 + p) * p) + 
                                                                                    mu/sigma) * sigma)^2))
      
      
    }
  }
   
  # burr7-XX----
  if (fd =="burr7"){
    if (sd == "arcsinh") {
       qant <- function(p, mu, sigma) 
          (-1 + mu - sigma *atanh(1 - 2* p) + sqrt(1 + mu^2 -2* mu *sigma* atanh(1 - 2* p) + sigma^2 *atanh( 1 - 2* p)^2))/( 2 *mu - 2 *sigma* atanh(1 - 2* p))

      }    
  
    if (sd == "burr7") {
       qant <- function(p, mu, sigma) (1/2 )*(1 + tanh(mu - sigma *atanh(1 - 2* p)))
    }
    
    if (sd == "burr8") {
       qant <- function(p, mu, sigma) (2* atan(exp(mu - sigma* atanh(1 - 2* p))))/(pi)
    }
    
    if (sd == "cauchy") {
       qant <- function(p, mu, sigma) 1/2 + atan(mu - sigma *atanh(1 - 2 *p))/(pi)
    }
    
    if (sd == "logistic") {
       qant <- function(p, mu, sigma) exp(mu)/(exp(mu) + exp(sigma* atanh(1 - 2*p)))
    }
    
    if (sd == "t2") {
       qant <- function(p, mu, sigma)  1/2 + ((-atanh(1 - 2* p) + mu/sigma)*sigma)/( 2* sqrt(2 + ((-atanh(1 - 2*p) + mu/sigma)*sigma)^2))
      }
  }
  
  # burr8-XX----
  if (fd == "burr8") {
    if (sd == "arcsinh") {
       qant <- function(p, mu, sigma) 
         1/(1 - mu + sigma* log(cot(((pi)*p)/2)) + sqrt( 1 + mu^2 - 2*mu * sigma*(log(cot(((pi)*p)/2))) + sigma^2 * (log( cot(((pi)*p)/2)))^2))

      }  
    
    if (sd == "cauchy") {
      qant <- function(p, mu, sigma) 1/2 + atan(mu - sigma * (log(cot((pi * p)/2))))/pi
    }
    
    
    if (sd == "logistic") {
       qant <- function(p, mu, sigma) 
         exp(mu)/(exp(mu)+ cot(((pi)* p)/2)^sigma)
      }  
    
    if (sd == "t2") {
      
      qant <- function(p, mu, sigma) 1/2 + (mu + sigma * log(tan(((pi) * p)/2)))/(2 * 
                                                                                    sqrt(2 + (mu + sigma * log(tan(((pi) * p)/2)))^2))
    }
    
    
    if (sd == "burr7") {
      qant <- function(p, mu, sigma) (1/2) * (1 + tanh(mu - sigma * (log(cot(pi * 
                                                                               p/2)))))
    }
    
    if (sd == "burr8") {
      qant <- function(p, mu, sigma) 2 * atan(exp(mu - sigma * log(cot(pi * p/2))))/pi
    }
    
  }
  
        #Cauchit-XX-----
  if (fd == "cauchit") {
    # cauchit-arcsinh
    if (sd == "arcsinh") {
  qant <- function(p, mu, sigma)   (-1 + mu - sigma* cot((pi)* p) + sqrt( 1 + mu^2 - 2* mu* sigma* cot((pi)* p) + sigma^2* cot((pi)* p)^2))/(2* mu - 2* sigma* cot((pi)* p))

    }
    
    # cauchit-burr7
    if (sd == "burr7") {
      qant <- function(p, mu, sigma) (1/2) *(1 + tanh(mu - sigma* cot((pi)* p)))
    }
    
    # cauchit-burr8
    if (sd == "burr8") {
      qant <- function(p, mu, sigma) (2* (atan(exp(mu - sigma* cot((pi)* p)))))/(pi)
    }
    
    # cauchit-Cauchy
    if (sd == "cauchy") {
      qant <- function(p, mu, sigma) 1/2 + atan((-(1/tan(pi*p)) + mu/sigma)*sigma)/pi
    }
    
     # cauchit-Logistic 
    if (sd == "logistic") {
    qant <- function(p, mu, sigma) 1/(1 + exp(-mu + sigma*tan(0.5*pi - pi*p)))
    }
    
    if (sd == "t2") {
    qant <- function(p, mu, sigma) {
     x <- (-cot(pi*p) + (mu/sigma))*sigma
     1/2 + x/(2*sqrt(2 + x^2)) 
    }
    }
  }
  
  
  
  
  # logit-XX----
  if (fd == "logit") {
    if (sd == "arcsinh") {
      qant <- function(p, mu, sigma) 1/(1 + exp(-asinh((-log(-1 + 1/p) + mu/sigma)*sigma)))
    }
    
    if (sd == "logistic") {
      qant <- function(p, mu, sigma) 1/(1 + exp((-(-log(-1 + 1/p) + mu/sigma)) * 
        sigma))
    }
    
    if (sd == "cauchy") {
      qant <- function(p, mu, sigma) 1/2 + atan((-log(-1 + 1/p) + mu/sigma) * 
        sigma)/pi
    }
    
    if (sd == "t2") {
      qant <- function(p, mu, sigma) (1/2 )*(1 + tanh(mu - sigma* (log(-1 + 1/p))))
      
    }
    
   if (sd == "burr7") {
      qant <- function(p, mu, sigma) (1/2 )*(1 + tanh(mu - sigma* (log(-1 + 1/p)))) 
   }
    
    if (sd == "burr8") {
      qant <- function(p, mu, sigma) (2 * (atan(exp(mu - sigma * (log(-1 + 1/p))))))/(pi)
    }
  }
  

  # t2-XX----
  if (fd == "t2") {
    
    if (sd == "arcsinh") {
       qant <- function(p, mu, sigma) {
        w <- function(p) {
         x <- p
    for (i in 1:length(p))
      {
      c1 <- ((p[i] >= 0) & (p[i] < 1/2))  #1st situation
      c2 <- (p[i] >= 1/2 & p[i] <= 1)  #2ed situation
      
      if (c1) { x[i]= -sqrt((1 - 2 * p[i])^2/(2 * (1 - p[i]) *p[i]))
      } else if (c2) {x[i]= sqrt((1 - 2 *p[i])^2/(2 * (1 - p[i]) * p[i]))
        }
      }
      x
        }
  
     Re(1/(1 + exp(-asinh(mu + sigma*w(p)))))
       }
    }
    
    if (sd == "cauchy") {
      qant <- function(p, mu, sigma) {
        w <- function(p) {
          x <- p
          for (i in 1:length(p)) {
          c1 <- ((p[i] >= 0) & (p[i] < 1/2))  #1st situation
          c2 <- (p[i] >= 1/2 & p[i] <= 1)  #2ed situation
          
          if (c1) {
            x[i] <- -sqrt((1 - 2 * p[i])^2/(2 * (1 - p[i]) * p[i]))
          } else if (c2) {
            x[i] <- sqrt((1 - 2 * p[i])^2/(2 * (1 - p[i]) * p[i]))
          }
          }
          x
        }
        
        Re(1/2 + atan((w(p) + mu/sigma) * sigma)/pi)
      }
    }
    
    
    if (sd == "t2") {
      qant <- function(p, mu, sigma) {
        w <- function(p) {
          x <- p
          for (i in 1:length(p)) {
          c1 <- ((p[i] >= 0) & (p[i] < 1/2))  #1st situation
          c2 <- ((p[i] >= 1/2) & (p[i] <= 1))  #2ed situation
          
          if (c1) {
            x[i] <- -sqrt((1 - 2 * p[i])^2/(2 * (1 - p[i]) * p[i]))
          } else if (c2) {
            x[i] <- sqrt((1 - 2 * p[i])^2/(2 * (1 - p[i]) * p[i]))
          }
          }
          x
        }
        
        1/2 + (w(p) * sigma + mu)/(2 * sqrt(2 + (w(p) * sigma + mu)^2))
      }
    }
    
    if (sd == "burr7") {
      qant <- function(p, mu, sigma) {
        
        w <- function(p) {
          x <- p
          for (i in 1:length(p)) {
          c1 <- ((p[i] >= 0) & (p[i] < 1/2))  #1st situation
          c2 <- (p[i] >= 1/2 & p[i] <= 1)  #2ed situation
          
          if (c1) {
            x[i] <- -sqrt((1 - 2 * p[i])^2/(2 * (1 - p[i]) * p[i]))
          } else if (c2) {
            x[i] <- sqrt((1 - 2 * p[i])^2/(2 * (1 - p[i]) * p[i]))
          }
          }
          x
        }
        
        (tanh((w(p) + mu/sigma) * sigma) + 1)/2
      }
      
      
    }
    
    if (sd == "burr8") {
      qant <- function(p, mu, sigma) {
        w <- function(p) {
          x <- p
          for (i in 1:length(p)) {
          c1 <- ((p[i] >= 0) & (p[i] < 1/2))  #1st situation
          c2 <- (p[i] >= 1/2 & p[i] <= 1)  #2ed situation
          
          if (c1) {
            x[i] <- -sqrt((1 - 2 * p[i])^2/(2 * (1 - p[i]) * p[i]))
          } else if (c2) {
            x[i] <- sqrt((1 - 2 * p[i])^2/(2 * (1 - p[i]) * p[i]))
          }
          }
          x
        }
        
        2 * (atan(exp((w(p) + mu/sigma) * sigma))/pi)
      }
    }
    
    if (sd == "logistic") {
      qant <- function(p, mu, sigma) {
        w <- function(q) {   
          x = q
    for (i in 1:length(q))
    {
      c1 <- ((q[i] >= 0) & (q[i] < 1/2))  #1st situation
      c2 <- (q[i] >= 1/2 & q[i] <= 1)  #2ed situation
      
      if (c1) { x[i]= -sqrt((1 - 2 * q[i])^2/(2 * (1 - q[i]) *q[i]))
      } else if (c2) {x[i]= sqrt((1 - 2 *q[i])^2/(2 * (1 - q[i]) * q[i]))
      }
    }
    x
        }
        
        Re(1/(1 + exp(-mu - sigma*w(p))))
      }
    }  
    
  }
  
 
  
  if (fd == "km" | sd == "km") {
    qant <- function(p, mu, sigma) (1 - (1 - p)^(1/sigma))^(1/mu)
  }
  
  qant(p, mu, sigma)
}
