#' @title The Family of Finite-Tailed Distributions
#' @aliases cdfft pdfft qqft rqft
#'
#' @description Density function, distribution function, quantile function, and random generation of variates for a specified cdf-quantile distribution.
#'
#' @param y vector of quantiles.
#' @param p vector of probabilities.
#' @param n Number of random samples.
#' @param q vector of quantiles.
#' @param sigma vector of standard deviations.
#' @param theta vector of skewness.
#' @param mu vector of means if 3-parameter case is used.
#' @param fd A string that specifies the parent distribution. At the moment, only "arcsinh", "cauchit" and "t2" can be used. See details. 
#' @param sd A string that specifies the child distribution. At the moment, only "arcsinh", "cauchy" and "t2" can be used. See details. 
#' @param inner A logic value that indicates if the inner (\code{inner = TRUE}) case or outer (\code{inner = FALSE}) will be used.
#' @param version A string indicates that which version will be used. "V" is the tilt parameter function while "W" indicates the Jones Pewsey transformation. 
#' @export
#' @import pracma
#' @return \code{pdfft} gives the density, \code{rqft} generates random variate, \code{qqft} gives the quantile function, and \code{cdfft} gives the cumulative density of specified distribution.


cdfft <-
  function(q,
           sigma,
           theta,
           fd,
           sd,
           mu = NULL,
           inner = TRUE,
           version) {
    if (is.null(mu)) {
      listfun <- .funlist(fd = fd,
                          sd = sd,
                          version = version)
    } else{
      listfun <- .funlist(
        fd = fd,
        sd = sd,
        npar = 3,
        version = version
      )
    }
    
    Func <- listfun$Func
    Hinv <- listfun$Hinv
   
    if (is.null(mu)){
      U <- listfun$U2
    } else{
      U <- listfun$U3
    }
    W <- listfun$W
    
    qcdf <- c(rep(0, length(q)))
    qcdf[q == 1] <- 1
    qcdf[q == 0] <- 0
    
    if(length(sigma) < length(q)){
      if(length(sigma) == 1){
        sigma <- rep(sigma, length(q))
        theta <- rep(theta, length(q))
        if(!is.null(mu)){
          mu <- rep(mu, length(q))
        }
      } else{
        stop("Parameter (e.g., sigma, theta, mu) should be either a single value or a vector that has the same length as q")
      }
    }
    
    if (is.null(mu)) {
      if (version == "V") {
        if (inner) {
          qcdf[q > 0 & q < 1] <-
            Func(U(Hinv(W(q[q > 0 & q < 1], theta[q > 0 & q < 1])),
                   sigma[q > 0 & q < 1]))
          
        } else{
          qcdf[q > 0 & q < 1] <-
            W(Func(U(Hinv(q[q > 0 & q < 1]),
                     sigma[q > 0 & q < 1])),
              theta[q > 0 & q < 1])
        }
      } else if (version == "W") {
        if (inner) {
          qcdf[q > 0 & q < 1] <-
            Func(U(W(Hinv(q[q > 0 &
                              q < 1]), theta[q > 0 & q < 1]), sigma[q > 0 & q < 1]))
        } else{
          qcdf[q > 0 & q < 1] <-
            Func(W(U(Hinv(q[q > 0 &
                              q < 1]), sigma[q > 0 & q < 1]), theta[q > 0 & q < 1]))
        }
        
      }
      
    } else{
      if (version == "V") {
        if (inner) {
          qcdf[q > 0 & q < 1] <-
            Func(U(Hinv(W(q[q > 0 & q < 1],
                          theta[q > 0 & q < 1])),
                   mu[q > 0 & q < 1],
                   sigma[q > 0 & q < 1]))
          
          
        } else{
          qcdf[q > 0 & q < 1] <-
            W(Func(U(Hinv(q[q > 0 & q < 1]),
                     mu[q > 0 & q < 1],
                     sigma[q > 0 & q < 1])),
              theta[q > 0 & q < 1])
          
        }
        
      } else if (version == "W") {
        if (inner) {
          qcdf[q > 0 & q < 1] <-
            Func(U(W(Hinv(q[q > 0 & q < 1]),
                     theta[q > 0 & q < 1]),
                   mu[q > 0 & q < 1],
                   sigma[q > 0 & q < 1]))
          
        } else{
          qcdf[q > 0 & q < 1] <-
            Func(W(U(Hinv(q[q > 0 &
                              q < 1]), mu[q > 0 &
                                            q < 1], sigma[q > 0 & q < 1]), theta[q > 0 & q < 1]))
          
        }
        
      }
      
    }
    
    
    qcdf
  }

#' @rdname cdfft
#' @export
pdfft <-
  function(y,
           sigma,
           theta,
           fd,
           sd,
           mu = NULL,
           inner = TRUE,
           version) {
    limts <- function(fd, sd, sigma, theta, bound, inner, version) {
      if( bound == 0){
        limtemp <- exp(theta)
      } else if( bound == 1){
        limtemp <- 1 / (exp(theta))
      }
      if (fd == "t2" & sd == "t2" & version == "W"){
        if( bound == 0){
          limtemp <- exp(2*theta)
        } else if( bound == 1){
          limtemp <- 1 / (exp(2*theta))
        }
      }
     
  
     
      if (fd == "arcsinh") {
        lim <- switch(sd,
                      arcsinh = sigma * limtemp,
                      cauchy  = pi * sigma * limtemp / 2)
      } else if (fd == "cauchit") {
        lim <- switch(sd,
                      arcsinh = 2 * sigma * limtemp / pi,
                      cauchy  = sigma * limtemp)
      } else if (fd == "t2" & sd == "t2") {
        lim <- limtemp * sigma ^ 2
        
      }
      lim
    }
    
    if(length(sigma) < length(y)){
      if(length(sigma) == 1){
        sigma <- rep(sigma, length(y))
        theta <- rep(theta, length(y))
        if(!is.null(mu)){
          mu <- rep(mu, length(y))
        }
      } else{
        stop("Parameter (e.g., sigma, theta, mu) should be either a single value or a vector that has the same length as y")
      }
    }
    qpdf <- c(rep(0, length(y)))
    qpdf[y == 0] <- limts(fd, sd, sigma = sigma[y == 0], 
                          theta = theta[y == 0], bound = 0, 
                          inner = inner, version = version)
    qpdf[y == 1] <- limts(fd, sd, sigma = sigma[y == 1], 
                          theta = theta[y == 1], bound = 1, 
                          inner = inner, version = version)
    
    if (is.null(mu)) {
      qpdf[y > 0 & y < 1] <-
        (
          cdfft(
            y[y > 0 & y < 1] + 0.000000001,
            sigma = sigma[y > 0 &
                            y < 1],
            theta = theta[y > 0 & y < 1],
            fd = fd,
            sd = sd,
            inner = inner,
            version = version
          ) -
            cdfft(
              y[y > 0 & y < 1],
              sigma = sigma[y > 0 &
                              y < 1],
              theta = theta[y > 0 & y < 1],
              fd = fd,
              sd = sd,
              inner = inner,
              version = version
            )
        ) / 0.000000001
    } else{
      qpdf[y > 0 & y < 1] <-
        (
          cdfft(
            y[y > 0 & y < 1] + 0.000000001,
            sigma = sigma[y > 0 & y < 1],
            theta = theta[y > 0 & y < 1],
            fd = fd,
            sd = sd,
            mu = mu[y > 0 & y < 1],
            inner = inner,
            version = version
          ) -
            cdfft(
              y[y > 0 & y < 1],
              sigma = sigma[y > 0 & y < 1],
              theta = theta[y > 0 & y < 1],
              fd = fd,
              sd = sd,
              mu = mu[y > 0 & y < 1],
              inner = inner,
              version = version
            )
        ) / 0.000000001
      
    }
    
    
    
    
    qpdf
  }

#' @rdname cdfft
#' @export
qqft <-
  function(p,
           sigma,
           theta,
           fd,
           sd,
           mu = NULL,
           inner = TRUE,
           version) {
    if (is.null(mu)) {
      listfun <- .funlist(fd = fd,
                          sd = sd,
                          version = version)
    } else{
      listfun <- .funlist(
        fd = fd,
        sd = sd,
        npar = 3,
        version = version
      )
    }
    
    Finv <- listfun$Finv
    H <- listfun$H
    if (is.null(mu)){
      Uinv <- listfun$Uinv2
    } else{
      Uinv <- listfun$Uinv3
    }
    Winv <- listfun$Winv
    
    qinvcdf <- c(rep(0, length(p)))
    qinvcdf[p == 1] <- 1
    qinvcdf[p == 0] <- 0
    
    if(length(sigma) < length(p)){
      if(length(sigma) == 1){
        sigma <- rep(sigma, length(p))
        theta <- rep(theta, length(p))
        if(!is.null(mu)){
          mu <- rep(mu, length(p))
        }
      } else{
        stop("Parameter (e.g., sigma, theta, mu) should be either a single value or a vector that has the same length as p")
      }
    }
    if (is.null(mu)) {
      if (inner & version == "V") {
        qinvcdf[p > 0 & p < 1] <-
          Winv(H(Uinv(Finv(p[p > 0 & p < 1]),
                      sigma[p > 0 & p < 1])), theta[p > 0 & p < 1])
        
      } else if(!inner & version == "V"){
        qinvcdf[p > 0 & p < 1] <-
          H(Uinv(Finv(Winv(p[p > 0 & p < 1], theta[p > 0 & p < 1])), sigma[p > 0 &
                                                                             p < 1]))
      } else if(inner & version == "W"){
        qinvcdf[p > 0 & p < 1] <-
        H(Winv(Uinv(Finv(p[p > 0 & p < 1]),sigma[p > 0 & p < 1]),theta[p > 0 & p < 1]))
        
      }else if(!inner & version == "W"){
        qinvcdf[p > 0 & p < 1] <-
          H(Uinv(Winv(Finv(p[p > 0 & p < 1]),theta[p > 0 & p < 1]),sigma[p > 0 & p < 1]))
      }
      
      
    } else{
      if (inner & version == "V") {
        qinvcdf[p > 0 & p < 1] <-
          Winv(H(Uinv(Finv(p[p > 0 & p < 1]),
                      mu[p > 0 & p < 1],
                      sigma[p > 0 & p < 1])),
               theta[p > 0 & p < 1])
        
    
      } else if(!inner & version == "V") {
        qinvcdf[p > 0 & p < 1] <-
          H(Uinv(Finv(Winv(p[p > 0 & p < 1],
                           theta[p > 0 & p < 1])),
                 mu[p > 0 & p < 1],
                 sigma[p > 0 & p < 1]))

      } else if(inner & version == "W") {
        qinvcdf[p > 0 & p < 1] <-
          H(Winv(Uinv(Finv(p[p > 0 & p < 1]),mu[p > 0 & p < 1],sigma[p > 0 & p < 1]),theta[p > 0 & p < 1]))
        
      } else if(!inner & version == "W") {
        qinvcdf[p > 0 & p < 1] <-
          H(Uinv(Winv(Finv(p[p > 0 & p < 1]),theta[p > 0 & p < 1]),mu[p > 0 & p < 1],sigma[p > 0 & p < 1]))
         
      }
      
      
      
      
    }
    qinvcdf
  }

#' @rdname cdfft
#' @export
rqft <-
  function(n,
           sigma,
           theta,
           fd,
           sd,
           mu = NULL,
           inner = TRUE,
           version) {
    nq <- runif(n)
    samp <- NULL
    if (length(theta)==1){
      for (i in 1:length(nq)){
        val <- qqft(nq[i], sigma, theta, fd, sd, mu, inner, version)
        samp <- c(samp, val)
      }
    }else if(length(theta)==nq & nq!=1){
      samp <- qqft(nq, sigma, theta, fd, sd, mu, inner, version)
    }
 
    samp
   
  }

.funlist <- function(fd, sd, npar = 2, version) {
  Func <- switch(
    fd,
    arcsinh = function(x) {
      1 / (exp(-asinh(x)) + 1)
    },
    cauchit = function(x) {
      (atan(x) / pi) + 0.5
    },
    t2 = function(x) {
      (x / (2 * sqrt((x) ^ 2 + 2))) + 0.5
    }
  )
  
  Finv <- switch(
    fd,
    arcsinh = function(x) {
      (1 - 2 * x) / (2 * x * (x - 1))
    },
    cauchit = function(x) {
      tan(pi * x - pi / 2)
    },
    t2 = function(x) {
      sign(x - 0.5) * (sqrt((1 - 2 * x) ^ 2)) / (sqrt(2) * sqrt((1 - x) * x))
    }
  )
  
  
  H <- switch(
    sd,
    arcsinh = function(x) {
      1 / (exp(-asinh(x)) + 1)
    },
    cauchy = function(x) {
      (atan(x) / pi) + 0.5
    },
    t2 = function(x) {
      (x / (2 * sqrt((x) ^ 2 + 2))) + 0.5
    }
  )
  
  
  Hinv <- switch(
    sd,
    arcsinh = function(x) {
      (1 - 2 * x) / (2 * x * (x - 1))
    },
    cauchy = function(x) {
      tan(pi * x - pi / 2)
    },
    t2 = function(x) {
      sign(x - 0.5) * (sqrt((1 - 2 * x) ^ 2)) / (sqrt(2) * sqrt((1 - x) * x))
    }
  )

    U2 <-  function(x, sigma) {
      x / sigma
    }
    Uinv2 <-  function(x, sigma) {
      sigma * x
    }
  
    U3 <-  function(x, mu, sigma) {
      (x - mu) / sigma
    }
    Uinv3 <-  function(x, mu, sigma) {
      sigma * x + mu
    }
  

  # if (npar == 2) {
  #   U <-  function(x, sigma) {
  #     x / sigma
  #   }
  #   Uinv <-  function(x, sigma) {
  #     sigma * x
  #   }
  # } else if (npar == 3) {
  #   U <-  function(x, mu, sigma) {
  #     (x - mu) / sigma
  #   }
  #   Uinv <-  function(x, mu, sigma) {
  #     sigma * x + mu
  #   }
  # }
  # 
  
  W <-  switch(
    version,
    V = function(x, theta) {
      (exp(theta) * x) / (1 - x + exp(theta) * x)
    },
    W = function(x, theta) {
      sinh(asinh(x) + theta)
    }
  )
  
  Winv <- switch(
    version,
    V = function(x, theta) {
      x / (x + exp(theta) - x * exp(theta))
    },
    W = function(x, theta) {
      sinh(asinh(x) - theta)
    }
  )
  
  list(
    Func = Func,
    Finv = Finv,
    H = H,
    Hinv = Hinv,
    U2 = U2,
    U3 = U3,
    W = W,
    Uinv2 = Uinv2,
    Uinv3 = Uinv3,
    Winv = Winv
  )
  # list(
  #   Func = Func,
  #   Finv = Finv,
  #   H = H,
  #   Hinv = Hinv,
  #   U = U,
  #   W = W,
  #   Uinv = Uinv,
  #   Winv = Winv
  # )
}

.qrLogLikFunFT <-
  function(h, y, sdata, tdata, fd, sd, mdata = NULL, inner = FALSE, version = "V") {
    # Parameters for estimating  sigma
    hx <- sdata %*% h[1:ncol(sdata)]
    sigma <- exp(hx)
    
    # Parameters for estimating theta
    gz <- tdata %*% h[ncol(sdata) + 1:ncol(tdata)]
    theta <- gz 
    
    if (!is.null(mdata)) {
      wk <- mdata %*% h[ncol(sdata) + ncol(tdata) + 1:ncol(mdata)]
      mu <- wk
    } else{
      mu <- NULL
    }
    
    # And get the log-Likelihood function for this value
    loglikv <- log(
      pdfft(
        y = y,
        sigma = sigma,
        theta = theta,
        fd = fd,
        sd = sd,
        mu = mu,
        inner = inner,
        version = version
      )
    )
    
    - sum(loglikv, na.rm = TRUE)
    
  }
