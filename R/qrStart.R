#' @title Starting Value Generation for CDF quantile Regressions
#' @aliases qrStart
#' @description \code{qrStart} is the function for generating starting values for a cdf-quantile GLM null model. 
#'
#' @param ydata The variable to be modeled
#' 
#' @param fd A string that specifies the parent distribution.
#' @param sd A string that specifies the sub-family distribution.
#' @param skew If ture, the starting values will be generated for the finited tailed distribution case.
#' @return A vector that consists initial values for mu and sigma.
#' @export
#' @details  The start values for the location parameter in a null model are 
#' the median of the empirical distribution, and a starting value for the dispersion parameter based on a 
#' specific quantile of the empirical distribution, specified according to the theoretical distribution
#' on which the model is based.
#' The start values for all new predictor coefficients in both the location and dispersion 
#' submodels are assigned the value 0.1. 
#' 
#' @examples
#' x <- rbeta(100, 1, 2)
#' qrStart(x, fd='t2', sd='t2')
#' #[1] -0.5938286  1.3996999

qrStart <- function(ydata,fd=NULL,sd=NULL, skew=FALSE) {
  ydata <- scaleTR(ydata)
  mu_0 = 0.1
  sigma_0 = 0.1

  m <- median(ydata) # Get the mean
  
  #The quantile used to generate qt depending on the parent distribution with some exceptions
  if (fd == "arcsinh") qt = quantile(ydata, 0.75);
  if (fd == "burr7")   qt = quantile(ydata, (1 - tanh(1))/2);
  if (fd == "burr8") qt = quantile(ydata,0.2244);
  if (fd == "cauchit") qt = quantile(ydata, 0.25);
  if (fd == "logit") qt = quantile(ydata, (1/(exp(1) + 1)))

  if (fd == "t2") qt = quantile(ydata, (3 - sqrt(3))/6)
  
  # The mu and sigma depending on the sub-family distribution some exceptions
if (sd == "arcsinh"){
  mu_0 <- (1 - 2 * m)/(2 * (-1 + m) * m)
  sigma_0 <- (1 - 2 * m)/(2 * (-1 + m) * m) - (1 - 2 * qt)/(2 * (-1 + qt) * qt) 
}

if (sd == "burr7"){
    mu_0 <- -atanh(1 - 2 * m)
    sigma_0 <- atanh(1 - 2 * qt) - atanh(1 - 2 * m) 
  }
  
  if (sd == "burr8")
  {
    mu_0 <- log(tan((m * pi)/2))
    sigma_0 <-  log(tan((m * pi)/2)) - log(tan((qt * pi)/2))
  }
  
  
  if (sd == "cauchy")
  {
    mu_0 <- tan((1/2)*(-1 + 2*m)*pi)
    sigma_0 <-  log(tan((m * pi)/2)) - log(tan((qt * pi)/2))
  }
  
  
  if (sd == "logit")
  {
    mu_0 <- log(-(m/(-1 + m)))
    sigma_0 <-  log((qt/(1 - qt))) - log(m/(1 - m))
  }
  
  if (sd == "t2")
  {
    mu_0 <- -sqrt((-1 + 4*m - 4*m^2)/(2*(-m + m^2)))
    sigma_0 <-  ((-sqrt(2))*sqrt(qt - 5*qt^2 + 8*qt^3 - 4*qt^4) - 2*qt*mu_0 +2*qt^2*mu_0)/(2*(-qt + qt^2))
    if (fd =="burr8" ) sigma_0 <- ((-sqrt(2))*sqrt(qt - 5*qt^2 + 8*qt^3 - 4*qt^4) + 2*qt*mu_0 -2*qt^2*mu_0)/(2*(-qt + qt^2))
  }
  
  if (sd == "weibull")
  {
    mu_0 <- log(-log(-m + 1))/2
    qt <- quantile(ydata, 1/(1 + exp(1)))
    sigma_0 <- mu_0 - log(-log(-qt + 1))/2
    
  }
  


  # Special cases******************
if (sd == "arcsinh" & fd == "t2")
{
  qt = quantile(ydata,(2 - sqrt(2))/2)
}

if (sd == "logit" & fd == "logit")
{
  mu_0 <- -log((m/(1 - m)))
  q25 <- quantile(ydata, 0.25)
  q75 <- quantile(ydata, 0.75)
  O25 <- q25/(1 - q25)
  O75 <- q75/(1 - q75)
  sigma_0 <- log(O75/O25)/(2 * log(3))
}


if (sd == "burr8" & fd == "arcsinh")
{
  mu_0 <- log(tan((m * pi)/2))
  qt <- quantile(ydata, 1/sqrt(2))
  sigma_0 <- log(tan((pi * qt)/2))-log(tan((m * pi)/2))
}

if (sd == "cauchy" & fd == "arcsinh")
{
  mu_0 <- tan((1/2) * (-1 + 2 * m) * pi)
  sigma_0 <- (-(3/4)) * (tan((1/2) * (-1 + 2 * m) * pi) - tan((1/2) * pi * (-1 + 2 * qt)))
}


if ( fd == "logit" & sd == "t2")
{
  mu_0 <- -sqrt((-1 + 4 * m - 4 * m^2)/(2 * (-m + m^2)))
  sigma_0 <- (sqrt(2)*sqrt((-(-1 + qt))*qt*(-1 + 2*qt)^2) + 2*qt*mu_0 -2*qt^2*mu_0)/(2*qt- 2*qt^2)
}


  if (sd == "km" & fd == "km")
  {
    sigma_0 = 0.1
    mu_0= 0.1
  }
  
  start = as.numeric(c(mu_0, sigma_0))
  
  if(skew){
    theta_0 <- (mean(ydata) - median(ydata))/sd(ydata)
    sigma_0 <- log(sd(ydata))
    start = as.numeric(c(mu_0, sigma_0))
    
  }
  start
  
} 
