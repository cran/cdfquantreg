#' @title Influence Diagnosis For Fitted Cdfqr Object
#' @description Influence Diagnosis (dfbetas) For Fitted Cdfqr Object
#' @aliases influence.cdfqr
#' @param model A cdfqr model object 
#' @param type A string that indicates whether the results for all parameters are to be returned, or only the submodel's parameters returned.
#' @param what for influence statistics based on coefficient values, indicate the predictor variables that needs to be tested. 
#' @param plot if plot is needed.
#' @param method Currently only 'dfbeta' method is available.
#' @param id for plot only, if TRUE, the case ids will be displayed in the plot.
#' @param ... Pass onto other functions or currently ignored
#' @return A matrix, each row of which contains the estimated influence on parameters when that row's observation is removed from the sample.
#' 
#' @examples
#' data(cdfqrExampleData)
#'fit <- cdfquantreg(crc99 ~ vert | confl, 't2', 't2', data = JurorData)
#'#It takes some time especially the data is large.
#'influcne <- influence(fit)
#'plot(influcne[,2])
#'
#'\dontrun{
#' # Same as influence(fit)
#'dfbetval <- dfbetas(fit)
#'}
#'
#' @seealso \code{\link[stats]{lm.influence}}, \code{\link[stats]{influence.measures}}
#' @import stats
#' @method influence cdfqr
#' @import graphics
#' @export
influence.cdfqr <- function(model, method = "dfbeta",
                            type = c("full","location", "dispersion","skew",
                                     "zero", "one"), 
                            what = "full",plot=FALSE, id = FALSE, ...) {
  # - dfbeta residuals (for influence check)
 object <- model
  
  method <- match.arg(method)
  influcne <- switch(method, dfbeta = {
    infl.stats <- dfbeta(object, type, what)
  })
  
  plotfun <- function(ddd, id=FALSE){
    if (!is.null(ncol(ddd))) {
      par(mfcol=c(2,2), mar = c(2,2,2,1), oma = c(2, 1, 0, 0))
      
      for (i in 1:ncol(ddd))
      {
        plot(ddd[,i],xlab = "", ylab = "Estimated")
        if(id) {text(ddd[,i], labels=1:nrow(ddd),cex=1, pos=2, col="red")}
        title(colnames(ddd)[i], line = 0.5)
        mtext("Cases", side = 1, outer = T, line = 0)
      }
      
    }
  }
  if (plot){
    infplot <- plotfun(influcne)
    return(list(influcne, infplot))
  }else{
    return(influcne)
  }


}

#' @method dfbeta cdfqr
#' @export
#' @rdname influence.cdfqr
dfbeta.cdfqr <- function(model, type = c("full","location", "dispersion","skew",
                                         "zero", "one"), 
                          what = "full", ...) {
  dfbetas(model, type = type, what = what)
  }

#' @method dfbetas cdfqr
#' @export
#' @rdname influence.cdfqr
dfbetas.cdfqr <- function(model, type = c("full","location", "dispersion","skew", "zero", "one"), what = "full",...) {
  object <- model
  # - dfbeta residuals (for influence check)
  type <- match.arg(type)
  
  # The original model estimation call
  call <- object$call
  fd1 <- object$family$fd
  sd1 <- object$family$sd
  
  # Get the data from the call
  dat <- eval(call$data)
  n <- nrow(dat)
  
  # the number of location parameters
  k_lm <- nrow(object$coefficients$location)
  
  #Get the original coefficients estimates (mean and se)
  coefficients <-  do.call(rbind, object$coefficients)
  coef0 <- coef(object)
  bse <- coefficients[, 2]
  
  # Casewise deletion
  betas <- NULL
  for (i in 1:nrow(dat)) {
    dat1 <- dat[-i, ]

    #modify the call for cdfqr function with the new dataset
    mod <- update(object, .~., fd=fd1, sd=sd1, data = dat1, start = coef0)
    
    dfbeta <- (coef(mod)[names(coef0)] - coef0)/bse
    betas <- rbind(betas, dfbeta)
  }
  
  colnames(betas) <- names(coef0)
  betas_dispersion <- betas[, grep("sigma", colnames(betas))] 
  
  if(!is.null(object$coefficients$location)){
    betas_location <- betas[, grep("mu", colnames(betas))] 
    colnames(betas_location) <- paste("lm.",colnames(betas_location),sep="")
    
  }
  if("cdfqrH"%in%class(object)){
    betas_zero <- betas[, grep("zero", colnames(betas))]
    betas_one <- betas[, grep("one", colnames(betas))]
  }

  if("cdfqrFT"%in%class(object)){
    betas_skew <- betas[, grep("theta", colnames(betas))]
  }
  
  # If only a specific subset of parameters is needed, extract such subsets
  if (what != "full") {
    if(!is.null(object$coefficients$location)){
      locind <- pmatch(colnames(betas_location), what, duplicates.ok = TRUE) 
      if(length(na.omit(locind))!=0) {
        betas_location<- betas_location[, what]
      }else{betas_location <- NULL}
    }else{betas_location <- NULL}
   
    precind <- pmatch(colnames(betas_dispersion), what, duplicates.ok = TRUE) 
    
    if(length(na.omit(precind))!=0) {
      betas_dispersion<- betas_dispersion[, what]
    }else{betas_dispersion <- NULL}
   
    
    if("cdfqrH"%in%class(object)){
      zeroind <- pmatch(colnames(betas_zero), what, duplicates.ok = TRUE) 
      if(length(na.omit(zeroind))!=0) {
        betas_zero<- betas_zero[, what]
      }else{betas_zero <- NULL}
      
      oind <- pmatch(colnames(betas_one), what, duplicates.ok = TRUE) 
      if(length(na.omit(oind))!=0) {
        betas_one<- betas_one[, what]
      }else{betas_one <- NULL}
      
    } 
  
   if("cdfqrFT"%in%class(object)){
     sind <- pmatch(colnames(betas_skew), what, duplicates.ok = TRUE) 
     if(length(na.omit(sind))!=0) {
       betas_skew<- betas_skew[, what]
     }else{betas_skew <- NULL}
   }
    
  }
  
  #Rename the two submodels' paraters to distinguish the two submodels
  colnames(betas_dispersion) <- paste("pm.",colnames(betas_dispersion),sep="")
  
  betas<-cbind(betas_location, betas_dispersion)
  
  betas <- switch(type, full = {
    betas
  }, location = {
    betas_location
  },  dispersion = {
    betas_dispersion
  },  skew = {
    betas_skew
  },  zero = {
    betas_zero
  },  one = {
    betas_one
  })
  
  return(betas)
  
} 

