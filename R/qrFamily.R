#' @title Overview of the family of distributions
#' 
#' @aliases cdfqrFamily
#' @description The cdfquantreg family consists of the currently available distributions that can be used to fit quantile regression models via the cdfquantreg() function.  
#' 
#' @details The cdfquantreg package includes a two-parameter family of distributions for 
#' modeling random variables on the (0, 1) interval by applying the cumulative 
#' distribution function (cdf) of one \dQuote{parent} distribution to the 
#' quantile function of another. \cr
#' The naming of these distributions are \dQuote{parent - child} or 
#' \dQuote{fd - sd}, where \dQuote{fd} is the parent distribution, and \dQuote{sd} 
#' is the child distribution. \cr
#' The distributions have four characteristic shapes: Logit-logistic, bimodal, trimodal, and fat-tailed. 
#' Here is the list of currently available distributions.
#' 
#'   \tabular{lllc}{
#'   \bold{Distribution}  \tab \bold{R input} \tab \bold{Alternative Input}  \tab \bold{Shape}\cr
#'   ArcSinh - Burr VIII \tab \code{fd = "arcsinh", sd = "burr8"} \tab \code{family = "arcsinh-burr8"} \tab trimodal\cr
#'   ArcSinh - Cauchy    \tab \code{fd = "arcsinh", sd = "cauchy"}\tab \code{family = "arcsinh-cauchy"} \tab fat-tailed\cr
#'   ArcSinh - T2        \tab \code{fd = "arcsinh", sd = "t2"}    \tab \code{family = "arcsinh-t2"} \tab trimodal\cr
#'   Burr VIII - Burr VII \tab \code{fd = "burr8", sd = "burr7"}   \tab \code{family = "burr8-burr7"} \tab fat-tailed\cr
#'   Burr VIII - Burr VIII\tab \code{fd = "burr8", sd = "burr8"}   \tab \code{family = "burr8-burr8"} \tab fat-tailed\cr
#'   Burr VIII - Cauchy  \tab \code{fd = "burr8", sd = "cauchy"}  \tab \code{family = "burr8-cauchy"} \tab bimodal\cr
#'   Burr VIII - T2      \tab \code{fd = "burr8", sd = "t2"}      \tab \code{family = "burr8-t2"} \tab bimodal\cr
#'   Logit - Burr VIII \tab \code{fd = "logit", sd = "burr8"}   \tab \code{family = "logit-burr8"} \tab bimodal\cr
#'   Logit - Cauchy \tab \code{fd = "logit", sd = "cauchy"}   \tab \code{family = "logit-cauchy"} \tab bimodal\cr
#'   Logit - Logistic  \tab \code{fd = "logit", sd = "logistic"}  \tab \code{family = "logit-logistic"} \tab logit-logistic\cr
#'   Logit - T2      \tab \code{fd = "logit", sd = "t2"}      \tab \code{family = "logit-t2"} \tab bimodal\cr
#'   T2 - Burr VII \tab \code{fd = "t2", sd = "burr7"}   \tab \code{family = "t2-burr7"} \tab trimodal\cr
#'   T2 - Burr VIII \tab \code{fd = "t2", sd = "burr8"}   \tab \code{family = "t2-burr8"} \tab trimodal\cr
#'   T2 - Cauchy  \tab \code{fd = "t2", sd = "cauchy"}  \tab \code{family = "t2-cauchy"} \tab bimodal\cr
#'   T2 - T2      \tab \code{fd = "t2", sd = "t2"}      \tab \code{family = "t2-t2"} \tab fat-tailed\cr
#'   Kumaraswamy      \tab \code{fd = "km", sd = "km"}      \tab \code{family = "km-km"} \tab  \cr
#'   }
#'   
#' @return A list of distributions that are available in the current version of package.
#' @export
#'
#' @examples
#' cdfqrFamily()
cdfqrFamily <- function() {
  db <- c('ArcSinh-BurrVIII', 'ArcSinh-Cauchy', 'ArcSinh-T2',
          'BurrVIII-BurrVIII', 'BurrVIII-BurrVIII', 'BurrVIII-Cauchy','BurrVIII-T2',
          'Logit-BurrVIII', 'Logit-Cauchy', 'Logit-Logistic', 'Logit-T2',
          'T2-BurrVIII', 'T2-BurrVIII', 'T2-Cauchy', 'T2-T2', 'Kumaraswamy')
  
  fds <- c("arcsinh","arcsinh","arcsinh","burr8","burr8","burr8","burr8",
           "logit","logit","logit","logit",
           "t2","t2","t2","t2","km")
  sds <- c("burr8","cauchy","t2","burr7","burr8","cauchy","t2",
           "burr8","cauchy","logistic","t2","burr7","Burr8","cauchy","t2","km")
  
  dist <- data.frame(Distributions = db,
                     fd = fds, sd = sds)
  cat("Overview cdfquant distributions: \n \n")
  print(dist, row.names=F)
}
 