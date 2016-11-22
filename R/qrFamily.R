#' @title Overview of the family of distributions
#' 
#' @aliases cdfqrFamily
#' @description The cdfquantreg family consists of the currently available distributions that can be used to fit quantile regression models via the cdfquantreg() function.  
#' @param shape To show all distributions or the set of distribution for a specific type of shape. Can be \code{BM}, \code{TM},\code{LL} or \code{FT} for Bimodal, Trimodal, Logit-logistic or Finite-tailed shapes, respectively.
#'  
#' @details The cdfquantreg package includes a two-parameter family of distributions for 
#' modeling random variables on the (0, 1) interval by applying the cumulative 
#' distribution function (cdf) of one \dQuote{parent} distribution to the 
#' quantile function of another. \cr
#' The naming of these distributions is \dQuote{parent - child} or 
#' \dQuote{fd - sd}, where \dQuote{fd} is the parent distribution, and \dQuote{sd} 
#' is the child distribution. \cr
#' The distributions have four characteristic shapes: Logit-logistic, bimodal, trimodal, and finite-tailed. 
#' Here is the list of currently available distributions.
#' 
#' \bold{Bimodal Shape Distributions}
#'   \tabular{lllc}{
#'   \bold{Distribution}  \tab \bold{R input} \tab \bold{Alternative Input}  \tab \bold{Shape}\cr
#'   Burr VII-ArcSinh	\tab \code{fd = "burr7", sd = "arcsinh"}	\tab \code{family = "burr7-arcsinh"} \tab	Bimodal	\cr
#'   Burr VII-Cauchy	\tab \code{fd = "burr7", sd = "cauchy"}	\tab \code{family = "burr7-cauchy"} \tab	Bimodal	\cr
#'   Burr VII-T2	\tab \code{fd = "burr7", sd = "t2"}	\tab \code{family = "burr7-t2"} \tab	Bimodal	\cr
#'   Burr VIII-ArcSinh	\tab \code{fd = "burr8", sd = "arcsinh"}	\tab \code{family = "burr8-arcsinh"} \tab	Bimodal	\cr
#'   Burr VIII-Cauchy	\tab \code{fd = "burr8", sd = "cauchy"}	\tab \code{family = "burr8-cauchy"} \tab	Bimodal	\cr
#'   Burr VIII-T2	\tab \code{fd = "burr8", sd = "t2"}	\tab \code{family = "burr8-t2"} \tab	Bimodal	\cr
#'   Logit-ArcSinh	\tab \code{fd = "logit", sd = "arcsinh"}	\tab \code{family = "logit-arcsinh"} \tab	Bimodal	\cr
#'   Logit-Cauchy	\tab \code{fd = "logit", sd = "cauchy"}	\tab \code{family = "logit-cauchy"} \tab	Bimodal	\cr
#'   Logit-T2	\tab \code{fd = "logit", sd = "t2"}	\tab \code{family = "logit-t2"} \tab	Bimodal	\cr
#'   T2-ArcSinh	\tab \code{fd = "t2", sd = "arcsinh"}	\tab \code{family = "t2-arcsinh"} \tab	Bimodal	\cr
#'   T2-Cauchy	\tab \code{fd = "t2", sd = "cauchy"}	\tab \code{family = "t2-cauchy"} \tab	Bimodal	\cr
#'   }
#'   
#'   
#' \bold{Trimodal Shape Distributions}
#'   \tabular{lllc}{
#'   \bold{Distribution}  \tab \bold{R input} \tab \bold{Alternative Input}  \tab \bold{Shape}\cr
#'   	ArcSinh-Burr VII	\tab \code{fd = "arcsinh", sd = "burr7"}	\tab \code{family = "arcsinh-burr7"} \tab	Trimodal	\cr
#'   	ArcSinh-Burr VIII	\tab \code{fd = "arcsinh", sd = "burr8"}	\tab \code{family = "arcsinh-burr8"} \tab	Trimodal	\cr
#'   	ArcSinh-Logistic	\tab \code{fd = "arcsinh", sd = "logistic"}	\tab \code{family = "arcsinh-logistic"} \tab	Trimodal	\cr
#'   	ArcSinh-T2	\tab \code{fd = "arcsinh", sd = "t2"}	\tab \code{family = "arcsinh-t2"} \tab	Trimodal	\cr
#'   	Cauchit-Burr VII	\tab \code{fd = "cauchit", sd = "burr7"}	\tab \code{family = "cauchit-burr7"} \tab	Trimodal	\cr
#'   	Cauchit-Burr VIII	\tab \code{fd = "cauchit", sd = "burr8"}	\tab \code{family = "cauchit-burr8"} \tab	Trimodal	\cr
#'   	Cauchit-Logistic	\tab \code{fd = "cauchit", sd = "logistic"}	\tab \code{family = "cauchit-logistic"} \tab	Trimodal	\cr
#'   	Cauchit-T2	\tab \code{fd = "cauchit", sd = "t2"}	\tab \code{family = "cauchit-t2"} \tab	Trimodal	\cr
#'   	T2-Burr VII	\tab \code{fd = "t2", sd = "burr7"}	\tab \code{family = "t2-burr7"} \tab	Trimodal	\cr
#'   	T2-Burr VIII	\tab \code{fd = "t2", sd = "burr8"}	\tab \code{family = "t2-burr8"} \tab	Trimodal	\cr
#'   	T2-Logistic	\tab \code{fd = "t2", sd = "logistic"}	\tab \code{family = "t2-logistic"} \tab	Trimodal	\cr
#'   }
#' 
#' \bold{Logit-logistic Shape Distributions}
#'   \tabular{lllc}{
#'   \bold{Distribution}  \tab \bold{R input} \tab \bold{Alternative Input}  \tab \bold{Shape}\cr
#'   	Burr VII-Burr VII	\tab \code{fd = "burr7", sd = "burr7"}	\tab \code{family = "burr7-burr7"} \tab	Logit-logistic	\cr
#'   	Burr VII-Burr VIII	\tab \code{fd = "burr7", sd = "burr8"}	\tab \code{family = "burr7-burr8"} \tab	Logit-logistic	\cr
#'   	Burr VII-Logistic	\tab \code{fd = "burr7", sd = "logistic"}	\tab \code{family = "burr7-logistic"} \tab	Logit-logistic	\cr
#'   	Burr VIII-Burr VII	\tab \code{fd = "burr8", sd = "burr7"}	\tab \code{family = "burr8-burr7"} \tab	Logit-logistic	\cr
#'   	Burr VIII-Burr VIII	\tab \code{fd = "burr8", sd = "burr8"}	\tab \code{family = "burr8-burr8"} \tab	Logit-logistic	\cr
#'   	Burr VIII-Logistic	\tab \code{fd = "burr8", sd = "logistic"}	\tab \code{family = "burr8-logistic"} \tab	Bimodal	\cr
#'   	Logit-Burr VII	\tab \code{fd = "logit", sd = "burr7"}	\tab \code{family = "logit-burr7"} \tab	Logit-logistic	\cr
#'   	Logit-Burr VIII	\tab \code{fd = "logit", sd = "burr8"}	\tab \code{family = "logit-burr8"} \tab	Logit-logistic	\cr
#'   	Logit-Logistic	\tab \code{fd = "logit", sd = "logistic"}	\tab \code{family = "logit-logistic"} \tab	Logit-logistic	\cr
#'   }  
#' 
#' \bold{Finite-tailed Shape Distributions}
#'   \tabular{lllc}{
#'   \bold{Distribution}  \tab \bold{R input} \tab \bold{Alternative Input}  \tab \bold{Shape}\cr
#'   	ArcSinh-ArcSinh	\tab \code{fd = "arcsinh", sd = "arcsinh"}	\tab \code{family = "arcsinh-arcsinh"} \tab	Finite-tailed	\cr
#'   	ArcSinh-Cauchy	\tab \code{fd = "arcsinh", sd = "cauchy"}	\tab \code{family = "arcsinh-cauchy"} \tab	Finite-tailed	\cr
#'   	Cauchit-ArcSinh	\tab \code{fd = "cauchit", sd = "arcsinh"}	\tab \code{family = "cauchit-arcsinh"} \tab	Finite-tailed	\cr
#'   	Cauchit-Cauchy	\tab \code{fd = "cauchit", sd = "cauchy"}	\tab \code{family = "cauchit-cauchy"} \tab	Finite-tailed	\cr
#'   	T2-T2	\tab \code{fd = "t2", sd = "t2"}	\tab \code{family = "t2-t2"} \tab	Finite-tailed	\cr
#'   }  
#'   
#' \bold{Kumaraswamy Distribution}
#'   \tabular{lllc}{
#'   \bold{Distribution}  \tab \bold{R input} \tab \bold{Alternative Input}  \tab \bold{Shape}\cr
#'   	Kumaraswamy	\tab \code{fd = "", sd = ""}	\tab \code{family = "-"} \tab		\cr
#'   } 
#'   
#' @return A list of distributions that are available in the current version of package.
#' @export
#'
#' @examples
#' cdfqrFamily()
#' 
cdfqrFamily <- function(shape = "all") {
  db <- c('ArcSinh-ArcSinh', 'ArcSinh-BurrVII', 'ArcSinh-BurrVIII', 
          'ArcSinh-Cauchy', 'ArcSinh-Logistic', 'ArcSinh-T2',
          'BurrVII-ArcSinh', 'BurrVII-BurrVII', 'BurrVII-BurrVIII', 
          'BurrVII-Cauchy', 'BurrVII-Logistic', 'BurrVII-T2',
          'BurrVIII-ArcSinh', 'BurrVIII-BurrVII', 'BurrVIII-BurrVIII', 
          'BurrVIII-Cauchy', 'BurrVIII-Logistic', 'BurrVIII-T2',
          'Cauchit-ArcSinh', 'Cauchit-BurrVII', 'Cauchit-BurrVIII', 
          'Cauchit-Cauchy', 'Cauchit-Logistic', 'Cauchit-T2',
          'Logit-ArcSinh', 'Logit-BurrVII', 'Logit-BurrVIII', 
          'Logit-Cauchy', 'Logit-Logistic', 'Logit-T2', 
          'T2-ArcSinh', 'T2-BurrVII', 'T2-BurrVIII', 'T2-Cauchy', 
          'T2-Logistic', 'T2-T2','Kumaraswamy')
  
  fds <- c(rep(c("arcsinh","burr7","burr8","cauchit","logit","T2"), each = 6),"km")
  
  sds <- c(rep(c("arcsinh","burr7","burr8","cauchy","logistic","T2"), 6),"km")
  
  shapes <- c("Finite-tailed","Trimodal","Trimodal","Finite-tailed","Trimodal",
             "Trimodal","Bimodal","Logit-logistic","Logit-logistic","Bimodal",
             "Logit-logistic","Bimodal","Bimodal","Logit-logistic",
             "Logit-logistic","Bimodal","Bimodal","Logit-logistic","Finite-tailed",
             "Trimodal","Trimodal", "Finite-tailed","Trimodal","Trimodal",
             "Bimodal","Logit-logistic","Logit-logistic","Bimodal",
             "Logit-logistic","Bimodal","Bimodal", "Trimodal","Trimodal",
             "Bimodal","Trimodal","Finite-tailed","")

  dist <- data.frame(Distributions = db,
                     fd = fds, sd = sds, shape=shapes, stringsAsFactors = F)
  
  if (shape =="FT")  dists <- subset(dist, shape =="Finite-tailed") else
    if  (shape =="BM")  dists <- subset(dist, shape =="Bimodal") else
      if (shape =="TM")  dists <- subset(dist, shape =="Trimodal") else
        if (shape =="LL")  dists <- subset(dist, shape =="Logit-logistic") else
          dists <- dist
  
  cat("Overview cdfquantreg distributions:")
  knitr::kable(dists, row.names=F)
}
 
