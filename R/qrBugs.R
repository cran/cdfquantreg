#' @title Running cdf quantile regression in OpenBUGS
#' @aliases qrBugs
#' @description Function used for running cdf quantile regression in OpenBUGS and performing 
#' Bayesian MCMC estimation.
#' In addition, a simple random effects model is allowed to be estimated in this function.
#' 
#' @param formula A formula object, with the DV on the left of an ~ operator, and predictors on the right. For the part on the right of '~', the specification of the submodels can be seperated by '|'. So \code{y ~ X1 | X2} means the DV is \code{y},\code{X1} is the term in the mean submodel, and \code{X2} is the term in the dispersion submodel.
#' @param data The data file. Can be either a data.frame or a list of with specifying the names of variable in the list. The structure of the list and list component (e.g., matrix) should follow BUGS data format (see XXX for more details about preparing for list format of the data) 
#' @param fd A string that specifies the parent distribution.
#' @param sd A string that specifies the sub-family distribution.
#' @param bugs A logical value to indicate whether to use JAGS or OpenBUGS (NOTE: currently only OpenBUGS supports the customized likelihood functions; the package will update to allow JAGS as soon as JAGS provides that function).
#' @param random Character or vector of characters that indicate the random effect factors.
#' @param nodes Character or vector of characters that indicate the parameters to be estimated. The default nodes are the coefficients in both the mean and dispersion submodel.
#' @param inits A list of values that serve as initial values for MCMC chain procedure.
#' @param modelname Name of the model (optional).
#' @param working.directory the directory for generating temporary BUGS files.
#' @param n.iter The number of MCMC samples to be drawn from the posterior distribution.
#' @param n.burnin The number of samples to be discarded when summarizing the MCMC simulation results.
#' @param ... further arguments to \code{\link[R2OpenBUGS]{bugs}}.
#' @return A bugs object (See more details)
#' 
#' @import Formula
#' 
#' @seealso \code{\link[R2OpenBUGS]{bugs}}
#' @examples
#' data(cdfqrExampleData)
#'\dontrun{
#' # Need to OpenBUGS has been installed, and R2OpenBUGS has been loaded first. 
#' library(R2OpenBUGS)
#' bugfit <- qrBugs(crc99 ~ vert | confl, data = JurorData, 't2', 't2',clearWD=TRUE)
#' bugfit
#' Inference for Bugs model at "bugmodel.txt", 
#' # Current: 2 chains, each with 10000 iterations (first 5000 discarded)
#' # Cumulative: n.sims = 10000 iterations saved
#' # mean  sd  2.5%   25%   50%   75% 97.5% Rhat n.eff
#' # b_0        0.8 0.1   0.6   0.7   0.8   0.9   1.0    1   830
#' # d_0       -0.2 0.1  -0.4  -0.3  -0.2  -0.1   0.1    1  1700
#' # b_vert     0.1 0.1  -0.1   0.0   0.1   0.2   0.3    1   170
#' # d_confl    0.0 0.1  -0.3  -0.1   0.0   0.0   0.2    1  4000
#' # deviance -49.3 2.8 -52.8 -51.3 -49.9 -47.9 -42.1    1  7400
#' # 
#' # For each parameter, n.eff is a crude measure of effective sample size,
#' # and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
#' # 
#' # DIC info (using the rule, pD = Dbar-Dhat)
#' # pD = 4.0 and DIC = -45.3
#' # DIC is an estimate of expected predictive error (lower deviance is better).
#'}

# qrBugs <- function(formula, data, fd, sd, bugs=TRUE, random = NULL, 
#                    nodes = NULL, inits = NULL,n.iter = 10000, 
#                    n.burnin = n.iter/2, modelname = "bugmodel",
#                    working.directory = NULL, ...) {
#   
#   # Processing input to identify the distribution
#   fd <- tolower(fd)
#   sd <- tolower(sd)
#   
#   # check the formula
#   formula <- Formula::as.Formula(formula)
#   if (length(formula)[2] < 2) {
#     # if the user only specifies the mean submodel, add the dispersion submodel intercept in to
#     # the formula
#     formula <- deparse(formula)
#     formula <- paste(formula, "| 1", sep = "")
#     formula <- Formula::as.Formula(formula)
#   }
#   formula.string <- deparse(formula)  # Turn the formula to a string
#   
#   # Process the data file
#   if (is.data.frame(data)) {
#     
#     for (i in 1:ncol(data)) {
#       if (!is.numeric(data[, i])) 
#         data[, i] <- factor(data[, i])
#     }  #turn character variables to factor 
#     
#     data_temp <- model.frame(formula, data = data)  # get the data with variables in the model only
#     data_temp <- data.frame(data_temp)
#     
#     for (i in 1:ncol(data_temp)) {
#       # processing non-numeric type variables to factors
#       if (!is.numeric(data_temp[, i])) {
#         data_temp[, i] <- factor(data_temp[, i])  # turn character variables to factors
#         
#         # Get the names dummy variables generated by dummy coding the factor variables
#         # (newname = oldname + factor level)
#         var.name <- names(data_temp[, i])
#         var.name.new <- paste(var.name, levels(data_temp[, i]), sep = "")
#         var.name.new <- paste("(", paste(var.name.new, collapse = " + "), 
#           ")", sep = "")
#         
#         # For those variables, replace the single variable name in the formula with the
#         # anmes of dummy variables
#         formula.string <- sub(var.name, formula.string)
#       }
#     }
#     
#     # Generate data set that factor variables have been dummy coded
#     ydata <- as.matrix(model.frame(formula, data, rhs = 0))
#     if (any(ydata <= 0) | any(ydata >= 1)) {
#       warning(paste("The values of the dependent variable is rescaled into (0, 1) interval", "\n", "see scaleTR() for details"))
#       ydata <- scaleTR(ydata)
#     }
#     
#     data_lm <- model.matrix(formula, data = data_temp, rhs = 1)  # datamatrix for location model without intercept
#     data_pm <- model.matrix(formula, data = data_temp, rhs = 2)  # datamatrix for dispersion model
#    
#     data_temp <- cbind(data_lm, data_pm)
#     
#     data_temp <- data_temp[, c(-1, -(ncol(data_lm)+1))]
#     data_temp <- cbind(as.numeric(ydata), data_temp)
#     colnames(data_temp)[1] <- 'y'
#     names <- unique(colnames(data_temp))
#     data_temp <- data_temp[, names]
#     
#     data_temp <- data.frame(data_temp, check.names = F)
#     names(data_temp) <- names
#     # If there are random factors:
#     random_level <- NULL
#     
#     if (!is.null(random)) {
#       # Get the variables that indexing the random factor
#       data_random <- data.frame(data[, random])
#       
#       # Get the number of levels in the random factors for each of the random factor
#       random_level <- matrix(0, nrow = 1, ncol = length(random))
#       for (i in 1:length(random)) {
#         data_random[, i] <- factor(data_random[, i])
#         random_level[, i] <- length(levels(data_random[, i]))
#         
#         # Replace the factor type variables with numerica types of variable
#         data_random[, i] <- as.numeric(data_random[, i])
#       }
#       
#       # Prepare the indexing variables for Bugs data file
#       names(random_level) <- paste("J_", random, sep = "")
#       
#       random_level <- as.list(random_level)
#       
#       data_random <- data.frame(data_random)
#       names(data_random) <- random
#       
#       # Combine the dataset generated in model.matrix method and the random factor
#       # variables
#       data_temp <- cbind(data_temp, data_random)
#     }
#     
#     names(data_temp) <- gsub(":", "_", x = names(data_temp))
#     data_list <- as.list(data_temp)  # Generate the list data for Bugs processing
#     
#     # Put case number into the dataset
#     data_list <- c(data_list, N = nrow(data_temp))
#     
#     # Put index number of random factors into the dataset
#     data_list <- c(data_list, random_level)
#     data <- data_list
#     
#     formula <- Formula::as.Formula(formula.string)
#   }
#   
#   
#   # set up working directory
#   if (is.null(working.directory)) 
#     working.directory <- getwd()
#   
#   
#   # Generate a default model
#   bugs.model <- bugsModel(formula, fd, sd, random = random, modelname = modelname, wd = working.directory)
#   
#   # if sampling node is not provided, use all nodes in the model
#   if (is.null(nodes)) {
#     nodes <- bugs.model$nodes_sample
#   }
#   
#   # if init is not provided, use the default init values
#   if (is.null(inits)) {
#     init1 <- bugs.model$init1
#     init2 <- bugs.model$init2
#   }
#   
#   if(!bugs){
#     # bugs.m <- R2jags::jags(data = data, 
#     #                        inits = list(init1, init2), 
#     #                        n.iter = n.iter,
#     #                n.burnin = n.burnin, n.thin = n.thin, parameters.to.save = nodes,
#     #                working.directory = wd,
#     #                model.file = paste(modelname, ".txt", sep = ""),
#     #                 n.chains = 2)
#     
#   
#     stop(paste("Sorry, currently only OpenBUGS can be used."))
#     
#   }else if (!requireNamespace("R2OpenBUGS", quietly = TRUE)) {
#     stop("R2OpenBUGS and OpenBUGS are needed for using the current version of this function. Please install and load R2OpenBUGS before using this function.",
#          call. = FALSE)
#   }else{
# 
#     bugs.m <- R2OpenBUGS::bugs(data = data, inits = list(init1, init2),  
#                    n.iter = n.iter, parameters.to.save = nodes, n.chains = 2, 
#                    model.file = paste(modelname, ".txt", sep = ""),
#                    n.burnin = n.burnin,  
#                    working.directory = working.directory, ...)
# 
#   }
# 
#   return(bugs.m)
# } 
