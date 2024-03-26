#' @title Simulate (generalized) linear network model data
#' 
#' @description
#' Generates data based on a (generalized) linear network model, where the
#' dependent variable can be either continuous or dichotomous. All the 
#' independent variables are dichotomous.
#' 
#' @param n The size of the vertex set (|V(G)|) for the random graphs.
#' @param m The number of graphs to generate.
#' @param family Family of the generalized linear model used to generate the 
#' data. The available families are "gaussian" and "binomial".
#' @param intercept Logical; should the intercept be simulated? 
#' \code{FALSE} by default.  
#' @param beta Values of the true coefficients used to generate the data. If 
#' \code{intercept} is \code{TRUE}, the first value in \code{beta} is considered
#' to be the intercept. If not provided, the coefficients are randomly sampled 
#' from a standard normal distribution and the intercept is defined as the maximum 
#' between \eqn{m/2} and the median of the beta coefficients. \code{NULL} by default.
#' @param red.var Value of the residual variance used when the family is "gaussian". 
#' If not provided, the residual variance is set to 1. 
#' @param seed Integer used as the seed in the generation process. Allows the 
#' user to generate the same data consistently.
#' @param ... other arguments passed on to the function \code{rgraph} of the 
#' \code{sna} package. 
#' 
 

gen.netglm <- function(n, m, family = "gaussian", intercept = FALSE, 
                       beta = NULL, red.var = NULL, seed = NULL,
                       ...) {
  # We need sna to generate random networks. 
  # If we want this function to be available to the public, we would need to 
  # import or suggest sna in the description file. Probably adding it to suggest
  # is the best option. I leave it here for now.
  # Same with abind
  require(sna)
  require(abind)
  
  # Set seed if provided
  if (!is.null(seed)) {set.seed(seed)}
  
  # Check whether family is gaussian or binomial
  if (!family %in% c("gaussian", "binomial")) {
    stop("This function only supports the following families:",
         " 'gaussian' and 'binomial'.")
  }
  
  # Generate beta coefficients if not provided.
  if (is.null(beta)) {
    beta <- rnorm(m)
    if (intercept) {
      beta <- c(-max(m/2, median(beta)), beta)
    }
  }
  
  if (intercept) {
    names(beta) <- c("intercept", paste0("x", 1:m))
  } else {
    names(beta) <- paste0("x", 1:m)
  }
  
  # Define residual variance for gaussian model
  if (is.null(red.var) && family == "gaussian") {
    red.var <- red.sd <- 1
  } else {
      red.sd <- sqrt(red.var)
      }
  if (!is.null(red.var) && family == "binomial") {
    red.var <- red.sd <- NULL
    message("Dispersion parameter for binomial family taken to be 1, ",
            "The argument 'red.var' = ", red.var, " is ignored.")
  }
  
  # Generate networks given the function rgraph of sna package
  gen.rgraph <- rgraph(n, m, ...)
  gen.rgraph.original <- gen.rgraph
  
  if (intercept) {
    gen.rgraph <- abind::abind(array(1, dim = c(1, n, n)), gen.rgraph, along = 1)
  }
  
  # Compute linear model (predicted values)
  xb <- gen.rgraph * beta
  xb <- apply(xb, c(2, 3), sum)
  
  # Compute dependent variable according to the family
  if (family == "gaussian") {
    y <- xb + rnorm(n * n, 0, red.sd)
  }
  
  if (family == "binomial") {
    y <- rbinom(n * n, 1, plogis(xb))
    y <- matrix(y, nrow = n, ncol = n)
  }
  
  # Bind dependent and independent variables into one array. 
  data <- abind::abind(gen.rgraph.original, y, along = 1)
  dimnames(data) <- list(vars     = c(paste0("x", 1:m), "y"),
                         sender   = paste0("sdr", 1:n),
                         receiver = paste0("rcr", 1:n))
  
  # Organize generated data and true parameters into a list.   
  out <- list(data         = data,
              beta.true    = beta,
              red.var.true = red.var,
              red.sd.true  = red.sd)
  
  # Drop unused elements.
  out <- out[!sapply(out, is.null)]
  
  return(out)
}



