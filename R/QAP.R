#' MRQAP for unnested data
#'
#' Estimates a Multiple Regression Quadratic Assignment Proccedure model 
#' (MRQAP; Krackhardt, 1988). MRQAPs allow investigating associations between
#'  characteristics of dyads in networks (e.g., the level of homophily 
#'  between two actors) and a binary or continuous
#'   tie variable (e.g., friendship, amount of time spent together). 
#'
#' @param dv a matrix with n * m dimensions and cells indicating the presence of a tie (1 = tie, 0 = no tie for binary variables) 
#'  or the weight of a tie (for continous tie variables) charcterizing the dependent variable
#' @param iv1 a list of matrices with n * m dimensions characterizing the independent variables
#' @param iv.list.per lists in the ivs argument are should be nested by group and independent matrices, if this is not the case (grouped by independent matrices and then groups)
#' the argument iv.list.per = "iv" can be used to restructure the data.
#' @param family family of the generalized linear model. default is "gaussian" for continuous dependent varaibles. F
#' or binday dependent variables "binomial" is advised.
#' @param iv.names names of the independent variables for the output object
#' @param mode permutation method to be applied. default is "yQAP" for permuting the Y / dv
#'  variables."dspQAP" applies Dekker's semi partialing method (Dekker, Krackhard, & Snijders, 2007) 
#' @param samples  number of permutations, default is 1000.
#' @param diag  boolean for using the diagonal values of matrices in the estimation. default is FALSE
#' @param directed "directed" if the dependent network is directed (ties from A to B and B to A are possible), 
#' "undirected" if the dependent network is undirected (ties from A to B are identical to B to A). Default is "directed".
#' @param round.to numeric, numer of digits in output table
#' @param cpu number of cpu's to be used for estimation, default is 1
#' @param logfilename name of log file printing intermediate reports during the estimation procedure.
#' @param verbose reports of what is happening under the hood during the call of the function, default is TRUE
#' @param global.deltas during "dspQAP" estimation, should global or local delta values be used. default is TRUE
#' @param return.perms should permuted networks be part of the output? default is FALSE
#'
#' @examples 
#' # create test data #
#' inspired by the example funciton in sna::netlm
#' ivnet1<-sna::rgraph(20,4)
#'
#' dv1<-ivnet1[1,,]+4*ivnet1[2,,]+2*ivnet1[3,,]   # Note that the fourth graph is unrelated
#' dv1 <- dv1 + rnorm(400,mean = 1, sd = 1)

#' iv1 <- list(ivnet1[1,,],ivnet1[2,,],ivnet1[3,,], ivnet1[4,,])
#'  QAP.MG(list(dv1), list(iv1), iv.names = c("intercept",paste0("IV",1:4)), samples = 3000)

#' @seealso \code{\link{QAP.MG}}
#' 
#' @references 
#' Dekker, D.; Krackhardt, D.; Snijders, T.A.B.  (2007).  \dQuote{Sensitivity of MRQAP Tests to Collinearity and Autocorrelation Conditions.}  \emph{Psychometrika}, 72(4), 563-581.
#' 
#' Krackhardt, D.  (1987).  \dQuote{QAP Partialling as a Test of Spuriousness.} \emph{Social Networks}, 9 171-186.
#' 
#' Krackhardt, D.  (1988).  \dQuote{Predicting With Networks: Nonparametric Multiple Regression Analyses of Dyadic Data.}  \emph{Social Networks}, 10, 359-382.
#'
#' @export
QAP <- function(dv, iv1,  iv.names, mode = "yQAP" ,samples = 1000, diag = F, directed = F){
  pb <- txtProgressBar(min = 0, max = samples, style = 3) # set progress bar
  
  if(!diag){
    diag(dv) <- NA
    for(IV in 1:length(iv1)){
      diag(iv1[[IV]]) <- NA
    }
    
  }
  
  
  matrixSubset <- matrix(T, nrow(dv), ncol(dv))
  if(!diag) diag(matrixSubset) <- F
  if(!directed) matrixSubset[upper.tri(matrixSubset)] <- F
  
  # get a vector for all variables
  
  observedEstimates <- lm(dv[matrixSubset] ~ sapply(iv1, function(x) x[matrixSubset]))$coefficients
  # prepare output file
  output <- data.frame(Estimates = observedEstimates)
  rownames(output) <- iv.names
  
  
  
  if(mode == "yQAP"){
    sampledEstimates <- data.frame()
    for(sampleNr in 1:samples){
      sampledEstimates <- rbind(sampledEstimates,lm(sample(dv[matrixSubset])
                                                    ~ sapply(iv1, function(x) x[matrixSubset]))$coefficients)
      setTxtProgressBar(pb, sampleNr)
    }
    
    for(est in 1:length(observedEstimates)){
      sampledEstimates[,est] <= observedEstimates[est]
      
      output[est,"p(1sided)"] <- ecdf(sampledEstimates[,est])(observedEstimates[est])
    }
    
  }
  
  if(mode == "dspQAP"){
    sampledEpsilon <- data.frame()
    observedEpsilons <- c()
    
    
    #add intercept as first element of IV
    iv.new <- list(matrix(1, nrow(dv), ncol(dv)))
    if(!diag) diag(iv.new[[1]]) <- NA
    for(ivNumber in 1:length(iv1)){
      iv.new[[ivNumber+1]] <- iv1[[ivNumber]]
    }
    
    for(sampleNr in 1:samples){
      for(xIV in 1:(length(iv1)+1)){
        
        
        # change variable names for easier comparison with the equations in Dekker et al. (2007)
        Y.m <- dv
        X.m <- iv.new[[xIV]]
        Z.m <- iv.new[-xIV]
        
        if(!directed){
          Y.m[upper.tri(Y.m, diag = !diag)] <- NA
          X.m[upper.tri(X.m, diag = !diag)] <- NA
          for(Zi in 1:length(Z.m)){
            Z.m[[Zi]][upper.tri(Z.m[[Zi]], diag = !diag)] <- NA
          }
        }
        deltas.m <- lm(as.numeric(X.m) ~ sapply(Z.m,function(x) x) -1)$coefficients # get the deltas (eq. 16 in paper)
        epsilon.m <- X.m - Reduce("+",lapply(1:length(Z.m), function(x) Z.m[[x]] * deltas.m[x])) # eq. 15
        
        if(!directed){
          epsilon.m[upper.tri(epsilon.m)] <- t(epsilon.m)[upper.tri(epsilon.m)]
        }
        
        # permutation
        require(sna)
        epsilon.perm <- rmperm(epsilon.m)
        
        perm.fit <- lm(as.numeric(Y.m) ~ as.numeric(epsilon.perm)  +  sapply(Z.m,function(x) x) -1)$coefficients # eq. 17
        sampledEpsilon[sampleNr, xIV] <- perm.fit[1]
        
        # compare sampledEpsilon with with epsilon!
        
        obs.fit <- lm(as.numeric(Y.m) ~ as.numeric(epsilon.m)  +  sapply(Z.m,function(x) x) -1)$coefficients
        observedEpsilons[xIV] <- obs.fit[1]
      }
    }
    
    for(est in 1:length(observedEpsilons)){
      sampledEpsilon[,est] <= observedEpsilons[est]
      
      output[est,"p(1sided)"] <- ecdf(sampledEpsilon[,est])(observedEpsilons[est])
    }
  } # end of dspQAP
  
  return(output)
}