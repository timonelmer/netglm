# mutligroup QAP

#' Multigroup MRQAP
#'
#' Estimates a MRQAP model taking the multilevel/grouped nature of the into account. 
#' An application and detailed description of the multigroup extension 
#' of MRQAP can be found in Elmer and Stadtfeld (2020).
#'
#' @param dvs a list of matrices where each matrix represents the dependent network of one group. Cells of the matrix indicating the presence of a tie (1 = tie, 0 = no tie for binary variables) 
#'  or the weight of a tie (for continuous tie variables) characterizing the dependent networks
#' @param ivs a list of lists with where each sets of independent matrices are grouped by group and then independent variables. 
#' If the list is organized per independent network the argument iv.list.per = "iv" can be used to restructure the data.
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
#' ivnet2<-sna::rgraph(20,4)
#'
#' dv1<-ivnet1[1,,]+4*ivnet1[2,,]+2*ivnet1[3,,]   # Note that the fourth graph is unrelated
#' dv1 <- dv1 + rnorm(400,mean = 1, sd = 1)
#' 
#' dv2 <- 2*ivnet2[1,,]+3*ivnet2[2,,]+3*ivnet2[3,,]
#' dv2 <- dv2 + rnorm(400,mean = 1, sd = 1)
#' dvs <- list(dv1, dv2)

#' iv1 <- list(ivnet1[1,,],ivnet1[2,,],ivnet1[3,,], ivnet1[4,,])
#' iv2 <- list(ivnet2[1,,],ivnet2[2,,],ivnet2[3,,], ivnet2[4,,])
#' ivs <- list(iv1, iv2)
#' iv.names = c("intercept",paste0("IV",1:4))
#'  QAP.MG(dvs, ivs, iv.names = c("intercept",paste0("IV",1:4)), samples = 3000)
#' @seealso \code{\link{QAP}}
#'
#' @references 
#' Dekker, D., Krackhardt, D., & Snijders, T.A.B.  (2007).  \dQuote{Sensitivity of MRQAP Tests to Collinearity and Autocorrelation Conditions.}  \emph{Psychometrika}, 72(4), 563-581.
#' 
#' Elmer, T., & Stadtfeld, C. (2020). \dQuote{Depressive symptoms are associated with social isolation in face-to-face interaction networks}. \emph{Scientific Reports}, 1–12. https://doi.org/10.1038/s41598-020-58297-9
#' 
#' Krackhardt, D.  (1987).  \dQuote{QAP Partialling as a Test of Spuriousness.} \emph{Social Networks}, 9 171-186.
#' 
#' Krackhardt, D.  (1988).  \dQuote{Predicting With Networks: Nonparametric Multiple Regression Analyses of Dyadic Data.}  \emph{Social Networks}, 10, 359-382.
#'
#'
#' @export
QAP.MG <- function(dvs, ivs, iv.list.per = "group", family = "gaussian",
                   iv.names = iv.names, mode = "yQAP" ,samples = 1000, 
                   diag = FALSE, directed = TRUE, cpu = 1, round.to = 5, 
                   logfilename = "QAP.log", verbose = TRUE,
                   global.deltas = FALSE){
  
  # #testing
  # ivnet1<-sna::rgraph(20,4)
  # ivnet2<-sna::rgraph(20,4)
  # 
  # dv1<-ivnet1[1,,]+4*ivnet1[2,,]+2*ivnet1[3,,]   # Note that the fourth graph is unrelated
  # dv1 <- dv1 + rnorm(400,mean = 1, sd = 1)
  # 
  # dv2 <- 2*ivnet2[1,,]+3*ivnet2[2,,]+3*ivnet2[3,,]
  # dv2 <- dv2 + rnorm(400,mean = 1, sd = 1)
  # dvs <- list(dv1, dv2)
  # 
  # iv1 <- list(ivnet1[1,,],ivnet1[2,,],ivnet1[3,,], ivnet1[4,,])
  # iv2 <- list(ivnet2[1,,],ivnet2[2,,],ivnet2[3,,], ivnet2[4,,])
  # ivs <- list(iv1, iv2)
  # iv.names <- 0:4
  # iv.list.per = "group"
  # diag = F
  # directed = T
  # family = "gaussian"
  # cpu = 1
  # round.to = 5
  # global.deltas =T
  # samples = 100
  # return.perms =F
  # #mode = "dspQAP"
  # mode = "yQAP"
  # logfilename = "QAP.log"
  # verbose = T
  
  # TODO: automatic check for directedness of DV matrix
  # TODO: check dimensions of IVs and DVs
  
  if(!diag){ # get rid of diagonal values
    if(verbose) cat("\n replace diagonal values with NAs")
    for(DV in 1:length(dvs)){
      diag(dvs[[DV]]) <- NA
    }
    
    for(IV in 1:length(ivs)){  
      IV.dat <- ivs[[IV]]
      for(subIV in 1:length(ivs[[IV]])){
        diag(ivs[[IV]][[subIV]]) <- NA
      }
    }
  }
  
  #rearrange list of ivs so its iv list per group 
  if(iv.list.per == "iv"){
    if(verbose) cat("\n rearrange list of ivs so its iv list per group \n")
    ivs.new <- list()
    n.ivs <- length(ivs)
    n.groups <- length(ivs[[1]])
    
    for(group in 1:n.groups){
      tmp.list <- list()
      for(iv in 1:n.ivs){
        tmp.list[[iv]] <- ivs[[iv]][[group]]
      }
      ivs.new[[group]] <- tmp.list
    }
    ivs <- ivs.new
  }
  
  
  # get the observed estimates
  if(verbose) cat("\n get the observed estimates")
  if(directed){ # directed
    
    if(family == "gaussian"){
      observedLm <- lm(unlist(dvs) ~ Reduce(rbind,lapply(1:length(ivs),
                                                         function(x) sapply(ivs[[x]], function(y) y)
      )))
    }else if(family == "inverse-Gaussian"){  
      observedLm <- glm(unlist(dvs) ~ Reduce(rbind,lapply(1:length(ivs),
                                                          function(x) sapply(ivs[[x]], function(y) y)
      )), family = tweedie(var.power = 3, link.power = 0))  
    }else{ # non-gaussian
      observedLm <- glm(unlist(dvs) ~ Reduce(rbind,lapply(1:length(ivs),
                                                          function(x) sapply(ivs[[x]], function(y) y)
      )), family = family)
    }
  }else{ #undirected
    if(family == "gaussian"){
      observedLm <- lm(unlist(lapply(dvs, function(x) x[lower.tri(x, diag = diag)])) ~
                         Reduce(rbind,lapply(1:length(ivs),function(x) sapply(ivs[[x]], function(y) y[lower.tri(y, diag = diag)]))),
      )
    }else if(family == "inverse-Gaussian"){  
      observedLm <- glm(unlist(lapply(dvs, function(x) x[lower.tri(x, diag = diag)])) ~
                          Reduce(rbind,lapply(1:length(ivs),function(x) sapply(ivs[[x]], function(y) y[lower.tri(y, diag = diag)]))),
                        family = tweedie(var.power = 3, link.power = 0))  
    }else{# non-gaussian
      observedLm <- glm(unlist(lapply(dvs, function(x) x[lower.tri(x, diag = diag)])) ~
                          Reduce(rbind,lapply(1:length(ivs),function(x) sapply(ivs[[x]], function(y) y[lower.tri(y, diag = diag)]))),
                        family = family)
    }
    
    
  }
  
  
  if(mode == "linearregression"){
    names(observedLm$coefficients) <- iv.names
    return(observedLm) 
  }
  
  # prepare output file
  observedEstimates <- observedLm$coefficients
  output <- data.frame(Estimates = observedEstimates)
  rownames(output) <- iv.names
  
  r.squared <- c(summary(observedLm)$r.squared, summary(observedLm)$adj.r.squared)
  if(is.null(r.squared)) r.squared <- c("TODO for non-gaussian",NA)
  names(r.squared) <- c("r.squared","adj.r.squared")
  if(family == "binomial") r.squared <- DescTools::PseudoR2(observedLm)
  
  
  #if(cpu == 1) pb <- txtProgressBar(min = 0, max = samples, style = 3) # set progress bar
  require(foreach)
  require(doParallel)
  cl <- makeCluster(cpu)
  registerDoParallel(cl)
  
  
  ##### YQAP #######
  
  if(!directed){ dvs <- lapply(dvs, function(x) {
      x[lower.tri(x, diag = FALSE)] <- NA
      return(x)
    })
  }
  
  if(mode == "yQAP"){
    if(verbose) cat("\n estimating permuted networks with mode yQAP \n")
    sampledEstimates <- data.frame()
    #for(sampleNr in 1:samples){
    sampledEstimates <- foreach(sampleNr=1:samples, .combine=rbind) %dopar% { # for loop using parallel processing
      if(family=="gaussian"){
        tmp <- lm(unlist(lapply(dvs, function(x)
        {
          s <- sample(1:ncol(x))
          return(c(x[s, s]))
        })) ~
          Reduce(rbind, lapply(1:length(ivs), function(x)
            sapply(ivs[[x]], function(y)
              y))))$coefficients
        
      }else if(family == "inverse-Gaussian"){  
        tmp <- glm(unlist(lapply(dvs, function(x)
        {
          s <- sample(1:ncol(x))
          return(c(x[s, s]))
        })) ~
          Reduce(rbind, lapply(1:length(ivs), function(x)
            sapply(ivs[[x]], function(y)
              y))), family = tweedie(var.power = 3, link.power = 0))$coefficients  
        
      }else{#non-gaussian
        tmp <- glm(unlist(lapply(dvs, function(x)
        {
          s <- sample(1:ncol(x))
          return(c(x[s, s]))
        })) ~
          Reduce(rbind, lapply(1:length(ivs), function(x)
            sapply(ivs[[x]], function(y)
              y))), family = family)$coefficients
      }
      
      print(tmp)
      #cat(paste0("\r", sampleNr," out of ", samples))
      #if(cpu == 1) setTxtProgressBar(pb, sampleNr) # update progress bar
      #write.table(paste0(sampleNr," out of ", samples),file = logfilename)
    }
    stopCluster(cl)
    
    ecdf.plots <- list()
    for(est in 1:length(observedEstimates)){
      
      #output[est,"p(1sided)"] <- mean(sampledEstimates[,est] <= observedEstimates[est])# ecdf(sampledEstimates[,est])(observedEstimates[est])
      #output[est,"p(2sided)"] <- mean(sampledEstimates[,est] >= abs(observedEstimates[est]))
      
      output[est,"Pr(<=b)"] <- mean(sampledEstimates[,est] <= observedEstimates[est])
      output[est,"Pr(>=b)"] <- mean(sampledEstimates[,est] >= observedEstimates[est])
      
      output[est,"Pr(>=|b|)"] <- mean(abs(sampledEstimates[,est]) >= abs(observedEstimates[est])) # because see https://stats.stackexchange.com/questions/35000/why-do-we-take-the-absolute-value-in-a-hypothesis-test
      
      #output[est,"adj.d"] <- abs(mean(sampledEstimates[,est]) - observedEstimates[est])/sd(sampledEstimates[,est])
      #output[est,"Exp.V"] <- mean(sampledEstimates[,est], na.rm = T)
      #output[est,"Exp.V.sd"] <- sd(sampledEstimates[,est], na.rm = T)
      output[est,"2.5th P"] <- quantile(sampledEstimates[,est],probs=c(.025))
      output[est,"97.5th P"] <- quantile(sampledEstimates[,est],probs=c(.975))
      

      
      
    }
    
  }
  
  
  ##### Dekker semi partialing #######
  if(mode == "dspQAP"){
    if(verbose) cat("\n estimating permuted networks with mode Dekker semi partialing (dspQAP) \n")
    sampledEpsilon <- data.frame()
    observedEpsilons <- c()
    
    #add intercept as first element of each IV
    ivs.new <- list()
    for(IV in 1:length(ivs)){
      ivs.new[[IV]] <- list(matrix(1, nrow(dvs[[IV]]), nrow(dvs[[IV]])))
      for(ivNumber in 1:length(ivs[[IV]])){
        ivs.new[[IV]][[ivNumber+1]] <- ivs[[IV]][[ivNumber]]
      }
    }
    
    
    #for(sampleNr in 1:samples){
    Epsilon.list <- foreach(sampleNr=1:samples, .combine=rbind,
                            .multicombine = T, .init = list()) %dopar% { # for loop using parallel processing
                              for(xIV in 1:(length(ivs[[1]])+1)){ # for each independent variable
                                
                                # storage places
                                epsilon.perm.list <- list() # object to store the permuted epsilon matrices
                                epsilon.list <- list() # object to store the observed epsilon matrices
                                Z.m.list <- list() # object to store all the list of Z matrices
                                
                                for(set in 1:length(ivs)){
                                  # change variable names for easier comparison with the equations in Dekker et al. (2007)
                                  Y.m <- dvs[[set]]
                                  X.m <- ivs.new[[set]][[xIV]]
                                  Z.m <- ivs.new[[set]][-xIV]
                                  
                                  if(!directed){
                                    Y.m[upper.tri(Y.m, diag = !diag)] <- NA
                                    X.m[upper.tri(X.m, diag = !diag)] <- NA
                                    for(Zi in 1:length(Z.m)){
                                      Z.m[[Zi]][upper.tri(Z.m[[Zi]], diag = !diag)] <- NA
                                    }
                                  }
                                  
                                  
                                  if(global.deltas){ # default
                                    # TODO: ask Tom. If in eq. 16 the deltas (associations between Z and X) must 
                                    # be estimate within each group (I think so, thus global.deltas = F) or for all groups.
                                    
                                    #global deltas
                                    X.m.global <- list()
                                    Z.m.global <- list()
                                    for(obs in 1:length(dvs)){
                                      X.m.global[[obs]] <- ivs.new[[obs]][[xIV]]
                                      Z.m.global[[obs]] <- ivs.new[[obs]][-xIV]
                                    }
                                    
                                    deltas.m <- lm(unlist(X.m.global) ~ Reduce(rbind,lapply(1:length(Z.m.global),
                                                                                            function(x) sapply(Z.m.global[[x]], function(y) y))) -1)$coefficients # get the deltas (eq. 16 in paper)
                                    
                                    
                                    if(family != "gaussian") stop("to implement for non gaussian")
                                  }else{
                                    #local deltas
                                    deltas.m <- lm(as.numeric(X.m) ~ sapply(Z.m,function(x) x) -1)$coefficients # get the deltas (eq. 16 in paper)
                                    #if(family != "gaussian") stop("to implement for non gaussian")
                                    
                                  }
                                  
                                  epsilon.m <- X.m - Reduce("+",lapply(1:length(Z.m), function(x) Z.m[[x]] * deltas.m[x])) # eq. 15
                                  
                                  if(!directed){
                                    epsilon.m[upper.tri(epsilon.m)] <- t(epsilon.m)[upper.tri(epsilon.m)]
                                  }
                                  
                                  # permutation
                                  require(sna)
                                  epsilon.m.perm <- rmperm(epsilon.m)
                                  epsilon.m.perm[upper.tri(epsilon.m.perm, diag = !diag)] <- NA
                                  epsilon.perm.list[[set]] <- epsilon.m.perm 
                                  epsilon.m[upper.tri(epsilon.m, diag = !diag)] <- NA
                                  epsilon.list[[set]] <- epsilon.m
                                  Z.m.list[[set]] <- Z.m
                                  
                                }
                                
                                if(family == "gaussian"){
                                  perm.fit <- lm(unlist(dvs) ~ 
                                                   unlist(epsilon.perm.list)  + 
                                                   Reduce(rbind,lapply(1:length(Z.m.list),function(x) sapply(Z.m.list[[x]], function(y) y)
                                                   )) -1)$coefficients # eq. 17
                                }else{#non-gaussian
                                  perm.fit <- glm(unlist(dvs) ~ 
                                                    unlist(epsilon.perm.list)  + 
                                                    Reduce(rbind,lapply(1:length(Z.m.list),function(x) sapply(Z.m.list[[x]], function(y) y)
                                                    )) -1, family = family)$coefficients # eq. 17
                                }
                                
                                
                                
                                
                                sampledEpsilon[sampleNr, xIV] <- perm.fit[1]
                                if(sampleNr %in% 1:2){ # only compute the observedEpsilons once (and compare with second sample for safety)
                                  # compare sampledEpsilon with epsilon!
                                  if(family == "gaussian"){
                                    obs.fit <- lm(unlist(dvs) ~ 
                                                    unlist(epsilon.list)  + 
                                                    Reduce(rbind,lapply(1:length(Z.m.list),function(x) sapply(Z.m.list[[x]], function(y) y)
                                                    )) -1)$coefficients # eq. 17
                                  }else{#non-gaussian
                                    obs.fit <- glm(unlist(dvs) ~ 
                                                     unlist(epsilon.list)  + 
                                                     Reduce(rbind,lapply(1:length(Z.m.list),function(x) sapply(Z.m.list[[x]], function(y) y)
                                                     )) -1, family = family)$coefficients # eq. 17
                                  }
                                  observedEpsilons[xIV] <- obs.fit[1]
                                }
                              }
                              print(list(sampledEpsilon[sampleNr,],observedEpsilons))
                              #setTxtProgressBar(pb, sampleNr) # update progress bar
                              #write.table(paste0(sampleNr," out of ", samples),file = logfilename)
                              
                            }
    stopCluster(cl)
    
    sampledEpsilon <- data.frame(matrix(unlist(Epsilon.list[1:samples]), nrow=length(Epsilon.list[1:samples]), byrow=TRUE))
    if(!all(Epsilon.list[samples+1][[1]] == Epsilon.list[samples+2][[1]])) stop("observed Epsilons are not identical per sample")
    observedEpsilons <- Epsilon.list[samples+1][[1]]
    
    ecdf.plots <- list()
    for(est in 1:length(observedEpsilons)){
      sampledEpsilon[,est] <= observedEpsilons[est] 
      
      output[est,"p(1sided)"] <- ecdf(sampledEpsilon[,est])(observedEpsilons[est])
      #output[est,"p(Percentile,1sided)"] <- (sum(sampledEpsilon[,est] <= observedEpsilons[est])/nrow(sampledEpsilon))    }
      output[est,"abs(p)"] <- output[est,"p(1sided)"] 
      output[est,"abs(p)"][output[est,"p(1sided)"] > 0.5] <- 1-output[est,"abs(p)"][output[est,"p(1sided)"] > 0.5]
    }
  } # end of dspQAP
  
  output <- round(output,round.to)
  
  # TODO: how to report the sampledEpsilons
  if(mode == "dspQAP"){sampledEstimates <- NULL}else{
    colnames(sampledEstimates) <- iv.names
    rownames(sampledEstimates) <- 1:nrow(sampledEstimates)
    output[,"significance"] <- sapply(output$`p(1sided)`, stars)
  }
  if(mode == "yQAP"){sampledEpsilon <- NULL}else{
    colnames(sampledEpsilon) <- iv.names
    rownames(sampledEpsilon) <- 1:nrow(sampledEpsilon)
  }
  
  out <- list(mode = c(mode, samples), sampledEpsilon = sampledEpsilon, sampledEstimates = sampledEstimates, output = output, r.squared = r.squared)
  class(out) <- "netglm"
  return(out)
  
}
