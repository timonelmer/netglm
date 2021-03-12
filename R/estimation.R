##### netglm ######
#### estimation ###



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


# mutligroup QAP
QAP.MG <- function(dvs, ivs, iv.list.per = "iv", family = "gaussian",
                   iv.names = iv.names, mode = "yQAP" ,samples = 1000, diag = F, directed = T, 
                   global.deltas = T, ecdf.plot = F, return.perms = F){
  
  
  
  if(!diag){ # get rid of diagonal values
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
    cat("\n rearrange list of ivs so its iv list per group \n")
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
  pb <- txtProgressBar(min = 0, max = samples, style = 3) # set progress bar
  
  if(directed){
    
    observedLm <- glm(unlist(dvs) ~ Reduce(rbind,lapply(1:length(ivs),
                                                        function(x) sapply(ivs[[x]], function(y) y)
    )), family = family)
    
  }else{
    if(!(family == "gaussian"))stop("non-gaussian estimation not yet implemented for undirected networks")
    
    observedLm <- glm(unlist(lapply(dvs, function(x) x[lower.tri(x, diag = diag)])) ~
                        Reduce(rbind,lapply(1:length(ivs),function(x) sapply(ivs[[x]], function(y) y[lower.tri(y, diag = diag)]))),
                      family = family)
    
    
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
  
  
  ##### YQAP #######
  if(mode == "yQAP"){
    sampledEstimates <- data.frame()
    for(sampleNr in 1:samples){
      if(directed){
        
        sampledEstimates <- rbind(sampledEstimates, glm(unlist(lapply(dvs, function(x) sample(x))) ~ Reduce(rbind,lapply(1:length(ivs),
                                                                                                                         function(x) sapply(ivs[[x]], function(y) y)
        )), family = family)$coefficients)
      }else{
        sampledEstimates <- rbind(sampledEstimates, glm(unlist(lapply(dvs, function(x) sample(x[lower.tri(x, diag = diag)]))) ~ 
                                                          Reduce(rbind,lapply(1:length(ivs), function(x) sapply(ivs[[x]], function(y) y[lower.tri(y, diag = diag)])
                                                          )), family = family)$coefficients)
      }
      
      setTxtProgressBar(pb, sampleNr)
    }
    
    ecdf.plots <- list()
    for(est in 1:length(observedEstimates)){
      
      output[est,"p(1sided)"] <- ecdf(sampledEstimates[,est])(observedEstimates[est])
      output[est,"abs(p)"] <- output[est,"p(1sided)"] 
      output[est,"abs(p)"][output[est,"p(1sided)"] > 0.5] <- 1-output[est,"abs(p)"][output[est,"p(1sided)"] > 0.5]
      output[est,"adj.d"] <- abs(mean(sampledEstimates[,est]) - observedEstimates[est])/sd(sampledEstimates[,est])
      output[est,"Exp.V"] <- mean(sampledEstimates[,est], na.rm = T)
      output[est,"Exp.V.sd"] <- sd(sampledEstimates[,est], na.rm = T)
      output[est,"2.5th P"] <- quantile(sampledEstimates[,est],probs=c(.025))
      output[est,"97.5th P"] <- quantile(sampledEstimates[,est],probs=c(.975))
      
      
      if(ecdf.plot){
        
        require(ggplot2)
        
        ecdf.plots[[est]] <- plot(ecdf(sampledEstimates[,est])(observedEstimates[est]))
        
      }
    }
    
  }
  
  
  ##### Dekker semi partialing #######
  if(mode == "dspQAP"){
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
    
    
    for(sampleNr in 1:samples){
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
            #global deltas
            X.m.global <- list()
            Z.m.global <- list()
            for(obs in 1:length(dvs)){
              X.m.global[[obs]] <- ivs.new[[obs]][[xIV]]
              Z.m.global[[obs]] <- ivs.new[[obs]][-xIV]
            }
            
            deltas.m <- lm(unlist(X.m.global) ~ Reduce(rbind,lapply(1:length(Z.m.global),
                                                                    function(x) sapply(Z.m.global[[x]], function(y) y))) -1)$coefficients # get the deltas (eq. 16 in paper)
            
            
            
          }else{
            #local deltas
            deltas.m <- lm(as.numeric(X.m) ~ sapply(Z.m,function(x) x) -1)$coefficients # get the deltas (eq. 16 in paper)
            if(family != "gaussian") stop("to implement for non gaussian")
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
        
        perm.fit <- glm(unlist(dvs) ~ 
                          unlist(epsilon.perm.list)  + 
                          Reduce(rbind,lapply(1:length(Z.m.list),function(x) sapply(Z.m.list[[x]], function(y) y)
                          )) -1, family = family)$coefficients # eq. 17
        
        
        
        
        sampledEpsilon[sampleNr, xIV] <- perm.fit[1]
        
        # compare sampledEpsilon with with epsilon!
        
        obs.fit <- glm(unlist(dvs) ~ 
                         unlist(epsilon.list)  + 
                         Reduce(rbind,lapply(1:length(Z.m.list),function(x) sapply(Z.m.list[[x]], function(y) y)
                         )) -1, family = family)$coefficients # eq. 17
        observedEpsilons[xIV] <- obs.fit[1]
      }
      setTxtProgressBar(pb, sampleNr) # update progress bar
      
    }
    
    ecdf.plots <- list()
    for(est in 1:length(observedEpsilons)){
      sampledEpsilon[,est] <= observedEpsilons[est]
      
      output[est,"p(1sided)"] <- ecdf(sampledEpsilon[,est])(observedEpsilons[est])
      #output[est,"p(Percentile,1sided)"] <- (sum(sampledEpsilon[,est] <= observedEpsilons[est])/nrow(sampledEpsilon))
      
      if(ecdf.plot){
        
        require(ggplot2)
        # g <- qplot(sampledEpsilon[,est], bins = 100) + 
        #   geom_vline(xintercept = observedEpsilons[est]) + 
        #   annotate("text", label = paste0(round(output[est,"p(1sided)"],3),"%"), x = observedEpsilons[est], y = 50) + xlab(iv.names[est])
        # ecdf.plots[[est]] <- g
        ecdf.plots[[est]] <- plot(ecdf(sampledEpsilon[,est])(observedEpsilons[est]))
        #rm(g)
      }
    }
  } # end of dspQAP
  
  
  output <- round(output,5)
  
  stars <- function(p.value){
    ifelse(p.value > 0.90 | p.value < 0.10,
           ifelse(p.value > 0.95 | p.value < 0.05, 
                  ifelse(p.value > 0.99 | p.value < 0.01, 
                         ifelse(p.value > 0.999 | p.value < 0.001, "***", "**"), "*"), "x"), "")
  }
  output[,"significance"] <- sapply(output$`p(1sided)`, stars)
  if(return.perms) return(list(sampledEstimates, observedEstimates))
  return(list(mode = c(mode, samples), plots = ecdf.plots,output = output, r.squared = r.squared))
  
}


