#' Multilevel MRQAP
#'
#' Estimates a MRQAP model taking the multilevel/grouped nature of the into account. 
#'
#'#' @param family family of the generalized linear model. default is "gaussian" for continuous dependent variables. F
#' or binary dependent variables "binomial" is advised. (currently only gaussian or binomial possible)
#' @export
QAP.ML <- function(formula, data, iv.list.per = "group", family = "gaussian",
                   iv.names = iv.names, mode = "yQAP" ,samples = 1000, diag = F, directed = T, 
                   cpu = 1, round.to = 5, logfilename = "QAP.log",
                   verbose = T,
                   global.deltas = T, return.perms = F){
  
  ##testing
  # formula <- as.formula(paste0(dv," ~ ",paste0(iv.names[-1], collapse = "+"),"+ (1|",nestVar,")"))
  # iv.list.per = "group"
  # diag = F
  # directed = T
  # family = "binomial"
  # cpu = 1
  # round.to = 5
  # global.deltas =T
  # samples = 100
  # return.perms =F
  # #mode = "dspQAP"
  # mode = "yQAP"
  # logfilename = "QAP.log"
  # verbose = T
  nestVar = "schools"
  dv = "network"
  
  require(lme4)
  
  # TODO: automatic check for directedness of DV matrix
  # TODO: check dimensions of IVs and DVs
  
  if(!diag){ # get rid of diagonal values
    if(verbose) cat("\n replace diagnoal values with NAs")
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
  
  ### make formula ####
  #TODO: integrate formula in function call
  
  
  
  
  # get the observed estimates
  if(verbose) cat("\n get the observed estimates")
  if(directed){ # directed
    
    if(family == "gaussian"){
      observedLm <- lmer(formula, data = data)
    }else{ # binomial 
      observedLm <- glmer(formula, data = data, family = binomial(logit))
    }
  }else{ #undirected
    stop("undirected must currently be handeled when creating the data object")
    
    
  }
  
  # prepare output file
  observedEstimates <- summary(observedLm)$coefficients[,1]
  output <- data.frame(Estimates = observedEstimates)
  #rownames(output) <- iv.names
  
  r.squared <- c("TODO") # TODO: implement r squared calculation / extraction for multilevel version
  
  if(cpu == 1) pb <- txtProgressBar(min = 0, max = samples, style = 3) # set progress bar
  #require(foreach)
  #require(doParallel)
  #cl <- makeCluster(cpu)
  #registerDoParallel(cl)
  
  
  ##### YQAP #######
  if(mode == "yQAP"){
    if(verbose) cat("\n estimating permuted networks with mode yQAP \n")
    sampledEstimates <- data.frame()
    for(sampleNr in 1:samples){
      #sampledEstimates <- foreach(sampleNr=1:samples, .combine=rbind) %dopar% { # for loop using parallel processing
      
      
      # permute dependent variable
      data.permuted <- data
      for(i in unique(data[,nestVar])){
        data.permuted[data.permuted[,nestVar] %in% i,dv] <- sample(data.permuted[data.permuted[,nestVar] %in% i,dv])
      }
      
      if(directed){
        if(family == "gaussian"){
          tmp <- lmer(formula, data = data.permuted)
        }else{ # binomial 
          tmp <- glmer(formula, data = data.permuted, family = binomial(logit))
        }
      }else{ # non-directed
        stop("undirected must currently be handeled when creating the data object")
      }
      
      #cat(paste0("\r", sampleNr," out of ", samples))
      if(cpu == 1) setTxtProgressBar(pb, sampleNr) # update progress bar
      #write.table(paste0(sampleNr," out of ", samples),file = logfilename)
      
      sampledEstimates <- rbind(sampledEstimates,summary(tmp)$coefficient[,1])
      
    }
    #stopCluster(cl)
    
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
      
      
    }
    
  }
  
  
  ##### Dekker semi partialing #######
  if(mode == "dspQAP"){
    stop("TODO: implement dspQAP for multilevel models")#TODO: implement for multilevel model 
  } # end of dspQAP
  
  output <- round(output,round.to)
  
  
  stars <- function(p.value){
    ifelse(p.value > 0.90 | p.value < 0.10,
           ifelse(p.value > 0.95 | p.value < 0.05, 
                  ifelse(p.value > 0.99 | p.value < 0.01, 
                         ifelse(p.value > 0.999 | p.value < 0.001, "***", "**"), "*"), "x"), "")
  }
  output[,"significance"] <- sapply(output$`p(1sided)`, stars)
  if(return.perms) return(list(sampledEstimates, observedEstimates))
  out <- list(mode = c(mode, samples), plots = ecdf.plots,output = output, r.squared = r.squared)
  class(out) <- "netglm"
  attr(out,"mode") <- mode
  attr(out,"model") <- "multilevel"
  return(out)
  
}