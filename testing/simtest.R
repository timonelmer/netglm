# Simulation to compare netglm against sna

# CONTENTS
# 0.0 Setup
# 1.0 Conditions
# 2.0 Simulation
# 3.0 Results

# 0.0 Setup ----
# Load required packages
library(sna)
devtools::load_all() # load current version of netglm
library(ggplot2)

# Set theme for ggplot
theme_set(theme_bw() +
            theme(panel.border     = element_blank(), 
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), 
                  axis.line        = element_line(colour = "black")))

# 1.0 Conditions ----

# Goal: Compare performance (running time) and parameter recovery of QAP and
# Dekker semi partialing as implemented in sna and netglm.

# In this test, we compare a network model with 4 independent variables. For
# simplicity, we fix the true beta parameters. We vary the number of nodes,
# the type of linear network (linear or logistic), and the type of permutation 
# method.

# Fixed factors for the simulation
n_iv     <- 4                     # Number of independent variables.
truebeta <- c(2.0, 2.5, 3.0, 0.0) # True coefficients for each IV.
truebeta <- c(-3.0, truebeta)     # Add true intercept. 
n_samp   <- 3000                  # Number of samples for permutations.

# Varying factors for the simulation
N_nodes <- c(5, 10, 20, 50)          # Number of nodes
lmtype  <- c("gaussian", "binomial") # Type of regression model
pmethod <- c("QAP", "DSP")           # Type of permutation

# Expand conditions
Cond <- expand.grid(N_nodes, lmtype, pmethod)
names(Cond) <- c("N_nodes", "lmtype", "pmethod")

# Prepare matrices to save output.
R <- 20 # Number of replications

outresults <- matrix(NA, nrow = nrow(Cond) * R, ncol = 8)
colnames(outresults) <- c("cond", "r", 
                          "time_sna", "bias_sna", 
                          "time_netglm", "bias_netglm", 
                          "time_netglm_par", "bias_netglm_par")

# Clear environment
rm(N_nodes, lmtype, pmethod)

# 2.0 Simulation ----

for (cond in 1:nrow(Cond)) {
  for (r in 1:R) {
    cat(sprintf("This is the %d-th replication of condition %d.\n", r, cond))
    
    # Set seed for this replication
    seed <- cond * 1000 + r 
    
    # Retrieve conditions for this replication
    n       <- Cond[cond, "N_nodes"]
    family  <- as.character(Cond[cond, "lmtype"])
    pmethod <- as.character(Cond[cond, "pmethod"])
    
    if (pmethod == "QAP") {pmode <- c("qapy", "yQAP")}
    if (pmethod == "DSP") {pmode <- c("qapspp", "dspQAP")}
    
    # Simulate data
    data_true <- gen.netglm(n = n, m = n_iv, 
                            family    = family, 
                            intercept = TRUE, 
                            beta      = truebeta,
                            seed      = seed)
    # Set diagonal to NA
    # for (i in 1:(n_iv + 1)) {
    #   diag(data_true$data[i, ,]) <- NA
    # }
    # rm(i)
    
    # Turn data into a list
    data_list <- asplit(data_true$data, 1)
    
    # Fit model in sna.
    if (family == "gaussian") {
      t0 <- proc.time()
      fit_sna <- netlm(data_list$y, 
                       data_list[1:n_iv], 
                       nullhyp = pmode[1],
                       reps    = n_samp)
      t_sna <- (proc.time() - t0)[3]
      
      bias_sna <- sum(coef(fit_sna) - truebeta)/(n_iv + 1)
    } else {
      t0 <- proc.time()
      fit_sna <- netlogit(data_list$y, 
                          data_list[1:n_iv], 
                          nullhyp = pmode[1],
                          reps    = n_samp)
      t_sna <- (proc.time() - t0)[3]
      
      bias_sna <- sum(coef(fit_sna) - truebeta)/(n_iv + 1)
    }
    
    # Fit model in netglm.
    t0 <- proc.time()
    fit_netglm <- QAP.MG(list(data_list$y),
                         list(data_list[1:n_iv]),
                         iv.names = c("int", paste0("iv", 1:n_iv)),
                         family   = family,
                         mode     = pmode[2],
                         samples  = n_samp) 
    t_netglm <- (proc.time() - t0)[3]
    
    bias_netglm <- sum(fit_netglm$output$Estimates - truebeta)/(n_iv + 1)
    
    # Fit model in netglm using parallel computing.
    t0 <- proc.time()
    fit_netpar <- QAP.MG(list(data_list$y),
                         list(data_list[1:n_iv]),
                         iv.names = c("int", paste0("iv", 1:n_iv)),
                         family   = family,
                         mode     = pmode[2],
                         samples  = n_samp,
                         cpu      = 4) 
    t_netpar <- (proc.time() - t0)[3]
    
    bias_netpar <- sum(fit_netglm$output$Estimates - truebeta)/(n_iv + 1)
    
    # Save output in output matrix.
    outresults[((cond - 1) * R) + r,] <- c(cond, r, t_sna, bias_sna, t_netglm, 
                                           bias_netglm, t_netpar, bias_netpar)
    
    rm(data_true, data_list, seed, n, family, pmethod, pmode, 
       fit_sna, t_sna, bias_sna, fit_netglm, t_netglm,
       bias_netglm, fit_netpar, t_netpar, bias_netpar)
    
  }
}

# 3.0 Results ----

outresults <- as.data.frame(outresults)

outsum <- sapply(outresults[, -(1:2)], function(x) {
  tapply(x, outresults$cond, mean)
  })

#### END #### 