### testing the QAP.MG function ####
rm(list = ls())
library(sna)
# library(netglm)
devtools::load_all()
# create test data #
# inspired by the example function in sna::netlm
set.seed(56)
ivnet1<-sna::rgraph(20,4)
ivnet2<-sna::rgraph(20,4)

dv1<-ivnet1[1,,]+4*ivnet1[2,,]+2*ivnet1[3,,]   # Note that the fourth graph is unrelated
dv1 <- dv1 + rnorm(400,mean = 1, sd = 1)

dv2 <- 2*ivnet2[1,,]+3*ivnet2[2,,]+3*ivnet2[3,,]
dv2 <- dv2 + rnorm(400,mean = 1, sd = 1)
dvs <- list(dv1, dv2)

iv1 <- list(ivnet1[1,,],ivnet1[2,,],ivnet1[3,,], ivnet1[4,,])
iv2 <- list(ivnet2[1,,],ivnet2[2,,],ivnet2[3,,], ivnet2[4,,])
ivs <- list(iv1, iv2)
iv.names = c("intercept",paste0("IV",1:4))
###### y QAP ########
# both groups together
QAP.MG(dvs, ivs, iv.names = c("intercept",paste0("IV",1:4)), samples = 3000)

# group 1 separately
set.seed(56)
QAP.MG(list(dv1), list(iv1), iv.names = c("intercept",paste0("IV",1:4)), diag = F, samples = 100000)

# comparison with the netlm function from the sna-package
set.seed(56)
netlm(dv1, iv1, nullhyp = "qapy", diag = F, reps = 100000, test.statistic = "beta")

# computation time with parallel processing
system.time(QAP.MG(dvs, ivs, iv.names = c("intercept",paste0("IV",1:4)), samples = 10000))
system.time(QAP.MG(dvs, ivs, iv.names = c("intercept",paste0("IV",1:4)), samples = 10000, cpu = 4)) # test speed of 4 cpus

#### Dekker semi partialing ######
# group 1 separately
QAP.MG(list(dv1), list(iv1), mode = "dspQAP", iv.names = c("intercept",paste0("IV",1:4)), samples = 300)
# comparison with the netlm function from the sna-package
netlm(dv1, iv1, nullhyp = "qapspp", reps = 300)

# computation time with parallel processing
system.time(QAP.MG(dvs, ivs, mode = "dspQAP", iv.names = c("intercept",paste0("IV",1:4)), samples = 1000))
system.time(QAP.MG(dvs, ivs, mode = "dspQAP", iv.names = c("intercept",paste0("IV",1:4)), samples = 1000, cpu = 4)) # test speed of 4 cpus

