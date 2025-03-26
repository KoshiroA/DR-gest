##############################################################
# Additional Simulation S2: evaluation of doubly robusteness #
##############################################################
source("./SubmitPrograms/0_Packages.R")
source("./SubmitPrograms/0_MyFunc.R")

rm(list = ls()[!str_detect(ls(), "result") & !str_detect(ls(), "seed") & 
                 !str_detect(ls(), "Est") & !str_detect(ls(), "psi1.true_")])

library(doParallel)
library(doRNG)
library(foreach)

# Number of cores to use
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

#====================================================#
# Evaluation under scenario under Wang et al. (2023) #
#====================================================#

# scenarios #
  # use functions for data generation developed by Wang et al. (2023)
  # specify these settings
#---------------------------------------------------------------------# 
source("2.0.0 MyFunc.R")
source("2.0.1 GetProb.R")
source("2.0.3 GetContrast.R")
source("2.1 DataGen.R")

{
  ee.assume.true.alpha = T; 
  
  # Nuisance model for y.hat.type for DR estimation; 
  y.hat.type = "MLE" # "saturated" or "MLE" or "logit" or "diy"
  prop.score.type = "correct" # "correct" or "incorrect"
  vectorize.get.prob = F;parallel.root = F
  l1.type ="correct" # "correct" or "incorrect" or "diy"
  gop.type = "correct" # "correct" or "incorrect"
}

# true parameter given by Wamg et al.
{
  alpha.true  = t(matrix(rep(c(0,0.7),5),2,5))
  beta.3 =  c(-0.5, 1)
  beta.1 = beta.2 = beta.4 = beta.5 =c(-0.5, 0.1);
  beta.true   = t(matrix(c(beta.1, beta.2, beta.3, beta.4, beta.5),2,5))
  gamma.true  = t(matrix(c(0.1,-0.5,NA,NA,0.1,-0.5,0.1,-0.5),4,2))
  
  params.true = list(alpha.true=alpha.true, beta.true=beta.true,
                     gamma.true=gamma.true)
  rownames(gamma.true) = c("a0","a1")
  
}
#---------------------------------------------------------------------# # discrete l0
# discrete l0
l0.type = "discrete"

# continuous l0
l0.type = "continuous"

# discrete and incorrect l0
l0.type = "discrete"
#---------------------------------------------------------------------# 

seed_summary$

result = foreach(i = 1:500, .combine = rbind,.packages = c("geeM", "tidyverse","brm","gesttools")) %dorng% {
  {
    datagen = DataGen(1000)
    L1 = datagen$l0.2.vec
    A1 = datagen$a0.vec
    # L2 = datagen$l1.vec # correct
    L2 = datagen$y.cov.tilde[,4] # incorrect
    A2 = datagen$a1.vec
    Z2 = W2 = cbind(1, L1, A1, L2, A1*L2)
    W1 = Z1 = cbind(1, L1)
    Y = datagen$y.vec
  }
  
  res = matrix(NA,nrow = 1,ncol = 11)
  colnames(res) = c(paste0("psi2.",1:5),paste0("psi1.",1:2),"EY0","EY1","EY0>1","EY1>1")
  
  ### 2nd stage ###
  pars2 = LSEst(Y,A2,Z2,W2)
  PS2.model = glm(A2~W2,family = binomial(link = "logit"))
  PS2 = predict(PS2.model,type="response")
  psi2 = DREst(Y,A2,Z2,W2,pars2,PS2)
  res[1,1:5] = psi2$par
  
  ### 1st stage ###
  H1 = Y*exp(-Z2 %*% psi2$par *A2)
  pars1 = LSEst(H1,A1,Z1,W1)
  PS1.model = glm(A1~W1,family = binomial(link = "logit"))
  PS1 = predict(PS1.model,type="response")
  psi1 = DREst(H1,A1,Z1,W1,pars1,PS1)
  res[1,6:7] = psi1$par
  
  EY0 = getProbRR(Z1 %*% psi1$par, W1 %*% pars1$alpha)[,1]
  EY1 = EY0 * exp(Z1 %*% psi1$par + Z2 %*% psi2$par)
  res[1,8] = mean(EY0)
  res[1,9] = mean(EY1)
  
  if(max(EY0) > 1) {res[1,10] = 1}else {res[1,10] = 0}
  if(max(EY1) > 1) {res[1,11] = 1}else {res[1,11] = 0}
  
  return(res)
}