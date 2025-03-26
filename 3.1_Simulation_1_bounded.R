###############################################
# Simulation 1_bounded baseline confounders + #
###############################################

source("./SubmitPrograms/0_Packages.R")
source("./SubmitPrograms/0_MyFunc.R")

# set the scenario
# =========================================
# effect: pos, neg, nul
psi.true = 0.7
psi.true = -0.7
psi.true = 0
# =========================================

# Number of cores to use
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# specify the seed value
.Random.seed = seed_summary$SIMS1_pos
  
result = foreach(i = 1:10, .combine = rbind,.packages = c("geeM", "tidyverse","brm","gesttools")) %dorng% {
    {
      n <- 1000
      id <- 1:n
      
      U <- runif(n, 0, 1); L1 <- runif(n, 0, 1 + U)
      A1 <- rbinom(n, 1, 1/(1 + exp(0.2 - 0.2*L1 )))
      L2 = rnorm(n, 0.5*L1 + 0.5*A1 + U,1)
      A2<-rbinom(n, 1, 1/(1 + exp(0.2 - 0.1*L1 - 0.2*L2 - 0.5*A1)))
      
      Z2 = Z1 = cbind(rep(1,n))
      W2 <- cbind(1, L1, A1, L2)
      W1 <- cbind(1, L1)
      
      L1.mis = runif(n, 0, 1 + U)
      L2.mis = rnorm(n, 0.5*L1.mis + 0.5 * A1 + U, 1)
      W2.mis <- cbind(1, L1.mis, A1, L2.mis)
      W1.mis <- cbind(1, L1.mis)
      
    X = cbind(1,L1,U)
    
    alpha.true = c(-1, 0.5, 0.5)
    p0p1.true = getProbRR(Z2 %*% psi.true,X %*% alpha.true)
    pA.true = p0p1.true[,1]
    pA.true[A2==1]=p0p1.true[A2==1,2]
    Y = rbinom(n, 1, pA.true)
    
  }
  
 res = matrix(NA,nrow = 1,ncol = 24)
 rescon = expand.grid(c("psi2","psi1","EY0","EY1","pred0","pred1"),
                    c(".PSmis",".OPmis",".BOTHmis",".NONEmis"))
 colnames(res) = paste0(rescon$Var1,rescon$Var2)
  
  ### PSmis
  pars2 = LSEst(Y,A2,Z2,W2)
  PS2.model.mis = glm(A2~W2.mis,family = binomial(link = "logit"))
  PS2.mis = predict(PS2.model.mis,type="response")
  psi2 = DREst(Y,A2,Z2,W2,pars2,PS2.mis)
  res[1,1] = psi2$minimum

  H0 = Y*exp(-psi2$minimum*A2)
  pars1 = LSEst(H0,A1,Z1,W1)
  PS1.model.mis = glm(A1~W1.mis,family = binomial(link = "logit"))
  PS1.mis = predict(PS1.model.mis,type="response")
  psi1 = DREst(H0,A1,Z1,W1,pars1,PS1.mis)
  res[1,2] = psi1$minimum
  
  EY0 = getProbRR(Z1 %*% psi1$minimum, W1 %*% pars1$alpha)[,1]
  EY1 = EY0 * exp(psi1$minimum + psi2$minimum)
  res[1,3] = mean(EY0)
  res[1,4] = mean(EY1)
  if(max(EY0) > 1) {res[1,5] = 1}else {res[1,5] = 0}
  if(max(EY1) > 1) {res[1,6] = 1}else {res[1,6] = 0}

  ### OPmis ###
  pars2.mis = LSEst(Y,A2,Z2,W2.mis)
  PS2.model = glm(A2~W2,family = binomial(link = "logit"))
  PS2 = predict(PS2.model,type="response")
  psi2 = DREst(Y,A2,Z2,W2.mis,pars2.mis,PS2)
  res[1,7] = psi2$minimum

  H0 = Y*exp(-psi2$minimum*A2)
  pars1.mis = LSEst(H0,A1,Z1,W1.mis)
  PS1.model = glm(A1~W1,family = binomial(link = "logit"))
  PS1 = predict(PS1.model,type="response")
  psi1 = DREst(H0,A1,Z1,W1.mis,pars1.mis,PS1)
  res[1,8] = psi1$minimum

  EY0 = getProbRR(Z1 %*% psi1$minimum, W1.mis %*% pars1.mis$alpha)[,1]
  EY1 = EY0 * exp(psi1$minimum + psi2$minimum)
  res[1,9] = mean(EY0)
  res[1,10] = mean(EY1)
  if(max(EY0) > 1) {res[1,11] = 1}else {res[1,11] = 0}
  if(max(EY1) > 1) {res[1,12] = 1}else {res[1,12] = 0}
# -----------------------------------------------------------------------------  
  
  #### BOTHmis ####
  psi2 = DREst(Y,A2,Z2,W2.mis,pars2.mis,PS2.mis)
  res[1,13] = psi2$minimum

  H0 = Y*exp(-psi2$minimum*A2)
  pars1.mis = LSEst(H0,A1,Z1,W1.mis)
  psi1 = DREst(H0,A1,Z1,W1.mis,pars1.mis,PS1.mis)
  res[1,14] = psi1$minimum
  
  EY0 = getProbRR(Z1 %*% psi1$minimum, W1.mis %*% pars1.mis$alpha)[,1]
  EY1 = EY0 * exp(psi1$minimum + psi2$minimum)
  res[1,15] = mean(EY0)
  res[1,16] = mean(EY1)
  if(max(EY0) > 1) {res[1,17] = 1}else {res[1,17] = 0}
  if(max(EY1) > 1) {res[1,18] = 1}else {res[1,18] = 0}
  
  #### BOTH:correct ####
  psi2 = DREst(Y,A2,Z2,W2,pars2,PS2)
  res[1,19] = psi2$minimum

  H0 = Y*exp(-psi2$minimum*A2)
  pars1 = LSEst(H0,A1,Z1,W1)
  psi1 = DREst(H0,A1,Z1,W1,pars1,PS1)
  res[1,20] = psi1$minimum

  EY0 = getProbRR(Z1 %*% psi1$minimum, W1 %*% pars1$alpha)[,1]
  EY1 = EY0 * exp(psi1$minimum + psi2$minimum)
  res[1,21] = mean(EY0)
  res[1,22] = mean(EY1)
  if(max(EY0) > 1) {res[1,23] = 1}else {res[1,23] = 0}
  if(max(EY1) > 1) {res[1,24] = 1}else {res[1,24] = 0}
  
  return(res)
}
