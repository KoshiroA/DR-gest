##################################################
# Simulation 2: evaluation in multiple time points
##################################################

source("./SubmitPrograms/0_Packages.R")
source("./SubmitPrograms/0_MyFunc.R")

# set the scenario
# =========================================
# effect: pos, neg, nul
#psi = rbind(0.2,0.3)
psi = rbind(-0.2,-0.3)
psi = rbind(0,0)

# prev: rare, low, common, middle, high
const = -1
const = -0.5
const = 0
const = 0.5
const = 1
# =========================================

# =============== #
# Proposed method #
# =============== #

.Random.seed = seed_summary$SIM2_neg1

result = 
  foreach(i = 1:5, .combine = rbind,.packages = c("geeM", "tidyverse","brm","gesttools")) %dorng% {
    {
      n = 1000
      id = 1:n
      U = rnorm(n, 0, 1); L1 = rnorm(n, U, 1)
      A1 = rbinom(n, 1, 1/(1 + exp(0.2 - 0.2*L1)))
      L2 = rnorm(n, 0.5*L1 + 0.5*A1 + U,1)
      A2 = rbinom(n, 1, 1/(1 + exp(0.2 - 0.1*L1 - 0.2*L2 - 0.5*A1)))
      L3 = rnorm(n, 0.5*L2 + 0.5*A2 + U,1)
      A3 = rbinom(n, 1, 1/(1 + exp(0.2 - 0.1*L1 - 0.2*L3 - 0.5*A2)))
      L4 = rnorm(n, 0.5*L3 + 0.5*A3 + U,1)
      A4 = rbinom(n, 1, 1/(1 + exp(0.2 - 0.1*L1 - 0.2*L4 - 0.5*A3)))
      
      Z4  = Z3 = Z2 = Z1 = cbind(rep(1, n))
      W4 = cbind(1, L1, A3, L4)
      W3 = cbind(1, L1, A2, L3)
      W2 = cbind(1, L1, A1, L2)
      W1 = cbind(1, L1)
      
      L1.mis = rnorm(n, U, 1)
      L2.mis = rnorm(n, 0.5*L1.mis + 0.5*A1 + U,1)
      L3.mis = rnorm(n, 0.5*L2.mis + 0.5*A2 + U,1)
      L4.mis = rnorm(n, 0.5*L3.mis + 0.5*A3 + U,1)
      
      W4.mis = cbind(1, L1.mis, A3, L4.mis)
      W3.mis = cbind(1, L1.mis, A2, L3.mis)
      W2.mis = cbind(1, L1.mis, A1, L2.mis)
      W1.mis = cbind(1, L1.mis)
      
      P = expit(cbind(1,L1,U) %*% rbind(const, 1, 0.5))
      Y = rbinom(n,1,eff(A3,A4,psi)*P)
    }
    
    res_i = matrix(NA, nrow = 1, ncol = 14)
    colnames(res_i) = c("psi4","psi3","psi2","psi1","EY0","EY1",
                        "pred0","pred0.2","pred1","pred1.2",
                        paste0("rss",4:1)) # 1:yes/no, 2:n of participants
    
    # prepare for the nuisance components
    pars4 = LSEst(Y, A4, Z4, W4)
    PS4.model = glm(A4 ~ W4, family = binomial(link = "logit"))
    PS4 = predict(PS4.model, type = "response")
    
    # solve the g-estimation equation
    psi4 = DREst(Y, A4, Z4, W4, pars4, PS4)
    
    # # re-update the g-estimator
    # psi4 = DREst_re(Y, A4, Z4, W4, pars4,psi4, PS4)
    
    res_i[1,1] = psi4$minimum
    res_i[1,11] = psi4$rss
    
    # calculate the pseudo-outcome
    H3 = Y * exp(-(Z4 %*% psi4$minimum) * A4)
    
    # obtain g-estimator at time 3 
    pars3 = LSEst(H3, A3, Z3, W3)
    PS3.model = glm(A3 ~ W3, family = binomial(link = "logit"))
    PS3 = predict(PS3.model, type = "response")
    psi3 = DREst(H3, A3, Z3, W3, pars3, PS3)
    res_i[1,2] = psi3$minimum
    res_i[1,12] = psi3$rss

    # repeat the same procedure at time 2 
    H2 = H3 * exp(-psi3$minimum * A3)
    pars2 = LSEst(H2, A2, Z2, W2)
    PS2.model = glm(A2 ~ W2, family = binomial(link = "logit"))
    PS2 = predict(PS2.model, type = "response")
    psi2 = DREst(H2, A2, Z2, W2, pars2, PS2)
    res_i[1,3] = psi2$minimum
    res_i[1,13] = psi2$rss
    
    # repeat the same procedure at time 1
    H1 = H2 * exp(-psi2$minimum * A2)
    pars1 = LSEst(H1, A1, Z1, W1)
    PS1.model = glm(A1 ~ W1, family = binomial(link = "logit"))
    PS1 = predict(PS1.model, type = "response")
    psi1 = DREst(H1, A1, Z1, W1, pars1, PS1)
    res_i[1,4] = psi1$minimum
    res_i[1,14] = psi1$rss
    
    EY0 = getProbRR(Z1 %*% psi1$minimum, W1 %*% pars1$alpha)[,1]
    EY1 = EY0 * exp(psi1$minimum + psi2$minimum + psi3$minimum + psi4$minimum)
    
    res_i[1,5] = mean(EY0)
    res_i[1,6] = mean(EY1)
    
    if(max(EY0) > 1) {
      res_i[1,7] = 1
      res_i[1,8] = sum(EY0 > 1)
    }else {res_i[1,7:8] = 0}
    
    if(max(EY1) > 1) {
      res_i[1,9] = 1
      res_i[1,10] = sum(EY1 > 1)
    }else {res_i[1,9:10] = 0}
    
    return(res_i)
  }

# ============================= #
# Dukes and Vansteelandt method #
# ============================= #

.Random.seed = seed_summary$SIM2_neg1

result = 
  foreach(i = 1:10, .combine = rbind,.packages = c("geeM", "tidyverse","brm","gesttools")) %dorng% {
    
    tryCatch({
      {
        n = 1000
        id = 1:n
        U = rnorm(n, 0, 1); L1 = rnorm(n, U, 1)
        A1 = rbinom(n, 1, 1/(1 + exp(0.2 - 0.2*L1)))
        L2 = rnorm(n, 0.5*L1 + 0.5*A1 + U,1)
        A2 = rbinom(n, 1, 1/(1 + exp(0.2 - 0.1*L1 - 0.2*L2 - 0.5*A1)))
        L3 = rnorm(n, 0.5*L2 + 0.5*A2 + U,1)
        A3 = rbinom(n, 1, 1/(1 + exp(0.2 - 0.1*L1 - 0.2*L3 - 0.5*A2)))
        L4 = rnorm(n, 0.5*L3 + 0.5*A3 + U,1)
        A4 = rbinom(n, 1, 1/(1 + exp(0.2 - 0.1*L1 - 0.2*L4 - 0.5*A3)))
        
        Z4 = Z3 = Z2 = Z1 = cbind(rep(1, n))
        W4 = cbind(1, L1, A3, L4)
        W3 = cbind(1, L1, A2, L3)
        W2 = cbind(1, L1, A1, L2)
        W1 = cbind(1, L1)
        
        # L1.mis = rnorm(n, U, 1)
        # L2.mis = rnorm(n, 0.5*L1.mis + 0.5*A1 + U,1)
        # L3.mis = rnorm(n, 0.5*L2.mis + 0.5*A2 + U,1)
        # L4.mis = rnorm(n, 0.5*L3.mis + 0.5*A3 + U,1)
        # 
        # W4.mis = cbind(1, L1.mis, A3, L4.mis)
        # W3.mis = cbind(1, L1.mis, A2, L3.mis)
        # W2.mis = cbind(1, L1.mis, A1, L2.mis)
        # W1.mis = cbind(1, L1.mis)
      }
      
      P = expit(cbind(1,L1,U) %*% rbind(const, 1, 0.5))
      Y = rbinom(n,1,eff(A3,A4,psi)*P)
      
      res_i = matrix(NA, nrow = 1, ncol = 14)
      colnames(res_i) = c("psi4","psi3","psi2","psi1","EY0","EY1",
                          "pred0","pred0.2","pred1","pred1.2", # 1:yes/no, 2:n of participants
                          paste0("rss",4:1)) 
      
      PS4.mod = glm(A4~W4,family=binomial)
      PS4 = predict(PS4.mod,type="response")
      gest4 = geem(Y ~ PS4 + A4 + L4 + A3 + L1,family=Gamma(link="log"),id = id)
      res_i[1,1] = gest4$beta[3]
      
      H3 = Y*exp(-gest4$beta[3]*A4)
      res_i[1,11] = sum((H3 - exp(cbind(1,as.vector(PS4),0,L4,A3,L1) %*% summary(gest4)$beta))^2)
      
      PS3.mod=glm(A3~W3,family=binomial)
      PS3=predict(PS3.mod,type="response")
      gest3 = geem(H3 ~ PS3 + A3 + L3 + A2 + L1,family=Gamma(link="log"),id = id)
      res_i[1,2] = gest3$beta[3]
      
      H2 = H3*exp(-gest3$beta[3]*A3)
      res_i[1,12] = sum((H2 - exp(cbind(1,as.vector(PS3),0,L3,A2,L1) %*% summary(gest3)$beta))^2)
      
      PS2.mod=glm(A2~W2,family=binomial)
      PS2=predict(PS2.mod,type="response")
      gest2 = geem(H2 ~ PS2 + A2 + L2 + A1 + L1 ,family=Gamma(link="log"),id=id)
      res_i[1,3] = gest2$beta[3]
      
      H1 = H2*exp(-gest2$beta[3]*A2)
      res_i[1,13] = sum((H1 - exp(cbind(1,as.vector(PS2),0,L2,A1,L1) %*% summary(gest2)$beta))^2)
      
      PS1.mod=glm(A1~L1,family=binomial)
      PS1=predict(PS1.mod,type="response")
      gest1 = geem(H1 ~ PS1 + A1 + L1,family=Gamma(link="log"),id=id)
      res_i[1,4] = gest1$beta[3]

      H0 = H1*exp(-gest1$beta[3]*A1)
      EY0 = exp(cbind(1,as.vector(PS1),0,L1) %*% summary(gest1)$beta)
      EY1 = EY0 * exp(gest1$beta[3] + gest2$beta[3] + gest3$beta[3] + gest4$beta[3])
      
      res_i[1,14] = sum((H0 - EY0)^2)
      res_i[1,5] = mean(EY0)
      res_i[1,6] = mean(EY1)
      
      if(max(EY0) > 1) {
        res_i[1,7] = 1
        res_i[1,8] = sum(EY0 > 1)
      }else {res_i[1,7:8] = 0}
      
      if(max(EY1) > 1) {
        res_i[1,9] = 1
        res_i[1,10] = sum(EY1 > 1)
      }else {res_i[1,9:10] = 0}
      
      return(res_i)
    }, error = function(e) {
      message("Error: ", e$message)
      return(res_i)
    }
    )
  }

# Stop the parallel cluster
stopCluster(cl)
