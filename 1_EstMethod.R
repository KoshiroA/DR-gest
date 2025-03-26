#########################
### Estimation method ###
#########################

# ==============================================================================
LSEst = function(Y,A,Z,W,max.step = 1000,thres = 1e-6){
  pa = dim(Z)[2]
  pb = dim(W)[2]
  psi.start = rep(0, pa)
  alpha.start = rep(0, pb)
  
  rss.psi = function(psi){
    p0p1 = getProbRR(Z %*% psi, W %*% alpha) # referred to brm package
    p0    = p0p1[,1];  p1 = p0p1[,2]
    return(sum((Y[A==0]-p0[A==0])^2)+sum((Y[A==1]-p1[A==1])^2))
  }
  
  rss.alpha = function(alpha){
    p0p1 = getProbRR(Z %*% psi, W %*% alpha) # referred to brm package
    p0    = p0p1[,1];  p1 = p0p1[,2]
    return(sum((Y[A==0]-p0[A==0])^2)+sum((Y[A==1]-p1[A==1])^2))
  }
  
  rss = function(psi,alpha){
    p0p1 = getProbRR(Z %*% psi, W %*% alpha) # referred to brm package
    p0    = p0p1[,1];  p1 = p0p1[,2]
    return(sum((Y[A==0]-p0[A==0])^2)+sum((Y[A==1]-p1[A==1])^2))
  }
  
  ## Optimization 
  Diff = function(x,y) {sum((x-y)^2)/sum(x^2+thres)}
  psi = psi.start; alpha = alpha.start
  diff = thres + 1; step = 0
  while(diff > thres & step < max.step){
    step = step + 1
    
    if(pa==1){
      opt1 = optimize(rss.psi,interval = c(-10,10),
                      tol = thres)
      diff1 = Diff(opt1$minimum,psi)
      psi = opt1$minimum
    }else{
      opt1 = stats::optim(psi,rss.psi,
                          control=list(maxit=max(100,max.step/10)))  
      diff1 = Diff(opt1$par,psi)
      psi = opt1$par
    }
    
    if(pb==1){
      opt2 = optimize(rss.alpha,interval = c(-3,3))
      diff  = max(diff1,Diff(opt2$minimum,alpha))
      alpha = opt2$minimum
    }else{
      opt2 = stats::optim(alpha,rss.alpha,
                          control=list(maxit=max(100,max.step/10)))  
      diff  = max(diff1,Diff(opt2$par,alpha))
      alpha = opt2$par
    }
    
  }
  pars2 = list(psi = psi,alpha=alpha, value = rss(psi,alpha),
               convergence = (step < max.step), step = step)
  return(pars2)
}
# ==============================================================================

# ==============================================================================
DREst = function(Y,A,Z,W,prepars,PS = NULL,
                 max.step = 1000,thres = 1e-6,startpars = NULL){
  
  if(is.null(PS)){
    PS.model = glm(A~W,family = binomial(link = "logit")) 
    PS = predict(PS.model,type="response")
  }
  
  if(is.null(startpars) == TRUE){
    psi = rep(0,dim(Z)[2])
  }else{
    psi = prepars$psi
  }
  
  Diff = function(a,b) {sum((a-b)^2)/sum(a^2+thres)}
  diff = thres + 1; step = 0
  while(diff > thres*1000 & step < max.step){
    step = step + 1
    
    if(dim(Z)[2] == 1){
      U.hat = getProbRR(Z %*% prepars$psi,W %*% prepars$alpha)[,1]
      gest = function(psi.dr){
        tmp = t(Z) %*% as.vector((A-PS) * (Y*exp(-A*psi.dr) - U.hat))
        return(sum(tmp^2))
      }
      opt = optimize(gest,c(-10,10),tol = thres)
      # opt$value = gest(opt$objective)
      diff = Diff(opt$minimum,psi)
      psi = opt$minimum
      opt$rss = sum((Y*exp(-A*(Z %*% opt$minimum)) - U.hat)^2)
      
    }else{
      U.hat = getProbRR(Z %*% prepars$psi,W %*% prepars$alpha)[,1]
      gest = function(psi.dr){
        tmp = t(Z) %*% as.vector((A-PS) * (Y*exp(-A*(Z %*% psi.dr)) - U.hat))
        return(sum(tmp^2))
      }
      opt = stats::optim(psi, gest, control=list(reltol=thres))
      diff = Diff(opt$par,psi)
      psi = opt$par
      opt$rss = sum((Y*exp(-A*(Z %*% opt$par)) - U.hat)^2)
    }
  }
  
  if(step < max.step){
    return(opt)
  }else{
    stop("The estimation equation is not properly be solved")
  }
}
# ==============================================================================
# ==============================================================================
DREst_re = function(Y,A,Z,W,prepars,psidr,PS = NULL,
                    max.step = 1000,thres = 1e-3,startpars = NULL){
  # re-update g-estimator
  Diff = function(a,b) {sum((a-b)^2)/sum(a^2+thres)}
  diff = thres + 1;  step = 0
  while (diff > thres & step < max.step) {
    step = step + 1
    prepars$psi = psidr$minimum
    psidr = DREst(Y, A, Z, W, prepars, PS)
    diff = Diff(psidr$minimum,prepars$psi)
  }
  
  if(step < max.step){
    psidr$step = step
    return(psidr)
  }else{
    stop("The estimation equation is not properly be solved")
  }
  
  
}
# ==============================================================================