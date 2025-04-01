
### Perform doubly robust g-estimation

# The INPUT of 'prepars' is obtained by performing LSEst() function in advance

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