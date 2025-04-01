
### Update the resulting g-estimator iteratively

# The INPUT of 'prepars' is obtained by performing LSEst() function in advance
# The INPUT of 'psidr' is the initial value of the iteration, and is the output DREst()

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