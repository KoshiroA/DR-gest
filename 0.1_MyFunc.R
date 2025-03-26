### MyFunc ###

expit = function(x){1/(1 + exp(-x))}
logit = function(x){log(x/(1 -x))}

eff = function(A3,A4,psi){
  exp(cbind(A3,A4) %*% psi)
}

P.func = function(lower,upper,out = c("a","P","summary")){
  a0 = -logit(lower)
  a2 = -(logit(lower) - logit(upper))/5
  a1 = a2*2
  #P = 1/(1 + exp(a0 - a1*L1 - a2*U))
  P = expit(-a0 + a1*L1 + a2*U)
  
  if(out == "a"){return(c(-a0,a1,a2))}
  else if(out == "summary"){list(nul = summary(P),neg = summary(P*exp(-0.5)),
                                 pos = summary(P*exp(0.5)))}
  #{c(EY0 = mean(P),pos = mean(P)*exp(0.5),neg = mean(P)*exp(-0.5))}
  else if(out == "P"){return(P)}
  #return(cbind(a0,a1,a2,P))
}