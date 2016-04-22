wblik13z_prop <- function(x, # vector of 2 weibull parameters and 1 cure rate := (k(shape), lambda(scale), r)
                          time, # vector
                          censor,# 1-censoring indicator
                          medimp=0.249, #initial median improvement(in years),
                          crimp=0.05 # initial cure rate improvement(in % eg: 5% = 0.05)
) {
  k1=exp(x[1])
  l1=exp(x[2])
  r1=(exp(x[3]))/(1+exp(x[3]))
  
  l2=l1+medimp*((log(2))^-(1/k1))
  r2=r1+crimp
  n=length(time)
  
  lcen <- sum((1-censor)*log(((1-r1)*dweibull(time,shape=k1,scale=l1))+((1-r2)*dweibull(time,shape=k1,scale=l2))))
  
  luncen <-sum(censor*log(r1+r2+((1-r1)*(1-pweibull(time,shape=k1,scale=l1)))+((1-r2)*(1-pweibull(time,shape=k1,scale=l2)))))
  
  nlnL <- (n*log(2)) - lcen - luncen
  nlnL
  
}

