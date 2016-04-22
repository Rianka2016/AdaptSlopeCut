scenariosim_prop <- function(j,
                             k,
                             par1=c(1.33,0.77,0.40), # true parameters for the control group: (k,l,r) where k and l are the shape and scale parameters for the weibull distribution and r is the cure rate for the control treatment. Default values: (1.33,0.77,0.40).,
                             maindata,
                             scenario,
                             tpred=0.25, # time of additional follow up (in years) after study ends. Default value: 0.25.    
                             tadd=0.25, # time of additional follow up (in years) after the first additional follow up time  in study. Default value: 0.25,
                             istop
){
  medimp1 = scenario[k,1]
  crimp1 = scenario[k,2]
  
  #wblik23 <- function(x){ wblik13(x,time=maindata$time,censor=maindata$censor,medimp=medimp1,crimp=crimp1)}
  wblik23z <- function(x){ wblik13z_prop(x,time=maindata$time,censor=maindata$censor,medimp=medimp1,crimp=crimp1)}
  
  ini=c(log(par1[1]), log(par1[2]), log(par1[3]/(1-par1[3])))
  fit1 = nlminb(start=ini,objective=wblik23z)
  tfitpar=fit1$par
  
  hesspar=hessian(wblik23z,tfitpar)
  vcovpar=solve(hesspar)
  
  par1f=c(exp(tfitpar[1]), exp(tfitpar[2]), ((exp(tfitpar[3]))/(1+exp(tfitpar[3]))))
  par2f = c(par1f[1], par1f[2]+medimp1*((log(2))^(-(1/par1f[1]))), par1f[3]+crimp1)
  
  #   #   fit1 = optim(par=par1,fn=wblik2,method="L-BFGS-B",lower=c(0,0,0),upper=c(Inf,Inf,1),
  #   #                hessian=TRUE)
  #   
  #   #   fit1 = optim(par=par1,fn=wblik2,method="BFGS",hessian=TRUE)
  #   
  #   #   fit1 = optgrid(wblik2, par=par1,lower=c(0,0,0),upper=c(Inf,Inf,1),verbose=1)
  #   
  #   
  #   par1f = fit1$par ## optimum parameter estimate for control
  #   hess1 = fit1$hessian ## hessian matrix
  #   # optimum parameter estimate for trt group for a fixed improvement in median and cure rate
  #   par2f = c(par1f[1], par1f[2]+(medimp1*log(2)^(-(1/par1f[1]))), par1f[3]+crimp1 )
  
  maindata$ctrsurvt = survcalc_prop(par1f,maindata$time) # survival for control at a time point
  maindata$trtsurvt = survcalc_prop(par2f,maindata$time) # survival for treatment at a time point
  
  ## ctrprob = S1/(S1+S2)  ## 1: control group, 2: treatment group both at current time point
  # prob of receiving control for patients at risk =P(patient i gets control|at risk)
  maindata$ctrprob <- maindata$ctrsurvt/(maindata$ctrsurvt+maindata$trtsurvt)
  
  #################################################################
  ## total hazard of the risk set at a future time point
  ctrhaztp = slopehaz_prop(par1f,maindata$time[maindata$censor==1],maindata,fu=tpred)
  
  #       ## total hazard of the risk set at an additional future time point
  #       ctrhaztpa = slopehaz(par1f,maindata$time,maindata,fu=tpred+tadd)
  
  ################################################################
  ## no of events between current time pt and a future time point
  ctrevetp = eventf_prop(par1f,maindata$time[maindata$censor==1],maindata,fu=tpred)
  
  ## no of events between current time pt and an additional future time point
  ctrevetpa = eventf_prop(par1f,maindata$time[maindata$censor==1],maindata,fu=tpred+tadd)
  
  ## maximum no of events between current time pt and an infinite future time point
  ctrevetmax = sum(maindata$ctrprob[maindata$censor==1]*(1-(par1f[3]/maindata$ctrsurvt[maindata$censor==1])))
  
  dato=data.frame(cutoff=j,scenario=k,param_k=round(par1f[1],4),
                  ctrparam_l=round(par1f[2],4),ctrparam_r=round(par1f[3],4),
                  trtparam_l=round(par2f[2],4),trtparam_r=round(par2f[3],4),
                  ctrhaztp,ctrevetp,ctrevetpa,ctrevetmax)
  #############################################################################
  
  return(list(dato=dato,vcovpar=vcovpar))
  
}
