scenariosimci_prop <- function(k,
                               cfj,
                               par,
                               vcovpara,
                               sc,
                               scenario,
                               tpred=0.25, # time of additional follow up (in years) after study ends. Default value: 0.25.    ,
                               tadd=0.25, # time of additional follow up (in years) after the first additional follow up time  in study. Default value: 0.25,
                               entry,
                               cencut=8
){
  #   ## resample time
  #   rd=sample(1:sized,size=sized,replace=TRUE,prob=NULL)
  #   dat=maindata[rd,c("time","censor")]
  
  ### simulate the parameters
  medimp2=scenario[k,1]
  crimp2=scenario[k,2]
  
  #par1c = c(par[k,1],par[k,2],par[k,3])
  #par2t = c(par[k,1],par[k,4],par[k,5])
  
  #   flag <- 0
  #   Par <- NULL
  #   while(flag < 1){
  #     Par <- rmvnorm(1, mean=par[k,],sigma=vcovpara[k][[1]])
  #     if(Par[1] > 0 & Par[2] > 0 & 0 < Par[3] & Par[3] < 1){
  #       flag <- flag + 1
  #       par1c <- Par
  #     }
  #   }
  
  Par <- rmvnorm(1,mean=c(log(par[k,][1]),log(par[k,][2]),log(par[k,][3]/(1-par[k,][3]))),
                 sigma=vcovpara[k][[1]])
  par1c <- c(exp(Par[1]),exp(Par[2]),((exp(Par[3]))/(1+exp(Par[3]))))
  par2t <- c(par1c[1], par1c[2]+medimp2*((log(2))^(-(1/par1c[1]))), par1c[3]+crimp2)
  
  ## simulate the event times
  
  cure1s <- rbinom((2*npat), size=1, par1c[3])
  x1s <- rweibull((2*npat), shape=par1c[1], scale=par1c[2])
  x1s[cure1s == 1] <- cencut # cured patients will be censored at yr 8 (cencut)
  
  #   cure2s <- rbinom(npat, size=1, par2t[3])
  #   x2s <- rweibull(npat,shape=par2t[1], scale=par2t[2])
  #   x2s[cure2s==1] <-cencut # cured patients will be censored at yr 8
  #  xts <- c(x1s,x2s)
  
  censs <- rep(0, length(x1s))
  censs[x1s + entry > cf[cfj]+tail(entry,1)] <- 1   
  x1s[censs ==1 ] <- cf[cfj] + tail(entry,1)-entry[censs ==1] 
  
  ### prob of getting control treatment : Pr(J=0)
  num=((((1-par1c[3])*dweibull(x1s,shape=par1c[1],scale=par1c[2]))^(1-censs)) * ((par1c[3] + ((1-par1c[3])*(1-pweibull(x1s,shape=par1c[1],scale=par1c[2]))))^censs))
  den=((((1-par2t[3])*dweibull(x1s,shape=par2t[1],scale=par2t[2]))^(1-censs)) * ((par2t[3] + ((1-par2t[3])*(1-pweibull(x1s,shape=par2t[1],scale=par2t[2]))))^censs))
  
  prJ0=num/(num+den)
  
  ind=unlist(lapply(1:length(prJ0),function(u) {rbinom(1,size=1,prob=prJ0[u])}))
  
  xts=x1s[ind==0]
  censs00=censs[ind==0]
  censs0=1-censs00
  
  dat=data.frame(time=xts,censor=censs0) 
  
  ########################################################
  dat$ctrsurvt = survcalc_prop(par1c,dat$time) # survival for control at a time point
  dat$trtsurvt = survcalc_prop(par2t,dat$time) # survival for treatment at a time point
  
  ## ctrprob = S1/(S1+S2)  ## 1: control group, 2: treatment group both at current time point
  # prob of receiving control for patients at risk =P(patient i gets control|at risk)
  dat$ctrprob <- dat$ctrsurvt/(dat$ctrsurvt+dat$trtsurvt)
  
  #################################################################
  ## total hazard of the risk set at a future time point
  ctrhaztps = slopehaz_prop(par1c,dat$time[dat$censor==1],dat,fu=tpred)
  
  #       ## total hazard of the risk set at an additional future time point
  #       ctrhaztpa = slopehaz(par1f,maindata$time,maindata,fu=tpred+tadd)
  
  ################################################################
  ## no of events between current time pt and a future time point
  ctrevetps = eventf_prop(par1c,dat$time[dat$censor==1],dat,fu=tpred)
  
  ## no of events between current time pt and an additional future time point
  ctrevetpas = eventf_prop(par1c,dat$time[dat$censor==1],dat,fu=tpred+tadd)
  
  ## maximum no of events between current time pt and an infinite future time point
  ctrevetmaxs = sum(dat$ctrprob[dat$censor==1]*(1-(par1c[3]/dat$ctrsurvt[dat$censor==1])))
  
  datos=data.frame(scenario=k,ctrhaztps,ctrevetps,ctrevetpas,ctrevetmaxs)
  #############################################################################
  
  datos
}

