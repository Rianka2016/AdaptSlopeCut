#' A function to calculate the optimal clinical cut-off
#' 
#' This command helps in calculating the cut-off point in a simulation study for cure rate survival data.
#' @param npat The total number of patients in the study.
#' @param par1 The true parametric value for the control group: (k,l,r) where k and l are the shape and scale parameters for the weibull distribution and r is the cure rate for the control treatment. Default values: (1.33,0.77,0.40).  
#' @param medianimp True improvement in median survival rate for uncured patients while using the experimental treatment. Default value: 0.249.    
#' @param cureimp True improvement in cure rate while using the experimental treatment. Default value: 0.05.    
#' @param cencut Time point (in years from the start of the study) after which a patient will no longer be followed. 
#' @param tpred Time of additional follow up (in years) after study ends. Default value: 0.25.    
#' @param tadd Time of additional follow up (in years) after the first additional follow up time  in study. Default value: 0.25.     
#' @param cf A vector representing all the time points (in years) after the last patient in, when the study can be stopped.   
#' @param scenario A dataframe representing different effect-size scenarios for diff combination of median and cure rate improvement. First column called "impmed" has different values for median improvement and the second column "impcure" represents the different values of cure-rate improvement.   
#' @param nsamp Number of resamples needed to calculate confidence interval. Default: 100
#' @param fCI A dummy to calculate confidence intervals. 1=yes, 0=no. Default: 0
#' @param lpi Last patient in (in years calculated from the start of the study). Default: 2.5.    
#' @param entry Entry time of parients, twice for control group and then treatment group. Default: rep(lpi/npat * 1:npat, 2)    
#' @param al (1-al)*100 \% CI is calculated. Default value: 0.05.   
#' 
#' @return It returns a dataframe.
#' @return cutoff: The position of the observation in the vector "cf" which is also the estimated cut-off point.   
#' @return param_k: The estimated parametric value of the shape parameter for the weibull distribution for both treatment groups.  
#' @return ctrparam_l: The estimated parametric value of the scale parameter for the weibull distribution for control group. 
#' @return ctrparam_r: The estimated parametric value of the cure rate of the control treatment. 
#' @return trtparam_l: The estimated parametric value of the scale parameter for the weibull distribution for treatment group. 
#' @return trtparam_r: The estimated parametric value of the cure rate of the experimental treatment.
#' @return ctrhaztp: Total hazard of the risk set of the control group at a future time point.   
#' @return ctrevetp: The estimated number of events to occur in the control group at additional follow-up time point "tpred" after the end of the study.  
#' @return ctrevetpa: The estimated number of events to occur in the control group at additional follow-up time point "tpred + tadd" after the end of the study. 
#' @return ctrevetmax: The estimated maximum number of events to occur in the control group with additional follow-up after the end of the study. 
#' @return mean_ctrhaz: Mean total hazard of the risk set of the control group at a future time point for the resampled data.  
#' @return med_ctrhaz: Median total hazard of the risk set of the control group at a future time point for the resampled data.
#' @return mean_ctrevep: Mean estimated number of events to occur in the control group at additional follow-up time point "tpred" after the end of the study for the resampled data.   
#' @return med_ctrevep: Median estimated number of events to occur in the control group at additional follow-up time point "tpred" after the end of the study for the resampled data.   
#' @return mean_ctrevepa: Mean estimated number of events to occur in the control group at additional follow-up time point "tpred + tadd" after the end of the study for the resampled data.    
#' @return med_ctrevepa: Median estimated number of events to occur in the control group at additional follow-up time point "tpred + tadd" after the end of the study for the resampled data.    
#' @return mean_ctrevemax: Mean estimated maximum number of events to occur in the control group with additional follow-up after the end of the study for the resampled data.   
#' @return med_ctrevemax: Median estimated maximum number of events to occur in the control group with additional follow-up after the end of the study for the resampled data.  
#' @return ctrhaz_l: Lower limit of the total hazard of the risk set of the control group at a future time point for the resampled data.
#' @return ctrhaz_u: Upper limit of the total hazard of the risk set of the control group at a future time point for the resampled data.   
#' @return ctrevep_l: Lower limit of the estimated number of events to occur in the control group at additional follow-up time point "tpred" after the end of the study for the resampled data.   
#' @return ctrevep_u: Upper limit of the estimated number of events to occur in the control group at additional follow-up time point "tpred" after the end of the study for the resampled data.   
#' @return ctrevepa_l: Lower limit of the estimated number of events to occur in the control group at additional follow-up time point "tpred + tadd" after the end of the study for the resampled data.    
#' @return ctrevepa_u: Upper limit of the estimated number of events to occur in the control group at additional follow-up time point "tpred + tadd" after the end of the study for the resampled data.    
#' @return ctrevemax_l: Lower limit of the estimated maximum number of events to occur in the control group with additional follow-up after the end of the study for the resampled data.  
#' @return ctrevemax_u: Upper limit of the estimated maximum number of events to occur in the control group with additional follow-up after the end of the study for the resampled data.  
#' 
#' @examples scenarios <- data.frame(impmed = rep(c(0.2,0.3, 0.4), rep(3, 3)), 
#'                       impcure = rep(c(0.1, 0.05, 0), 3))  
#' cuts <- c(1,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.7,1.8)  
#' optcutoff_prop(npat=100,scenario=scenarios,cf=cuts) # gives only estimates 
#' optcutoff_prop(npat=100,scenario=scenarios,cf=cuts,fCI=1) # gives estimates along with confidence intervals
#' 
#' @export
optcutoff_prop <- function(npat, # total number of patients.  
                           cf, # vector representing all the time points (in years) after the last patient in, when the study can be stopped.   
                           scenario, # dataframe representing different effect-size scenarios for diff combination of median and cure rate improvement.First column called "impmed" has different values for median improvement and the second column "impcure" represents the different values of cure-rate improvement.   
                           par1=c(1.33,0.77,0.40), # true parameters for the control group: (k,l,r) where k and l are the shape and scale parameters for the weibull distribution and r is the cure rate for the control treatment. Default values: (1.33,0.77,0.40).  
                           medianimp=0.249, # true improvement in median survival rate for uncured patients while using the experimental treatment. Default value: 0.249.    
                           cureimp=0.05, # true improvement in cure rate while using the experimental treatment. Default value: 0.05.    
                           cencut=8, # time point (in years from the start of the study) after which a patient will no longer be followed.    
                           tpred=0.25, # time of additional follow up (in years) after study ends. Default value: 0.25.    
                           tadd=0.25, # time of additional follow up (in years) after the first additional follow up time  in study. Default value: 0.25.     
                           nsamp=100, # no of resamples needed to calculate confidence interval. Default: 100
                           fCI=0, # dummy to calculate CI or not 1=yes, 0=no. Default: 0
                           lpi=2.5, # last patient in (in years calculated from the start of the study). Default: 2.5.    
                           entry = rep(lpi/npat * 1:npat, 2), # entry time of parients, twice for control group and then treatment group. Default: rep(lpi/npat * 1:npat, 2)    
                           al=0.05 # (1-al)*100 % CI is calculated. Default value: 0.05.      
){
  nscen <- dim(scenario)[1] # no of scenarios
  par2 = c(par1[1], par1[2]+(medianimp*((log(2))^(-(1/par1[1])))), par1[3]+cureimp )
  
  # generate the cure indicator for control
  # 1: cured, 0: not cured
  cure1 <- rbinom(npat, size=1, par1[3])
  # generate weibull survival data for control
  x1 <- rweibull(npat, shape=par1[1], scale=par1[2])
  # generate mixture distribution data for control
  x1[cure1 == 1] <- cencut # cured patients will be censored at yr 8 (cencut)
  
  # generate the cure indicator for treatment
  cure2 <- rbinom(npat, 1, par2[3])
  # generate weibull survival data for control
  x2 <- rweibull(npat,shape=par2[1], scale=par2[2])
  # generate mixture distribution data for control
  x2[cure2==1] <-cencut # cured patients will be censored at yr 8
  
  xt0 <- c(x1,x2)
  cure0 <- c(cure1, cure2)
  
  for(j in 1:length(cf)){
    xt = xt0
    censs <- rep(0, length(xt))
    censs[xt + entry > cf[j]+tail(entry,1)] <- 1   
    xt[censs ==1 ] <- cf[j] + tail(entry,1)-entry[censs ==1] 
    
    cind <- 1:length(x1)
    
    cens=1-censs
    
    maindata=data.frame(time=xt,censor=cens)
    
    cevtact <- data.frame(CUTOFF = cf[j], cact = sum(cens[cind] == 1 & cure1[cind] == 0 & (x1 + entry[cind] <= cf[j]+tpred + tail(entry[cind],1))),
                          cactadd = sum(cens[cind] == 1 & cure1[cind] == 0 & (x1 + entry[cind] <= cf[j]+tpred + tadd +tail(entry[cind],1))),
                          cactmax = sum(cens[cind] == 1 & cure1[cind] == 0)) 
    print(cevtact)
    istop <- rep(1, nscen) 
    
    scenariosum <- llply(1:nscen, function(k){ scenariosim_prop(j,k,par1,maindata,scenario,tpred,tadd,istop)})
    
    sc=matrix(0,nrow=nscen,ncol=11)
    for(i in 1:nscen){
      sc[i,]=unlist(scenariosum[i][[1]]$dato[,])
    }
    
    colnames(sc) <- c("cutoff","scenario","param_k","ctrparam_l","ctrparam_r",
                      "trtparam_l","trtparam_r","ctrhaztp","ctrevetp","ctrevetpa",
                      "ctrevetmax")
    
    istop=as.numeric(sc[,"ctrevetp"] <=2 & sc[,"ctrevetpa"] <=3 & sc[,"ctrevetmax"] <=5)
    
    #scenariosum$istop=istop
    scenariosum
    sc
    
    if(mean(istop) >=0.6) break
  }
  #scenariosum
  
  vcovpara=list()
  
  for(i in 1:nscen){
    vcovpara[[i]]=scenariosum[i][[1]]$vcovpar
  }
  
  
  if(fCI == 1){## 95% confidence interval calculation
    cfj=unique(sc[,"cutoff"])
    
    par=sc[,c("param_k", "ctrparam_l", "ctrparam_r")]
    
    
    msd=ldply(1:nscen, function(k) {mulsimci_prop(nsamp,k,cfj,par,vcovpara,sc,scenario,tpred,tadd,al,entry,cencut=8)})
    
    return(cbind(sc,msd))
    
  }
  
  else{return(sc)}
}

