mulsimci_prop <- function(nsamp,
                          k,
                          cfj,
                          par,
                          vcovpara,
                          sc,
                          scenario, # dataframe representing different effect-size scenarios for diff combination of median and cure rate improvement.First column called "impmed" has different values for median improvement and the second column "impcure" represents the different values of cure-rate improvement.    
                          tpred=0.25, # time of additional follow up (in years) after study ends. Default value: 0.25.    
                          tadd=0.25, # time of additional follow up (in years) after the first additional follow up time  in study. Default value: 0.25
                          al=0.05, # (1-al)*100 % CI is calculated. Default value: 0.05.
                          entry = rep(lpi/npat * 1:npat, 2), # entry time of parients, twice for control group and then treatment group.
                          cencut=8 # time point (in years) after which a patient will no longer be followed.
){
  output=replicate(nsamp, scenariosimci_prop(k,cfj,par,vcovpara,sc,scenario,tpred,tadd,entry,cencut=8),simplify=TRUE)
  output=matrix(unlist(output),nrow=5,ncol=nsamp,byrow=F)
  
  data.out=data.frame(mean_ctrhaz=mean(output[2,]),med_ctrhaz=median(output[2,]),
                      ctrhaz_l=quantile(output[2,],probs=al/2),ctrhaz_u=quantile(output[2,],probs=1-(al/2)),
                      mean_ctrevep=mean(output[3,]),med_ctrevep=median(output[3,]),
                      ctrevep_l=quantile(output[3,],probs=al/2),ctrevep_u=quantile(output[3,],probs=1-(al/2)),
                      mean_ctrevepa=mean(output[4,]),med_ctrevepa=median(output[4,]),
                      ctrevepa_l=quantile(output[4,],probs=al/2),ctrevepa_u=quantile(output[4,],probs=1-(al/2)),
                      mean_ctrevemax=mean(output[5,]),med_ctrevemax=median(output[5,]),
                      ctrevemax_l=quantile(output[5,],probs=al/2),ctrevemax_u=quantile(output[5,],probs=1-(al/2)))
  
  data.out

}

