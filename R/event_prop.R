eventf_prop <- function(cpar,
                        ctime,
                        maindat1,
                        fu=tpred # time of additional follow up (in years) after study ends. Default value: 0.25.    
){
  maindat=maindat1[maindat1$censor==1,]
  sum(maindat$ctrprob*((maindat$ctrsurvt - survcalc_prop(cpar,ctime+fu))/maindat$ctrsurvt))
}

