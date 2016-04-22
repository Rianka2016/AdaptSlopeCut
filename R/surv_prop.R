survcalc_prop <- function(par1,
                          ctime
){
  par1[3]+ (1-par1[3])*(1-pweibull(ctime, shape=par1[1], scale=par1[2]))
}

