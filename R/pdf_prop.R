pdfcalc_prop <- function(par1,
                         ctime
){
  (1-par1[3])*(dweibull(ctime, shape=par1[1], scale=par1[2]))
}

