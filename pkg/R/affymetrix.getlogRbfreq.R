affymetrix.getlogRbfreq <-
function(A,B,OFFSET=2){
    logR = log((A+B)/2)
    bfreq = B/(A+B)
    list(logR=logR,bfreq=bfreq)
}

