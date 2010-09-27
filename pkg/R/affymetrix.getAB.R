affymetrix.getAB <-
function(logR, bfreq, R.BASELINE=2){
    R = exp(logR)*R.BASELINE
    B=R*bfreq
    A=R-B
    list(A=A,B=B)
}

