illumina.getAB <-
function(logR, bfreq, R.BASELINE=2){
    R = exp(logR)
    R = R*R.BASELINE
    A=R/(1+tan(bfreq*pi/2))
    B=R-A
    list(A=A,B=B)
}

