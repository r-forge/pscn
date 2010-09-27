initialmixstate <-
function(cnvobj){
  n = length(cnvobj$rawdata$theta)
  mixstate = rep(4,n)
  id2 = which(cnvobj$rawdata$theta<0.99)
  mixstate[id2] = 3
  id3 = which(cnvobj$rawdata$theta<0.5)
  mixstate[id3] = 2
  id4 = which(cnvobj$rawdata$theta<0.01)
  mixstate[id4] = 1
  mixstate
}

