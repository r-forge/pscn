LengthEnough <-
function(vec, MinLength){
  n = length(vec)
  if (n <=1){
    return(vec)
  }else{
    keepchpt = rep(1,n)
    temp = vec[2:n]-vec[1:(n-1)]
    for (i in 1:(n-1)){
      if (temp[i] <= MinLength){
        keepchpt[i+1] = 0
      }
    }
    return(vec[which(keepchpt==1)])
  }
}

