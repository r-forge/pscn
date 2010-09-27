pscn.plot <-
function(cnvobj=NULL, which.plot="bfreq", region=FALSE, regionid=NULL, region.col=NULL, scatter=FALSE, Contour=FALSE, use.main=FALSE, cndata=NULL, chrid=NULL, use.pos=FALSE, loc=NULL, changepoint=TRUE, col.gain="red", col.loss="blue", col.normal="green", ...){
  if (!is.null(loc) && region && !is.null(regionid)){
    cat("Please only specify loc or regionid.")
  }
  if (is.null(cndata)){
    cat("Please specify copy number variation list (eg. shortlist or longlist.update).\n")
  }else{
    if (use.main){
      if (!is.null(cnvobj)){
        mainlab = cnvobj$label
      }else{
        mainlab = paste(cndata$Sample[1]," Chr",chrid,sep="")
      }
    }else{
      mainlab = NULL
    }
    if (which.plot=="copy number"){
      plot.cnv3(mydata=cndata, chrid=chrid, pos=use.pos, gaincol=col.gain, normalcol=col.normal, losscol=col.loss, mainlab=mainlab, locs=loc,...)
    }else{
      if(is.null(loc)){
        loc=c(1:length(cnvobj))
      }
      if (region || scatter || Contour){
  #      chpts = c(1, cnvobj$chpts,(length(cnvobj)+1))
        tempid = which(cndata$Chr==strsplit(cnvobj$label,"Chr")[[1]][2])
        cndata = cndata[tempid,] 
        chpts = c(cndata$SNP.begin)
        if (length(tempid)>0){
          tempn = length(tempid)
          for (kk in 1:tempn){
            if (is.na(match((cndata$SNP.end[kk]+1),chpts))){
              chpts = c(chpts,(cndata$SNP.end[kk]+1))
            }
          }
        }
        chpts = sort(chpts)
        if (is.na(match(1,chpts))){
          chpts = c(1,chpts)
        }
        id = which(!is.na(match(chpts,loc)))
        n = length(id)
        if (n>0){
          chpts = chpts[id]
          if (chpts[n]!=loc[length(loc)]){
            chpts = c(chpts,(loc[length(loc)]+1))
            n = n+1
          }
          if (chpts[1]!=loc[1]){
            chpts = c(loc[1], chpts)
            n = n+1
          }
          n = n-1
        }else{
          chpts = c(loc[1],(loc[length(loc)]+1))
          n = 1
        }
        if (is.null(regionid)){
          regionid = 1:n
        }else{
          chpts = chpts[c(regionid,(regionid[length(regionid)]+1))]
          loc = chpts[1]:(chpts[length(chpts)]-1)
        }
        reg = rep(0, ((chpts[length(chpts)])-chpts[1]))
        n = length(regionid)
        for (i in 1:n){
          reg[(chpts[i]-chpts[1]+1):(chpts[i+1]-chpts[1])]=i
        }
        if (is.null(region.col)){
          region.col = rainbow(n)
        }
        reg.col = region.col[reg]
        cols = reg.col
      }else{
        cols = cnvobj$newstate[loc]
        if (which.plot=="logR" || which.plot=="R"){
          cols = rep("black", length(cols))
        }
      }
      
      if (Contour){
        locs = which(reg==1)
        f1 = kde2d(cnvobj$rawdata$A[locs], cnvobj$rawdata$B[locs], n=100, h=c(0.7,0.7))
        contour(f1, levels=c(0.2), drawlabels=FALSE, col=region.col[1], lwd=3, xlab="A", ylab="B", cex.lab=1.5, xlim=c(0,5), ylim=c(0,5), main=mainlab)
        if (n>1){
          for(ind in 2:n){
            locs=which(reg==ind)
            f1 = kde2d(cnvobj$rawdata$A[locs],cnvobj$rawdata$B[locs], n=100, h=c(0.7,0.7))
            contour(f1,levels=c(0.2), drawlabels=FALSE, col=region.col[ind], add=TRUE, lwd=3)
          }
        }    
      }else if (scatter){
        if (region){
          plot(cnvobj$rawdata$A[loc], cnvobj$rawdata$B[loc], col=reg.col, ylim=c(0,max(5,max(cnvobj$rawdata$B[loc]))),xlim=c(0,max(5,max(cnvobj$rawdata$A[loc]))),ylab="B", xlab="A", cex.lab=1.5, main=mainlab)
        }else{
          plot(cnvobj$rawdata$A[loc], cnvobj$rawdata$B[loc], col=cnvobj$newstate[loc], ylim=c(0,max(5,max(cnvobj$rawdata$B[loc]))),xlim=c(0,max(5,max(cnvobj$rawdata$A[loc]))),ylab="B", xlab="A", cex.lab=1.5, main=mainlab)
        } 
      }else{    
        if (use.pos){
          if (which.plot=="A"){
            plot(cnvobj$pos[loc]/1000000, cnvobj$rawdata$A[loc], col=cols,xlab="Position (Mb)", ylab="A intensity", main=mainlab,...)
          }
          if (which.plot=="B"){
            plot(cnvobj$pos[loc]/1000000, cnvobj$rawdata$B[loc], col=cols,xlab="Position (Mb)", ylab="B intensity", main=mainlab,...)
          }
          if (which.plot=="logR"){
            plot(cnvobj$pos[loc]/1000000, cnvobj$rawdata$logR[loc], col=cols,xlab="Position (Mb)", ylab="logR", main=mainlab,...)
          }
          if (which.plot=="R"){
            plot(cnvobj$pos[loc]/1000000, exp(cnvobj$rawdata$logR[loc])*2, col=cols,xlab="Position (Mb)", ylab="R", main=mainlab,...)
          }
          if (which.plot=="bfreq"){
            plot(cnvobj$pos[loc]/1000000, cnvobj$rawdata$theta[loc], col=cols,xlab="Position (Mb)", ylab="BAF", main=mainlab,...)
          }
          if (changepoint && !region){
            chptall = sort(c(cnvobj$chpts,cnvobj$SNP.star))
            n = length(chptall)
            for (i in 1:n){
              abline(v=cnvobj$pos[chptall[i]])
            }
          }
        }else{
          if (which.plot=="A"){
            plot(loc, cnvobj$rawdata$A[loc], col=cols,xlab="SNP Index", ylab="A intensity", main=mainlab,...)
          }
          if (which.plot=="B"){
            plot(loc, cnvobj$rawdata$B[loc], col=cols,xlab="SNP Index", ylab="B intensity", main=mainlab,...)
          }
          if (which.plot=="logR"){
            plot(loc, cnvobj$rawdata$logR[loc], col=cols, xlab="SNP Index", ylab="logR", main=mainlab,...)
          }
          if (which.plot=="R"){
            plot(loc, exp(cnvobj$rawdata$logR[loc])*2, col=cols, xlab="SNP Index", ylab="R", main=mainlab,...)
          }
          if (which.plot=="bfreq"){
            plot(loc, cnvobj$rawdata$theta[loc], col=cols,xlab="SNP Index", ylab="BAF", main=mainlab,...)
          }
          if (changepoint & !region){
            tempid = which(cndata$Chr==strsplit(cnvobj$label,"Chr")[[1]][2])
            cndata = cndata[tempid,] 
            chptall = c(cndata$SNP.begin,cnvobj$SNP.star)
            if (length(tempid)>0){
              tempn = length(tempid)
              for (kk in 1:tempn){
                if (is.na(match((cndata$SNP.end[kk]+1),chptall)) && cndata$SNP.end[kk]!=length(cnvobj)){
                  chptall = c(chptall,(cndata$SNP.end[kk]+1))
                }
              }
            }
            chptall = sort(chptall)
            n = length(chptall)
            for (i in 1:n){
              abline(v=chptall[i])
            }
          }
        }
      }
    }
  }
}
