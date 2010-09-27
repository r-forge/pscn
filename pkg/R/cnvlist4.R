cnvlist4 <-
function(samplename, chrs=1:23, talpha1=0.01, talpha2 = 0.01, cutvalue=1e-7, pi1=0.5, mu1.initial=1.1, mu2.initial=0.9,sigma1.initial=1, sigma2.initial=1, cutdev=0.4, MIN.SNPS=20){
# cutdev=0.4; level0 = 0.8; level1 = 1.2; cut1=1.9; cut2=2.1; talpha1=0.02; talpha2 = 0.01; cutvalue=1e-7; pi1=0.5; pi2=0.5; mu1.initial=1.1; mu2.initial=0.9; sigma1.initial=1; sigma2.initial=1 ; MIN.SNPS=20
# talpha: the alpha for t-test
# cutvalue: delta for iteration stop
  pi2 = 1-pi1
  if (!is.na(match(23,chrs))){
    load(paste(samplename,".Chr23.pscn0.Rdata",sep=""))
    if (cnvobj$gender == "Male"){
      chrs = chrs[is.na(match(chrs,23))]
      print("The gender is Male, so Chromosome 23 will not be analyzed.")
    }
  }
  # calculate Rmedian
  Rvec = c()
  for (chrid in chrs){
    load(paste(samplename,".Chr",chrid,".pscn0.Rdata",sep=""))
    Rvec = c(Rvec, exp(cnvobj$rawdata$logR)*2)
  }
  Rmedian = median(Rvec)
  total0 = Rmedian-cutdev
  total1 = Rmedian+cutdev
  # determine the normal state
  norm.total = c()
  for (chrid in chrs){
    load(paste(samplename,".Chr",chrid,".pscn0.Rdata",sep=""))
    cnvobj$Rmedian = Rmedian
    save(cnvobj,file=paste(cnvobj$label,".pscn.Rdata",sep=""))
    if (!is.null(cnvobj$clus)){
      n = length(cnvobj$clus)
      normid = rep(0,n) 
      
      N = length(cnvobj$newstate)
     
      chptall = c(1,cnvobj$chpts,(N+1))
      SNP.star = cnvobj$SNP.star
      if (!is.na(SNP.star)){
        if (length(which(chptall==(SNP.star+1)))==0){
          chptall = sort(c(chptall,(SNP.star+1)))
      # 3/7/2009 if add SNP.star, there may be possible that some segments have less that MIN.SNPs snps.
          throw = match(SNP.star+1,chptall)
          keepid = rep(1,length(chptall)) 
          if ((chptall[throw]-chptall[throw-1]) < MIN.SNPS){  
            keepid[throw-1] = 0
          }
          if ((chptall[throw+1]-chptall[throw]) < MIN.SNPS){
            keepid[throw+1] = 0
          }
          chptall = chptall[which(keepid==1)]            
        }  
      }              
      med = ABmedian(cnvobj,chptall)
           
      for (i in 1:n){
        if (cnvobj$clus[i]==1 && med$total.median[i]>total0 && med$total.median[i]<total1){
          normid[i] = 1
        }
      }
      if (sum(normid)>0){
        norm2 = c()
        norm3 = c()
        for (i in 1:n){
          if (normid[i]==1){
            temp2 = which(cnvobj$newstate[med$beginpt[i]:med$endpt[i]]==2)+med$beginpt[i]-1
            temp3 = which(cnvobj$newstate[med$beginpt[i]:med$endpt[i]]==3)+med$beginpt[i]-1
            norm2 = c(norm2,temp2)
            norm3 = c(norm3,temp3)
          }
        }
        norm.total = c(norm.total,cnvobj$rawdata$A[norm2],cnvobj$rawdata$A[norm3], cnvobj$rawdata$B[norm2],cnvobj$rawdata$B[norm3])
      }
    }
  }
  mu0 = mean(norm.total)
  sigma0 = sd(norm.total)
  norm.length = length(norm.total)
  
  longlist = c()
  shortlist = c()
  
  for (chrid in chrs){
    load(paste(samplename,".Chr",chrid,".pscn.Rdata",sep=""))
    N = length(cnvobj$intensity[,1])
     
    chptall = c(1,cnvobj$chpts,(N+1))
  
    SNP.star = cnvobj$SNP.star
    if (!is.na(SNP.star)){
      if (length(which(chptall==(SNP.star+1)))==0){
        chptall = sort(c(chptall,(SNP.star+1)))
        # 3/7/2009 if add SNP.star, there may be possible that some segments have less that MIN.SNPs snps.
        throw = match(SNP.star+1,chptall)
        keepid = rep(1,length(chptall)) 
        if ((chptall[throw]-chptall[throw-1]) < MIN.SNPS){  
          keepid[throw-1] = 0
        }
        if ((chptall[throw+1]-chptall[throw]) < MIN.SNPS){
          keepid[throw+1] = 0
        }
        chptall = chptall[which(keepid==1)]     
      }
    }              
    med = ABmedian(cnvobj,chptall)
    SNP.star = cnvobj$SNP.star
    Effective.pos = cnvobj$Effective.pos            
    if (length(med$beginpt)!=length(cnvobj$clus)){     
      print(paste(samplename, "Chr", chrid))
    }
    if (is.null(cnvobj$clus)){
      longtemp = c()
      shorttemp = c()
    }else{
      longtemp = c()
      n = length(chptall)-1
      count = 0
      for (k in 1:n){
        temp = cnvobj$newstate[chptall[k]:(chptall[k+1]-1)]
        id2 = which(temp==2)+chptall[k]-1
        id3 = which(temp==3)+chptall[k]-1
        n1 = length(id2)
        n2 = length(id3) 
        if (n1>0 && n2>0){
          count = count+1
          if (cnvobj$clus[count]==2 || (cnvobj$clus[count]==1 && (med$total.median[count]>total1 || med$total.median[count]<total0))){
            y = c(cnvobj$rawdata$A[id2],cnvobj$rawdata$A[id3],cnvobj$rawdata$B[id2],cnvobj$rawdata$B[id3])
            ny = length(y)
            mu1 = mu1.initial
            mu2 = mu2.initial
            sigma1 = sigma1.initial
            sigma2 = sigma2.initial
            pi1.vec.new = rep(0,ny)
            pi2.vec.new = rep(0,ny)
            delta = 1
            for (i in 1:ny){
                pi1.vec.new[i] = exp(-(mu1-y[i])^2/(2*sigma1^2))*pi1/sigma1/(exp(-(mu1-y[i])^2/(2*sigma1^2))*pi1/sigma1+exp(-(mu2-y[i])^2/(2*sigma2^2))*pi2/sigma2)
                pi2.vec.new[i] = exp(-(mu2-y[i])^2/(2*sigma2^2))*pi2/sigma2/(exp(-(mu1-y[i])^2/(2*sigma1^2))*pi1/sigma1+exp(-(mu2-y[i])^2/(2*sigma2^2))*pi2/sigma2)
                  }
         #   deltavec = c()
            while (!is.na(delta) && delta>cutvalue){              
              mu1.new = sum(pi1.vec.new[!is.na(pi1.vec.new)]*y[!is.na(pi1.vec.new)])/sum(pi1.vec.new[!is.na(pi1.vec.new)])
              sigma1.new = sqrt(sum(pi1.vec.new[!is.na(pi1.vec.new)]*(y[!is.na(pi1.vec.new)]-mu1.new)^2)/sum(pi1.vec.new[!is.na(pi1.vec.new)]))
              mu2.new = sum(pi2.vec.new[!is.na(pi2.vec.new)]*y[!is.na(pi2.vec.new)])/sum(pi2.vec.new[!is.na(pi2.vec.new)])
              sigma2.new = sqrt(sum(pi2.vec.new[!is.na(pi2.vec.new)]*(y[!is.na(pi2.vec.new)]-mu2.new)^2)/sum(pi2.vec.new[!is.na(pi2.vec.new)]))
              delta = (mu1.new-mu1)^2+(mu2.new-mu2)^2+(sigma1.new-sigma1)^2+(sigma2.new-sigma2)^2
     # Hao: added Mar. 29, 2009.  when there are too few Het snps, there may be only one point for one category, which will cause sigma to be 0 and delta to be NaN           
           #   deltavec = rbind(deltavec, c(delta, mu1.new, mu2.new))
              if (is.na(delta)){
                break
              }else{
                mu1 = mu1.new
                mu2 = mu2.new
                sigma1 = sigma1.new
                sigma2 = sigma2.new
                pi1.vec = pi1.vec.new
                pi2.vec = pi2.vec.new
                if (delta<cutvalue){
                  break
                }
              }

              for (i in 1:ny){
                pi1.vec.new[i] = exp(-(mu1-y[i])^2/(2*sigma1^2))*pi1.vec[i]/sigma1/(exp(-(mu1-y[i])^2/(2*sigma1^2))*pi1.vec[i]/sigma1+exp(-(mu2-y[i])^2/(2*sigma2^2))*pi2.vec[i]/sigma2)
                pi2.vec.new[i] = exp(-(mu2-y[i])^2/(2*sigma2^2))*pi2.vec[i]/sigma2/(exp(-(mu1-y[i])^2/(2*sigma1^2))*pi1.vec[i]/sigma1+exp(-(mu2-y[i])^2/(2*sigma2^2))*pi2.vec[i]/sigma2)
                  }
            }
            pvalue1 = pt((mu1-mu0)/sqrt(sigma0^2/norm.length+sigma1^2/ max( (sum(pi1.vec[!is.na(pi1.vec)])-1), 1) ), (norm.length+sum(pi1.vec[!is.na(pi1.vec)])-2))
            pvalue2 = pt((mu2-mu0)/sqrt(sigma0^2/norm.length+sigma2^2/ max( (sum(pi2.vec[!is.na(pi2.vec)])-1), 1) ), (norm.length+sum(pi2.vec[!is.na(pi2.vec)])-2))
            if (pvalue1<pvalue2){
              pvaluetemp = pvalue1
              pvalue1 = pvalue2
              pvalue2 = pvaluetemp
            }
            if (pvalue1<talpha2/2){
              flag1 = -1
            }else if (pvalue1>1-talpha1/2){
              flag1 = 1
            }else{
              flag1 = 0
            }
            if (pvalue2<talpha2/2){
              flag2 = -1
            }else if (pvalue2>1-talpha1/2){
              flag2 = 1
            }else{
              flag2 = 0
            }
                            
            if (flag1 == 1){
              if (flag2 == 1){
                Type = "Gain/Gain"
              }else if(flag2 == 0 ){
                Type = "Gain/Normal"
              }else{
                Type = "Gain/Loss"
              }
            }else if(flag1 == 0){
              if (flag2 == 1){
                Type = "Gain/Normal"
              }else if(flag2 == 0){
                Type = "NA.0.0"
              }else{
                Type = "Normal/Loss"
              }
            }else{
              if (flag2 == 1){
                Type = "Gain/Loss"
              }else if(flag2 == 0){
                Type = "Normal/Loss"
              }else{
                Type = "Loss/Loss"
              }
            }
            # Hao add the compared of balanced/unbalanced G/L, Nov. 16, 2009
           # if (Type == "Gain/Loss" && GLbalance==TRUE){
#              ratio1 = 1/ratio0
#              if ( (mu1+mu2)/mu0>2-cutdev && (mu1+mu2)/mu0<2+cutdev && (mu1-mu0)/(mu0-mu2)>ratio0 && (mu1-mu0)/(mu0-mu2)<ratio1){
#                Type = "Gain/Loss_balanced"
#               }else{
#                Type = "Gain/Loss_unbalanced"
#               }
#            }  
            a = sort(c(mu1,mu2))[2]
            b = sort(c(mu1,mu2))[1]
            a = round(a,digits=3)
            b = round(b,digits=3)
            Value = paste(a,"/",b,sep="")
            SNP.begin = med$beginpt[count]
            SNP.end = med$endpt[count]
            SNP.length = SNP.end-SNP.begin+1
            SNP.percent = SNP.length/N*100
            SNP.percent = round(SNP.percent,digits=3)
            Pos.begin = cnvobj$pos[SNP.begin]
            Pos.end = cnvobj$pos[SNP.end]
            Pos.length = Pos.end-Pos.begin+1
            Pos.percent = Pos.length/Effective.pos*100
            Pos.percent = round(Pos.percent,digits=3)
            temp1 = c(samplename, chrid, Type, Value, mu0, SNP.begin, SNP.end,  Pos.begin, Pos.end,  c(pvalue1,pvalue2))
            names(temp1)=c("Sample", "Chr", "Type","Value", "Normal.copy", "SNP.begin", "SNP.end", "Pos.begin", "Pos.end","major.pvalue","minor.pvalue")
            longtemp = rbind(longtemp,temp1)
          }
        }
      }
  }
    longlist = rbind(longlist,longtemp)
  }
  return(longlist)
}

