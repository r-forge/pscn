pscnlist <-
function(samplename, chrs=1:23, MIN.SNPS=20, alpha1=0.01, alpha2=0.01, FWalpha=0.01, threshold=0.18, hard.threshold=0.125, GLBalance = TRUE, ratio0=0.8, cutdev=0.4, ...){
  longlist = cnvlist4(samplename=samplename, chrs=chrs, MIN.SNPS=MIN.SNPS, talpha1=alpha1, talpha2 = alpha2, ...)
  write.table(longlist, file=paste(samplename,".longlist.txt",sep=""),row.names=FALSE, sep="\t", quote=FALSE)
  mydata = read.table(paste(samplename,".longlist.txt",sep=""),header=TRUE, sep="\t")
  longlist2 = cnvlist.update(mydata, alpha=FWalpha, threshold=threshold, hard.threshold=hard.threshold, ...)
  longid = which(!is.na(longlist2$Type))
  longlist2 = longlist2[longid,]
  write.table(longlist2, file=paste(samplename,".longlist.update0.txt",sep=""),row.names=FALSE, sep="\t", quote=FALSE)
  mydata2 = read.table(paste(samplename,".longlist.update0.txt",sep=""),header=TRUE, sep="\t")
  mydata2 = mydata2[which(!is.na(mydata2$Type)),]
  shortlist = c()
  for (i in chrs){
    id = which(mydata2$Chr==i)
    if (length(id)>0){
      chrdata = mydata2[id,]
      if (length(id)==1){
        temp = chrdata
        temp$Chr = as.character(temp$Chr)
        temp$SNP.begin = as.character(temp$SNP.begin)
        temp$SNP.end = as.character(temp$SNP.end)
        temp$Pos.begin = as.character(temp$Pos.begin)
        temp$Pos.end = as.character(temp$Pos.end)
        temp$major.pvalue = as.character(temp$major.pvalue)
        temp$minor.pvalue = as.character(temp$minor.pvalue)
        shortlist = rbind(shortlist, temp )
      }else{
        chrdata2 = cnvlist.combine(chrdata, chrid=i, samplename=samplename, chrs=chrs, MIN.SNPS=MIN.SNPS, ...)
        shortlist = rbind(shortlist,chrdata2)
      }
    }
  }
  write.table(shortlist, file=paste(samplename,".shortlist.txt",sep=""),row.names=FALSE, sep="\t", quote=FALSE)
  count = 1
  while(dim(shortlist)[1]<dim(mydata)[1]){
    mydata = read.table(paste(samplename,".shortlist.txt",sep=""),header=TRUE, sep="\t")
    longlist2 = cnvlist.update(mydata, alpha=FWalpha, threshold=threshold, hard.threshold=hard.threshold, ...)
    longid = which(!is.na(longlist2$Type))
    longlist2 = longlist2[longid,]
    write.table(longlist2, file=paste(samplename,".longlist.update",count,".txt",sep=""),row.names=FALSE, sep="\t", quote=FALSE)
    mydata2 = read.table(paste(samplename,".longlist.update",count,".txt",sep=""),header=TRUE, sep="\t")
    mydata2 = mydata2[which(!is.na(mydata2$Type)),]
    shortlist = c()
    for (i in chrs){
      id = which(mydata2$Chr==i)
      if (length(id)>0){
        chrdata = mydata2[id,]
        if (length(id)==1){
          temp = chrdata
          temp$Chr = as.character(temp$Chr)
          temp$SNP.begin = as.character(temp$SNP.begin)
          temp$SNP.end = as.character(temp$SNP.end)
          temp$Pos.begin = as.character(temp$Pos.begin)
          temp$Pos.end = as.character(temp$Pos.end)
          temp$major.pvalue = as.character(temp$major.pvalue)
          temp$minor.pvalue = as.character(temp$minor.pvalue)
          shortlist = rbind(shortlist, temp )
        }else{
          chrdata2 = cnvlist.combine(chrdata, chrid=i, samplename=samplename, chrs=chrs, MIN.SNPS=MIN.SNPS, ...)
        shortlist = rbind(shortlist,chrdata2)
      }
    }
  }
  write.table(shortlist, file=paste(samplename,".shortlist.txt",sep=""),row.names=FALSE, sep="\t", quote=FALSE)
  count = count + 1
 }
 if (GLBalance == TRUE){
  shortlist2 = read.table(file=paste(samplename,".shortlist.txt",sep=""), header=TRUE)
  mu0 = shortlist2$Normal.copy[1]
  ratio1 = 1/ratio0
  for (k in 1:dim(shortlist2)[1]){
    shortlist2$Type = as.character(shortlist2$Type)
    if (shortlist2$Type[k]=="Gain/Loss"){
      M = as.numeric(strsplit(as.character(shortlist2$Value[k]),"/")[[1]][1])
      m = as.numeric(strsplit(as.character(shortlist2$Value[k]),"/")[[1]][2])      
      if ( (M+m)/mu0>2-cutdev && (M+m)/mu0<2+cutdev && (M-mu0)/(mu0-m)>ratio0 && (M-mu0)/(mu0-m)<ratio1){
        shortlist2$Type[k] = "Gain/Loss_balanced"
      }else{
        shortlist2$Type[k] = "Gain/Loss_unbalanced"
      }
    }
  }
  write.table(shortlist2, file=paste(samplename,".shortlist2.txt",sep=""),row.names=FALSE, sep="\t", quote=FALSE)       
 }
}
