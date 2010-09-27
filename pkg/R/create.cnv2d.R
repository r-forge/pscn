create.cnv2d <-
function(obs, samplename, chromosome=NULL, pos=NULL, label=NULL, platform="Unknown",genotype.freq=NULL, LargestCN = 20)
{        
    dimobs=dim(obs)
    n=dimobs[1]

    if(is.null(chromosome)) chromosome = rep(1,n)                                     
    if(is.null(pos)) pos = c(1:n)

    if (n != length(chromosome) || n != length(pos))
        stop("intensity, chromosome, and pos need to have same length.")

    isvalid = which(!(is.na(obs[,1]) | is.na(obs[,2])))    
    
    n=length(isvalid)
    nchroms = length(unique(chromosome[isvalid]))

    logR = obs[isvalid,1]
    theta = obs[isvalid,2]
    if (platform == "Illumina"){  
      AB = illumina.getAB(logR,theta)
    }else{
      AB = affymetrix.getAB(logR, theta)
    }
    
    outlier = which((AB$A>LargestCN | AB$B>LargestCN))
    A = AB$A
    B = AB$B
    A[outlier] = LargestCN + rnorm(length(outlier),0,0.1)
    B[outlier] = LargestCN + rnorm(length(outlier),0,0.1)
    rawdata = list(logR = logR, theta=B/(A+B), A=A, B=B)
    intensity=cbind(A, B) 

    if(is.null(label)){
      if(is.null(chromosome)){
        label= samplename
      }else{
        label = paste(samplename,".Chr",chromosome[1],sep="")
      }
    }
        
    # Initialize genotype frequencies.
    if(is.null(genotype.freq)){
        genotype.freq = matrix(nrow=nrow(intensity),ncol=4, data=0.25)
    } else {
        nas = is.na(genotype.freq)
        rownas = rowSums(nas)
        rowhasna = which(rownas>0)
        genotype.freq[rowhasna,] = 0.25
        genotype.freq = genotype.freq[isvalid,]
    }
        
    
    cnvobj <-
        list(intensity = intensity,
             chromosome = chromosome[isvalid],
             pos= pos[isvalid],
             nchroms = nchroms,
             samplename = samplename,
             label = label,
             statemat = vector("list",4),
             rawdata = rawdata,
             platform = platform,
             genotype.freq = genotype.freq   
            )
    cnvobj$statemat[[1]] = matrix(c(2,0,0,0),2,2,byrow=TRUE)
    cnvobj$statemat[[2]] = matrix(c(1,1,1,-1),2,2,byrow=TRUE)
    cnvobj$statemat[[3]] = matrix(c(1,1,-1,1),2,2,byrow=TRUE)
    cnvobj$statemat[[4]] = matrix(c(0,2,0,0),2,2,byrow=TRUE)
 
    class(cnvobj) <- "cnv2d"
    attr(cnvobj, "call") <- sys.call()    
    inherits(cnvobj,"cnv")
    
    cnvobj$intensity.unnormed=cnvobj$intensity
    cnvobj$normparams$scale.fac=c(1,1)
    cnvobj$normparams$cstar=c(1,1)

    cat(paste("CNV object ",label," created: Number of probes:", n, ", number of chromosomes:", nchroms, ".\n",sep=""))
    
    cnvobj
}
