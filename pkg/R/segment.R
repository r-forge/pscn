segment <-
function(cnvobj, verbose=FALSE, regroup.MIN=15, regroup.percent = 0.15,...){

    if(verbose) cat("\nSegmenting ",cnvobj$label,"...\n",sep="")
    
    cnvobj=segment.sub(cnvobj, verbose, start.from.seg=FALSE,...)
    cnvobj = suppressWarnings(regroup.theta(cnvobj, level0=level0, level1=level1, percent=regroup.percent, regroup.MIN=regroup.MIN))
    while(TRUE){
        if(verbose) cat("Round 1 error checking...")
        old.cp = cnvobj$chpts
        cnvobj = segment.sub(cnvobj, verbose, start.from.seg=TRUE,...)
        cnvobj = suppressWarnings(regroup.theta(cnvobj, level0=level0, level1=level1, percent=regroup.percent))
        new.cp = cnvobj$chpts
        if(length(old.cp)-length(new.cp)==0) break        
 # Hao 09/09/2008
        if(is.null(cnvobj$chpts)) break
    }
       
   if(verbose) cat("Done.\n\n")
    cnvobj
}
