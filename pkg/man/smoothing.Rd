\name{smoothing}
\alias{smoothing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{BCMIX Smoothing for allele-specific DNA copy number data.}
\description{
This function uses Bounded Complexity MIXture (BCMIX) algorithm to estimate parent-specify copy numbers.
}
\usage{
smoothing(samplename, inputdata,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samplename}{A name indicates where the data is sampled from, eg. patient ID}
  \item{inputdata}{A file with four columns "Chromosome ID" (eg. any number from 1 to 23), "SNP Position (bp)", "logR" and "Bfrequency". The column names should be "Chr", "Position", "logR" and "bfreq". }
  \item{platform}{Specifies where is data is from, eg. Illumina. (Default="Unknown")}
  \item{sm.plot}{If TRUE, then the interim results will be plotted.  (Default=FALSE) }
  \item{hyper}{A txt file specifies initial hyper parameters, if null (default), then the initial hyperparameters are set using moment estimators.}
  \item{mixstate}{A txt file specifying initial states of SNPs, if null (default), then the initial states are generated via K-means clustering ignoring position information. }
  \item{niters}{Specifies how many iterations to do for hyperparameter estimation. (Default=2)}
  \item{nBIMIXiters}{Specifies the number of iterations for BIMIX.  (Default=10) }
  \item{K}{'M' parameter for BCMIX.}
  \item{M}{'K' parameter for BCMIX.}
   \item{selectHyper}{Which hyperparameters to update via EM?  If "1", then update only the transition matrix.  If "2", then update everything but the transition matrix.  If "3", then update all hyperparameters.  If "4", then update only the error covariance matrices.  (Default=4)}
  \item{verbose}{if TRUE (default), print out process during the run }
  \item{save.interim}{if TRUE, interim cnvobj will be saved as .Rdata file.  (Default=FALSE)}
}
\details{
    This function takes 2-3 minutes for a typical file containing logR and Bfrequency of SNP data from one human choromosome.
}
\value{
It saves files samplename.Chr1.Smooth.Rdata etc.
}
\author{
Hao Chen, Haipeng Xing, and Nancy R. Zhang

Maintainer:  Hao Chen (haochen@stanford.edu)
}
\references{
Chen, H., Xing, H. and Zhang, N.R. (2009) Estimation of parent specific DNA copy number
in tumors using high-density genotyping arrays.  Under Review, manuscript can be downloaded from www-stat.stanford.edu/~nzhang/.

}
\seealso{
this function follows with \code{\link{segmentation}} and \code{\link{pscnlist}}. 
\code{\link{pscn.plot}} is the plot function to view a variety of results.  
}
\examples{
data(ex0)  # this file contain four columns "Chromosome ID", "SNP Position (bp)", "logR" and "Bfrequency" with column names "Chr", "Position", "logR" and "bfreq". 
# smoothing(samplename="TCGA.02.0043.01A.01D", inputdata=mydata, platform="Illumina")  # to excute, remove the "#" at the begining
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

