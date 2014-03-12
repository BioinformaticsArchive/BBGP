# Copyright (c) 2014, Hande TOPA
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
# 
#     Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in
#     the documentation and/or other materials provided with the
#     distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

runSample <-
function(dataFileName,L,timePoints,thr,pl=0) {

#################################################################################################
##  											       ##
## dataFileName: Name of the data file which contains the counts for bi-allelic SNPs.          ##
## L: number of lines in data file, including the header line.                                 ##
## timePoints: Vector containing the time points which will be used in GP models.              ##
## thr: Threshold for the BF such that if BF>thr the plot will be created for the model fit.   ##
##                                                                                             ##
##                                                                                             ##
## Examples usage for the sampleData:                                                          ##
##                                                                                             ##
## > runSample('sampleDataCounts',6,c(0,14,22,28,38,50,60))                                    ##
## For also getting the GP model fit plots if Bayes Factor > thr :                             ##
## > runSample('sampleDataCounts',6,c(0,14,22,28,38,50,60),thr)                                ##
##                                                                                             ## 
#################################################################################################

	# Install gptk package:

	install.packages("gptk") 
	library("gptk")

	#gptk_path=paste(getwd(),"/gptk/",sep="")  
	#install.packages(pkgs="gptk_1.08.tar.gz",lib=gptk_path,repos=NULL,type="source") # Bug report was sent for the original gptk package. Until it is fixed, let's use our 'fixed' gptk package locally.
	#library('gptk',lib.loc=gptk_path)

	source("fixedvarianceKernFunctions.R")
	source("calculateLbound.R")
	source("betabinomialModel.R")
	source("bbgp_test.R")
	source("constructModel.R")
	source("getInitParams.R")
	source("modelfitPlot.R")
	source("readCountsData.R")

	plots_path=paste(getwd(),"/plots/",sep="")
	results_path=paste(getwd(),"/results/",sep="")
	
	if (nargs()>3) {
		pl=1
	}

	setwd("data/")
	snpData=readCountsData(dataFileName,L) 
	COUNTS1=snpData$allele1_counts
	COUNTS2=snpData$allele2_counts
	SNP_ID=snpData$SNP_ID
	X=snpData$timeVector
 	N=length(SNP_ID) # number of SNPs

	BayesFactors=matrix(0,N,1)	
	SNP=matrix("",N,1)

	ind_selectedTimePoints=which(X %in% timePoints)
	X=as.matrix(X[ind_selectedTimePoints])



	for (i in 1:N) {

		counts1=as.matrix(COUNTS1[i,ind_selectedTimePoints])
		counts2=as.matrix(COUNTS2[i,ind_selectedTimePoints])
		bb_model=betabinomialModel(counts1,counts2,X)
		x=bb_model$timeVector
		y=bb_model$posteriorMean
		v=bb_model$posteriorVariance

		rslt=bbgp_test(x,y,v)
		BayesFactors[i]=rslt$BF
		SNP[i]=SNP_ID[i]


		if (pl==1) {
			if (BayesFactors[i] > thr) {
				model0=rslt$independentModel
				model1=rslt$dependentModel
				setwd(plots_path)
				modelfitPlot(model0,SNP_ID[i])
				modelfitPlot(model1,SNP_ID[i])
			}	
		}
	
	}

	d=data.frame(SNP,BayesFactors)
	names(d)=c("SNP_ID","Bayes Factor")
	filename=paste(dataFileName,"_BBGPsummary",sep="")
	setwd(results_path)
	write.table(d,file=filename,quote=F,sep='\t',col.names=FALSE,row.names=FALSE)

}
