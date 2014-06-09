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
function(timePoints="all",dataFileName,start_line,end_line,resultFileName,plots_path,thr=10) {

####################################################################################################
##  											          ##
## dataFileName: Name of the data file which contains the counts for bi-allelic SNPs.             ##
## L: number of lines in data file, including the header line.                                    ##
## timePoints: Vector containing the time points which will be used in GP models.                 ##
## thr: Threshold for the BF such that if BF>thr the plot will be created for the model fit.      ##
##                                                                                                ##
##                                                                                                ##
## Example usage for the sampleData:                                                              ##   
##                                                                                                ##
## > dataFileName="data/sampleDataCounts"                                                         ##
## > resultFileName="results/BBGP_summary"                                                        ##
## > plots_path=paste(getwd(),"/plots/",sep="")                                                   ##
##                                                                                                ##
## > runSample(c(0,14,22,28,38,50,60),dataFileName,1,5,resultFileName)                            ##
## For also getting the GP model fit plots if Bayes Factor > thr :                                ##
## > runSample(c(0,14,22,28,38,50,60),dataFileName,1,5,resultFileName,plots_path,thr)             ##
##                                                                                                ## 
####################################################################################################

	# Install gptk package:
	# install.packages("gptk") 
	library("gptk")
	# load BBGP functions:
	source("loadBBGP.R")
	loadBBGP()

	snpData=readCountsData(dataFileName,start_line,end_line)
	COUNTS=snpData$counts
	SEQ_DEPTH=snpData$seq_depth
	SNP_ID=snpData$ID
	X=snpData$timeVector
 	N=length(SNP_ID) # number of SNPs in the data file

	BayesFactors=matrix(0,N,1)	
	SNP=matrix("",N,1)

	if (timePoints[1]!="all") {
		ind_selectedTimePoints=which(X %in% timePoints)
		X=as.matrix(X[ind_selectedTimePoints])
		COUNTS=as.matrix(COUNTS[,ind_selectedTimePoints])
		SEQ_DEPTH=as.matrix(SEQ_DEPTH[,ind_selectedTimePoints])
	}


	for (i in 1:N) {

		counts=as.matrix(COUNTS[i,])
		seq_depth=as.matrix(SEQ_DEPTH[i,])

		bb_model=betabinomialModel(counts,seq_depth,X)
		x=bb_model$timeVector
		y=bb_model$posteriorMean
		v=bb_model$posteriorVariance

                if ((range(y)[2]-range(y)[1])==0) {
			BayesFactors[i]=NA
		} else {
			indModelCovTypes=c("white","fixedvariance")
			depModelCovTypes=c("rbf","white","fixedvariance")
			rslt=bbgp_test(x,y,v,indModelCovTypes,depModelCovTypes)
			BayesFactors[i]=rslt$BF
                }

		SNP[i]=SNP_ID[i]

		if (nargs()>5) {
			if (BayesFactors[i] > thr) {
				model0=rslt$independentModel
				model1=rslt$dependentModel
				plot_name=paste(plots_path,SNP_ID[i],"_model0",".pdf",sep="")
				modelfitPlot(model0,plot_name)
				plot_name=paste(plots_path,SNP_ID[i],"_model1",".pdf",sep="")
				modelfitPlot(model1,plot_name)
			}	
		}
	
	}

	d=data.frame(SNP,BayesFactors)
	names(d)=c("SNP_ID","Bayes Factor")
	writeOutputFile(resultFileName,d)


}
