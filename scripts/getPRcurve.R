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

getPRcurve <-
function(dataFileNames,L,dataFileName_ref,L_ref,trueDataFile,L_true,methodNames,colors,linewidths,filename_out,ind_all_common) {

	source("getPrecisionRecall.R")
	source("getCommonIndices.R")
	setwd("results/")

	AP=matrix(0,(length(dataFileNames)+1),1)
	true_vector=matrix(0,L_ref,1)
	true_ind=getCommonIndices(trueDataFile,L_true,dataFileName_ref,L_ref)
	true_vector[true_ind]=1

	d=read.table(dataFileName_ref,nrows=L_ref)
	unsorted_values=d$V2[ind_all_common]
	PR=getPrecisionRecall(unsorted_values,true_vector[ind_all_common])
	prec=PR$precision
	rec=PR$recall
	AP[1]=PR$AP

	library(Hmisc)
	FONTSIZE <- 10
	file_name=paste(filename_out,".pdf",sep="")
	pdf(file=file_name, width=86/25.4, height=70/25.4)
	par(ps=FONTSIZE, cex=1)
	par(mar=c(2, 2, 0, 0)+0.4)
	par(mgp=c(1.5, 0.5, 0))

	plot (c(0,1),c(0,1),type="n",xlab="Recall",ylab="Precision")

	lines(rec,prec,lty=1,col=colors[1],lwd=linewidths[1])
	
	
	for (i in 1:length(dataFileNames)) {

		filename=dataFileNames[[i]]
		ind_in_file=getCommonIndices(filename,L[i],dataFileName_ref,L_ref)
		selected_ind_in_file=which(ind_in_file %in% ind_all_common)

		d=read.table(filename,nrows=L[i])
		unsorted_values=d$V2[selected_ind_in_file]

		PR=getPrecisionRecall(unsorted_values,true_vector[ind_all_common])
		prec=PR$precision
		rec=PR$recall
		AP[i+1]=PR$AP

		lines(rec,prec,lty=1,col=colors[i+1],lwd=linewidths[i+1])
	}

	txt=c()
	for (i in 1:(length(dataFileNames)+1)) {
		txt=append(txt,paste(methodNames[i],", AP = ",round(AP[i],digits=3),sep=""))
	}

	legend('topright',txt,lty=1,col=colors[seq(1,(length(dataFileNames)+1))],lwd=linewidths[seq(1,(length(dataFileNames)+1))])
	dev.off()

	return(AP)

}
	


