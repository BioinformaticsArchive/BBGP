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

readCountsData <-
function(dataFileName,L,SNPinfoColumns=3) {

	# L: number of lines in the data file, including the header line.
 
	sample_data=read.table(dataFileName,nrows=L)
	X=as.matrix(sample_data[1,][,-seq(1,SNPinfoColumns)])
	X=as.matrix(as.numeric(X[1,]))
	sample_data=sample_data[-1,]	
	J=ncol(sample_data)-SNPinfoColumns # number of (time,rep) combinations in the data file.
	SNP_ID=as.matrix(paste(as.character(sample_data[,1]),"_",as.character(sample_data[,2]),sep=""))
	REF_ALLELES=as.matrix(as.character(sample_data[,3]))
	COUNTS=as.matrix(sample_data[,-seq(1,SNPinfoColumns)])

	allele1_counts=matrix(nrow=L-1,ncol=J)
	allele2_counts=matrix(nrow=L-1,ncol=J)
	for (i in 1:(L-1)) {
		d=COUNTS[i,]
		countsMatrix=matrix(nrow=6,ncol=J)
		for (j in 1:J) {
			d1=d[j]
			s=as.matrix(as.numeric(unlist(strsplit(d1,":"))))
			countsMatrix[,j]=s
		}
		if (0 %in% colSums(countsMatrix)) {
			print(sprintf("None of the alleles have been sequenced for the SNP on line %d, which results in zero sequencing depth. Check the data file.", i+1))
		} else {
			ind_nonzero=unique(which(countsMatrix!=0,arr.ind=TRUE)[,1])
			if (length(ind_nonzero)>2) {
				print(sprintf("SNP on line %d is not bi-allelic. Check the data file.", i+1))
			} else {
				allele_counts=countsMatrix[ind_nonzero,]
				allele1_counts[i,]=allele_counts[1,]
				if (length(ind_nonzero)<2) {
					allele2_counts[i,]=0	
				} else {
					allele2_counts[i,]=allele_counts[2,]
				}
			} 
		}
	}

	snpData=list("SNP_ID"=SNP_ID,"allele1_counts"=allele1_counts,"allele2_counts"=allele2_counts,"timeVector"=X)
	return(snpData)

}

