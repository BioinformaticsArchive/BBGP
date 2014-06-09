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
function(dataFileName,start_line,end_line,NoHeaderLines=1,NoInfoColumns=3,NoOptions=6,sep=":") {

	headerLines=read.table(dataFileName,skip=0,nrows=NoHeaderLines)
	X=as.matrix(as.numeric(headerLines[1,][,-seq(1,NoInfoColumns)]))
 	
	noLines=end_line-start_line+1
	noSkip=start_line+NoHeaderLines-1

	sample_data=read.table(dataFileName,skip=noSkip,nrows=noLines)

	J=ncol(sample_data)-NoInfoColumns # number of (time,rep) combinations in the data file.

	ID=as.matrix(paste(as.character(sample_data[,1]),"_",as.character(sample_data[,2]),sep=""))
	#REF_ALLELES=as.matrix(as.character(sample_data[,3]))
	COUNTS=as.matrix(sample_data[,-seq(1,NoInfoColumns)])

	counts=matrix(nrow=noLines,ncol=J)
	seq_depth=matrix(nrow=noLines,ncol=J)

	for (i in 1:noLines) {
		d=COUNTS[i,]
		countsMatrix=matrix(nrow=NoOptions,ncol=J)
		for (j in 1:J) {
			d1=d[j]
			s=as.matrix(as.numeric(unlist(strsplit(d1,sep))))
			countsMatrix[,j]=s
		}
#		if (0 %in% colSums(countsMatrix)) {
#			print(sprintf("None of the alleles have been sequenced for the SNP on line %d, which results in zero sequencing depth. Check the data file.", i+start_line-1))
#		} else {
			ind_nonzero=unique(which(countsMatrix!=0,arr.ind=TRUE)[,1])
			if (length(ind_nonzero)>2) {
				print(sprintf("SNP on line %d is not bi-allelic. Check the data file.", i+start_line-1))
			} else {
				allele_counts=countsMatrix[ind_nonzero,]
				counts[i,]=allele_counts[1,]
				if (length(ind_nonzero)<2) {
					seq_depth[i,]=counts[i,]+0	
				} else {
					seq_depth[i,]=counts[i,]+allele_counts[2,]
				}
			} 
		}
#	}

	snpData=list("ID"=ID,"counts"=counts,"seq_depth"=seq_depth,"timeVector"=X)
	return(snpData)

}


