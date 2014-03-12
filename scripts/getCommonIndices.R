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

getCommonIndices <-
function(dataFileNames,L,dataFileName_ref,L_ref) {

	# dataFileNames=c("dataFileName1","dataFileName2","dataFileName3","dataFileName4",...)
	# L=c(L1,L2,L3,L4,...)
	# dataFileName_ref=
	# L_ref=

	d_ref=read.table(dataFileName_ref,nrows=L_ref)
	SNP_ref=as.matrix(d_ref$V1)
	commonIndices=as.matrix(seq(1,L_ref))

	for (i in 1:length(dataFileNames)) {

		filename=dataFileNames[i]
		d=read.table(filename,nrows=L[i])
		SNP_in_file=as.matrix(d$V1)
		commonind_in_ref=which(SNP_ref %in% SNP_in_file)
		commonIndices=as.matrix(intersect(commonIndices,commonind_in_ref))

	}

	return(commonIndices)

}
	


