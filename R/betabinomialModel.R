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

betabinomialModel <-
function(counts,seq_depth,x,rising=1) {

	# Remove the SNPs with zero coverage, i.e. zero sequencing depth:
	ind_nonzero=which(seq_depth!=0)
	counts=as.matrix(counts[ind_nonzero])
	seq_depth=as.matrix(seq_depth[ind_nonzero])
	x=as.matrix(x[ind_nonzero])
	#

	if (is.unsorted(x)==TRUE) {	
		order_ind=order(x)
		counts=as.matrix(counts[order_ind])
		seq_depth=as.matrix(seq_depth[order_ind])
		x=as.matrix(x[order_ind])
	}

	rising_allele_counts=counts

	if (rising==1) {
		R=as.matrix(as.vector(table(x))) # number of replicates corresponding to time points in x

	        if (mean(head((rising_allele_counts/seq_depth),R[1]))>mean(tail((rising_allele_counts/seq_depth),tail(R,1)))) { 
	           rising_allele_counts=seq_depth-rising_allele_counts # adjustment for choosing the rising allele
	        }	
	}

	alpha=1;
	beta=1;

 	y=(alpha+rising_allele_counts)/(alpha+beta+seq_depth)
        v=((alpha+rising_allele_counts)*(1+seq_depth-rising_allele_counts))/((alpha+beta+seq_depth)^2*(alpha+beta+seq_depth+1))
        y=as.matrix(y)
	v=as.matrix(v)

	bb_list=list("posteriorMean"=y, "posteriorVariance"=v, "timeVector"=x)

	return(bb_list)

}
