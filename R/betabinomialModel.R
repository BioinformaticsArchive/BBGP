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
function(counts1,counts2,x) {

	if (is.unsorted(x)==TRUE) {	
		order_ind=order(x)
		counts1=as.matrix(counts1[order_ind])
		counts2=as.matrix(counts2[order_ind])
		x=as.matrix(x[order_ind])
	}

	R=as.matrix(as.vector(table(x))) # number of replicates corresponding to time points in x_1

	seq_dept=counts1+counts2
	rising_allele_counts=counts1
        if (mean(head((rising_allele_counts/seq_dept),R[1]))>mean(tail((rising_allele_counts/seq_dept),tail(R,1)))) { 
           rising_allele_counts=seq_dept-rising_allele_counts # adjustment for choosing the rising allele
        }	

	alpha=1;
	beta=1;

 	y=(alpha+rising_allele_counts)/(alpha+beta+seq_dept)
        v=((alpha+rising_allele_counts)*(1+seq_dept-rising_allele_counts))/((alpha+beta+seq_dept)^2*(alpha+beta+seq_dept+1))
        y=as.matrix(y)
	v=as.matrix(v)

	bb_list=list("posteriorMean"=y, "posteriorVariance"=v, "timeVector"=x)

	return(bb_list)

}
