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

biasKernParamInit <-
function (kern) {

        kern$variance=exp(-2)
        kern$nParams=1
        kern$transforms=list(list(index=c(1),type="positive"))
        kern$isStationary=TRUE
	kern$paramNames="variance"

        return(kern)
}


biasKernExtractParam <-
function (kern,only.values=TRUE,untransformed.values=TRUE) {

        params=c(kern$variance)
        if (!only.values)
        names(params)=c("variance")

        return(params)
}

biasKernExpandParam <-
function (kern, params) {

        kern$variance=params[1]

        return(kern)
}

biasKernCompute <-
function (kern, x, x2=NULL) {
        if ( nargs() < 3 ) {
                dims <- c(dim(as.array(x))[1], dim(as.array(x))[1])
        } else {
                dims <- c(dim(as.array(x))[1], dim(as.array(x2))[1])
        }

        k=matrix(kern$variance,nrow=dims[1],ncol=dims[2])


        return(k)
}



biasKernGradient <-
function (kern,x,x2,covGrad) {

        if ( nargs()==3 ) {
                covGrad <- x2
        }

        g=matrix(sum(covGrad),1,1)

        return(g)
}

biasKernDiagCompute <-
function (kern, x) {
        k <- matrix(kern$variance, dim(as.array(x))[1], 1)
        return (k)
}
