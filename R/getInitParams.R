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

getInitParams <-
function (model,sigma2n) {

	if (nargs()<2) {
		grid_size=10
		sigma2n_search_range=seq(-20,20,length=grid_size)
		LogLik=matrix(0,grid_size)

		for(i in 1:grid_size)
		{
			sigma2n=sigma2n_search_range[i]
			model1 = gpExpandParam(model, c(sigma2n))
			LogLik[i] = gpLogLikelihood(model1)
		}

		ind_maxLogLik=which.max(LogLik)
		sigma2n_updated=sigma2n_search_range[ind_maxLogLik]
		initial_params=c(sigma2n_updated)
	}
	else {
		grid_size=5
		iw_bound= model$kern$comp[[1]]$options$inverseWidthBounds[2]
		iw_search_range=seq((log(0.01)),(log(iw_bound)),length=grid_size)
		sigma2f_search_range=seq(-10,10,length=grid_size)
		sigma2n_search_range=seq((log(sigma2n)-5),(log(sigma2n)+5),length=grid_size)
		LogLik=matrix()
		I=matrix()
		K=matrix()
		J=matrix()
		vec_size=grid_size^3;
		for(k in 1:grid_size)
			{
			for(j in 1:grid_size)
			{	
				for(i in 1:grid_size)
				{
					
					iw=iw_search_range[k]
					sigma2f=sigma2f_search_range[j]
					sigma2n=sigma2n_search_range[i]
					model1 = gpExpandParam(model, c(iw, sigma2f,sigma2n))
					LogLik = rbind(LogLik,gpLogLikelihood(model1))
					K=rbind(K,k)
					J=rbind(J,j)
					I=rbind(I,i)
				
				}
			}
		}	

		LogL=LogLik[2:(vec_size+1)]
		ind_maxLogLik=which.max(LogL)
		i=I[(ind_maxLogLik+1)]
		j=J[(ind_maxLogLik+1)]
		k=K[(ind_maxLogLik+1)]

		iw_updated=iw_search_range[k]	
		sigma2f_updated=sigma2f_search_range[j]	
		sigma2n_updated=sigma2n_search_range[i]
	
		initial_params=c(iw_updated,sigma2f_updated,sigma2n_updated)
	}

	return(initial_params)

}

