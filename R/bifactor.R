
#' approximate negtive log-likelihood for bi-factor
#'@description approximate negtive loglikielihood for bi-factor copula model with
#'all linking copulas Frank
#'@param udata input data in the u-scale. N*D
#'@param th parameters in the bi-factor model, vector of 2*D length(if all the
#'linking copulas are one-parameters copula)
#'@param vlat proxies of the latent
#'@param grsize group size in the model
#'@return nllk negative log-likelihood
#'@export
f90bifctFrkProxynllk=function(udata,th,vlat,grsize){
	n=nrow(udata)
	npar=length(th)
	mgrp=length(grsize)
	dvar=sum(grsize)
	#if(!is.loaded("frkproxynllk"))  dyn.load("./libs/FactorModels.so")
	out=.Fortran("frkproxynllk",
							 as.integer(npar),as.double(th),as.integer(mgrp),as.integer(n),as.integer(dvar),
							 as.integer(grsize),as.double(udata),as.double(vlat),
							 nllk=as.double(0.),PACKAGE = "FactorModels")
	out$nllk
}










#' approximate negative log-likelihood for bi-factor
#'@description approximate negtive loglikielihood for bi-factor copula model with
#'all linking copulas Gumbel
#'@param udata input data in the u-scale. N*D
#'@param th parameters in the bi-factor model, vector of 2*D length(if all the
#'linking copulas are one-parameters copula)
#'@param vlat proxies of the latent variables in the model
#'@param grsize group size in the model
#'@return nllk negative log-likelihood
#'@export
#'
f90bifctGumbelProxynllk=function(udata,th,vlat,grsize){
	n=nrow(udata)
	npar=length(th)
	mgrp=length(grsize)
	dvar=sum(grsize)
	#if(!is.loaded("gumproxynllk"))  dyn.load("./libs/FactorModels.so")
	out=.Fortran("gumproxynllk",
							 as.integer(npar),as.double(th),as.integer(mgrp),as.integer(n),as.integer(dvar),
							 as.integer(grsize),as.double(udata),as.double(vlat),
							 nllk=as.double(0.),PACKAGE = "FactorModels")
	out$nllk
}







#' compute new proxies for bi-factor copula models
#'@description latent update for the bi-factor copulas with all linking
#'copulas Frank
#'@param udata input data in u-scale N*D
#'@param th parameters for linking copulas in the model
#'@param grsize group size in the bi-factor model
#'@param nq number of Gaussian points
#'@param xl Gaussian-Legendre points
#'@param wl Gaussian-Legendre weights
#'@return v0: proxies of the global latent variable
#'@return vg: proxies of the local latent variables
#'@export
latUpdateBifctFrk=function(udata,th,grsize,nq,xl,wl){
	n=nrow(udata)
	npar=length(th)
	mgrp=length(grsize)
	dvar=sum(grsize)
	nq=length(xl)
	npar=2*dvar
	#if(!is.loaded("latupdatebifact"))  dyn.load("./libs/FactorModels.so")
	out=.Fortran("latupdatebifact",
							 as.integer(npar),as.double(th),as.integer(mgrp),as.integer(dvar),as.integer(n),
							 as.integer(grsize),as.double(udata),as.integer(nq),as.double(xl),as.double(wl),
							 v0mat=rep(0,n),
							 vgmat=matrix(0.0,nrow=n,ncol=mgrp),PACKAGE = "FactorModels")

	return(list(v0=out$v0mat,vg=out$vgmat))
}






#'compute new proxies for bi-factor copula models
#'@description latent update for the bi-factor copulas with all linking
#'copulas Gumbel
#'@param udata input data in u-scale N*D
#'@param th parameters for linking copulas in the model
#'@param grsize group size in the bi-factor model
#'@param nq number of  Gaussian-Legendre points
#'@param xl  Gaussian-Legendre points
#'@param wl  Gaussian-Legendre weights
#'@return v0: proxies of the global latent variable
#'@return vg: proxies of the local latent variables
#'
#'@export
#'
latUpdateBifctGumbel=function(udata,th,grsize,nq,xl,wl){
	n=nrow(udata)
	npar=length(th)
	mgrp=length(grsize)
	dvar=sum(grsize)
	nq=length(xl)
	npar=2*dvar
	#if(!is.loaded("latupdatebifact2"))  dyn.load("./libs/FactorModels.so")
	out=.Fortran("latupdatebifact2",
							 as.integer(npar),as.double(th),as.integer(mgrp),as.integer(dvar),as.integer(n),
							 as.integer(grsize),as.double(udata),as.integer(nq),as.double(xl),as.double(wl),
							 v0mat=rep(0,n),
							 vgmat=matrix(0.0,nrow=n,ncol=mgrp),PACKAGE = "FactorModels")

	return(list(v0=out$v0mat,vg=out$vgmat))
}
