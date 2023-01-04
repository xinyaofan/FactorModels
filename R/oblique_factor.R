
#' approximate nllk for oblique factor copula model
#' @description This function computes the approxiamte nllk
#' for oblique-factor copula model with the proxies plugged in
#' @param param parameters in the oblique-factor copula model
#' @param dstruct list of input parameters:
#' @param iprfn print gradient and hessian if true
#' @return nllk: approximate nllk; grad: gradient of nllk;
#' hess: hessian matrix of nllk
#' @export
f90str1proxynllk<-function(param,dstruct,iprfn=FALSE){
	grsize=dstruct$grsize
	mgrp=length(grsize)
	udata=dstruct$udata
	lat=dstruct$lat
	th=param
	xl=dstruct$xl
	wl=dstruct$wl
	nq=length(wl)
	n=nrow(udata)
	dvar=ncol(udata)
	npar=dvar+mgrp
	#if(!is.loaded("strfrkgum1"))  dyn.load("./mylib.so")
	out=.Fortran("strfrkgum1",as.integer(npar),
		     as.double(th),as.integer(mgrp),
		     as.integer(n),as.integer(dvar), 
		     as.integer(grsize),
		     as.double(udata),as.double(lat),
		     as.integer(nq),as.double(wl),
		     as.double(xl),
		     nllk=as.double(0.),lgrad=rep(0.0,npar),
		     lhess=matrix(0.0,npar,npar),
		     PACKAGE = "FactorModels")
	return(list("fnval"=out$nllk,"grad"=out$lgrad,"hess"=out$lhess))
}






#'parameter estimation using proxy method
#'@description parameter estimation using mean version of proxies in
#'oblique factor copula model
#'@param udata N*D data in u-scale
#'@param grsize group size in the model
#'@param fam copula families of linking copulas
#'@param start starting values for paramters
#'@param LB upper bound for parameters, length d
#'@param UB lower bound for parameters, length d
#'@param xl  Gaussian-Legendre points
#'@param wl  Gaussian-Legendre weights
#'@param iprint print the gradient and hessian if T
#'@return mlpx1 estimated parameters from mean-proxy method
#'        mlpx2 estimated parameters from new proxy method
#'        proxyMean mean verison of proxies
#'        proxyNew  new proposed proxies
#'@export
proxyOblique<-function(udata,grsize,fam,start,LB,UB,xl,wl,iprint){
	n=nrow(udata)
	d=ncol(udata)
	mgrp=length(grsize)
	npar=d+mgrp*(mgrp-1)/2
	proxyMean=matrix(NA,n,mgrp)
	proxyNew=matrix(NA,n,mgrp)
	gp=c(0,cumsum(grsize))
	for(i in 1:mgrp){
		proxyMean[,i]=uscore(apply(udata[,((gp[i]+1):gp[i+1])],1,mean))
	}
	dstruct=list(udata=udata,grsize=grsize,lat=proxyMean,xl=xl,wl=wl)
	out1=pdhessminb(param=start,f90str1proxynllk,ifixed=rep(FALSE,npar), 
			dstruct,LB=LB,UB=UB, mxiter=30, eps=5.e-5,iprint=T)
 	mlpx1=out1$parmin

	for(i in 1:mgrp){
		th_tem=mlpx1[((gp[i]+1+mgrp):(gp[i+1]+mgrp))]
		proxyNew[,i]=uscore(latUpdateOnefct(th=th_tem,udata=udata[,((gp[i]+1):gp[i+1])],nq=25,xl=xl,wl=xl,family=fam))
	}

	dstruct=list(udata=udata,grsize=grsize,lat=proxyNew,xl=xl,wl=wl)
	out2=pdhessminb(param=mlpx1,f90str1proxynllk,
			ifixed=rep(FALSE,npar), dstruct,
			LB=LB, UB=UB, mxiter=30, eps=5.e-5,iprint=T)
  	mlpx2=out2$parmin
	return(list(mlpx1=mlpx1,mlpx2=mlpx2,proxyMean=proxyMean,proxyNew=proxyNew))
}

