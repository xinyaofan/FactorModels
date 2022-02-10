#' approximate nllk for 1-factor copula model
#' @description This function computes the approxiamte nllk
#' for 1-factor copula model with the proxies plugged in
#' @param param parameters of dimension d
#' in the 1-factor copula model
#' @param dstruct list of input parameters: $lat: proxies of
#' the latent variables, $udata: observed data; $fam: copula
#' families for D linking copulas
#' @param iprint if T print gradient and Hessian
#' @return nllk: approximate nllk; grad: gradient of nllk;
#' hess: hessian matrix of nllk
#' @export
f901fproxynllk<-function(param,dstruct,iprint=F){
	udata=dstruct$udata
	lat=dstruct$lat
	fam=dstruct$fam
	n=nrow(udata)
	d=ncol(udata)
  #call fortran function
	#if(!is.loaded("onefact"))  dyn.load("./libs/FactorModels.so")
	out=.Fortran("onefact",as.double(param),as.integer(d),
							 as.integer(n),as.double(udata),as.double(lat)
							 ,as.integer(fam),
							 nllk=as.double(0.),lgrad=rep(0.0,d),
							 lhess=matrix(0.0,d,d),PACKAGE = "FactorModels")
	return(list("fnval"=out$nllk,"grad"=out$lgrad,"hess"=out$lhess))
}


#' approximate nllk for 1-factor copula model with copulas from
#' different families(2-parameter familes,like bb1, student-t)
#' @description This function computes the approxiamte nllk
#' for 1-factor copula model with the proxies plugged in
#' @param param 2*d vector (par1,par2) for j-th linking copula
#' @param dstruct list of input parameters: $lat: proxies of
#' the latent variables, $udata: observed data; $fam: copula
#' families for D linking copulas
#' @param iprint if T print gradient and Hessian
#' @return nllk: approximate nllk; grad: gradient of nllk;
#' hess: hessian matrix of nllk
#' @export
f901fproxynllk2<-function(param,dstruct,iprint=F){
	param=matrix(param,2,d)
	udata=dstruct$udata
	lat=dstruct$lat
	fam=dstruct$fam
	n=nrow(udata)
	d=ncol(udata)
	#call fortran function
	#if(!is.loaded("onefact"))  dyn.load("./libs/FactorModels.so")
	out=.Fortran("onefact2",as.double(param),as.integer(d),
							 as.integer(n),as.double(udata),as.double(lat)
							 ,as.integer(fam),
							 nllk=as.double(0.),lgrad=rep(0.0,2*d),
							 lhess=matrix(0.0,2*d,2*d),PACKAGE = "FactorModels")
	return(list("fnval"=out$nllk,"grad"=out$lgrad,"hess"=out$lhess))
}








#' compute the new proxies corrected on the mean of observations
#' @description This function computes the proposed new proxies
#' in the 1-factor copula model
#' @param th the estimated parameters in the 1-factor copula model
#' @param udata input N*D udata
#' @param nq number of  Gaussian-Legendre points
#' @param xl Gaussian quadrature points
#' @param wl Gaussian quadrature weights
#' @param family copula familes for d linking copulas
#' @return lat: proxies of the latent variable
#' @export
latUpdateOnefct<- function(th,udata,nq,xl,wl,family) {
	n=dim(udata)[1]
	d=dim(udata)[2]
	nq=length(xl)
	#if(!is.loaded("latupdate"))  dyn.load("./libs/FactorModels.so")
	out=.Fortran("latupdate",as.double(th),as.integer(n),as.integer(d),
							 as.double(udata),as.integer(nq),as.double(xl),as.double(wl),
							 as.integer(family),
							 lat=rep(0,n),PACKAGE = "FactorModels")
	return(out$lat)
}







#' compute the new proxies corrected on the mean of observations
#' @description This function computes the proposed new proxies
#' in the 1-factor copula model with different copula families
#' @param th the estimated parameters in the 1-factor copula model
#' @param udata input N*D udata
#' @param nq number of  Gaussian-Legendre points
#' @param xl Gaussian quadrature points
#' @param wl Gaussian quadrature weights
#' @param family copula familes for d linking copulas
#' @return lat: proxies of the latent variable
#' @export
latUpdateOnefct2<- function(th,udata,nq,xl,wl,family) {
	n=dim(udata)[1]
	d=dim(udata)[2]
	nq=length(xl)
	#if(!is.loaded("latupdate"))  dyn.load("./libs/FactorModels.so")
	out=.Fortran("latupdate3",as.double(th),as.integer(n),as.integer(d),
							 as.double(udata),as.integer(nq),as.double(xl),as.double(wl),
							 as.integer(family),
							 lat=rep(0,n),PACKAGE = "FactorModels")
	return(out$lat)
}










#'parameter estimation using proxy method
#'@description parameter estimation using mean version of proxies in 1-factor
#'copula model
#'@param udata N*D data in u-scale
#'@param family copula families of D linking copulas(the same), a number
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
proxy1fct<-function(udata,family,start,LB,UB,xl,wl,iprint){
	d=ncol(udata)
	proxyMean=uscore(apply(udata,1,mean)) #mean proxy

	fam=rep(family,d)
	dstruct=list(udata=udata,lat=proxyMean,fam=fam)
	out1=pdhessminb(start,f901fproxynllk,ifixed = rep(F,d),dstruct=dstruct,
									LB=LB,UB=UB,mxiter=30,eps=1.e-4,iprint=F)
	mlpx1=out1$parmin
	#new proxy
	lat_update=latUpdateOnefct(th=mlpx1,udata=udata,nq=25,xl=xl,wl=wl,family=family)
	proxyNew=uscore(lat_update)
	dstruct=list(udata=udata,lat=proxyNew,fam=fam)
	out2=pdhessminb(start,f901fproxynllk,ifixed = rep(F,d),dstruct=dstruct,
									LB=LB,UB=UB,mxiter=30,eps=1.e-4,iprint=F)
	mlpx2=out2$parmin
	return(list(mlpx1=mlpx1,mlpx2=mlpx2,proxyMean=proxyMean,proxyNew=proxyNew))
}








#'parameter estimation using proxy method
#'@description parameter estimation using mean version of proxies in 1-factor
#'copula model with different copula families
#'@param udata N*D data in u-scale
#'@param fam copula families of D linking copulas, a vector of length d
#'now only allow 2: student-t, 4:gumbel  and 5: frank copulas
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
proxy1fct2<-function(udata,fam,start,LB,UB,xl,wl,iprint){
	d=ncol(udata)
	proxyMean=uscore(apply(udata,1,mean)) #mean proxy

	dstruct=list(udata=udata,lat=proxyMean,fam=fam)
	out1=pdhessminb(start,f901fproxynllk,ifixed = rep(F,d),dstruct=dstruct,
									LB=LB,UB=UB,mxiter=30,eps=1.e-4,iprint=F)
	mlpx1=out1$parmin
	#new proxy
	lat_update=latUpdateOnefct(th=mlpx1,udata=udata,nq=25,xl=xl,wl=wl,family=fam)
	proxyNew=uscore(lat_update)
	dstruct=list(udata=udata,lat=proxyNew,fam=fam)
	out2=pdhessminb(start,f901fproxynllk,ifixed = rep(F,d),dstruct=dstruct,
									LB=LB,UB=UB,mxiter=30,eps=1.e-4,iprint=F)
	mlpx2=out2$parmin
	return(list(mlpx1=mlpx1,mlpx2=mlpx2,proxyMean=proxyMean,proxyNew=proxyNew))
}
