#' nllk for 1-factor model with weak residual dep
#' @description negative log-likelihood with gradient and Hessian computed in f90 for
#' copula from 1-factor/1-truncated vine (tree for residual dependence
#' conditional on a latent variable);
#' models included normal,t,BB1(r),BB8(r),gumbel,frank,survival gumbel
#' families for the frist tree and also the residual tree
#' @param param parameter vector (length 2*d+2*(d-1))
#' @param dstruct = list with  data set $data on (0,1)^d
#'      $edg1, $edg2 (d-1)x1 vectors for node labels of edges 1,...,d-1
#'       for the tree for residual dependence,
#'       $fam for copula family labels for the first tree
#' @param iprfn =indicator for printing of function and gradient (within NR iterations)
#' @return nllk negative log-likelihood
#' @return grad gradient of nllk
#' @return hessian matrix of nllk
#' @export
#'
f901f1tnllk=function(param,dstruct,iprfn=FALSE){
	udata=dstruct$data
	edg1=dstruct$edg1
	edg2=dstruct$edg2
	xl=dstruct$xl;
	wl=dstruct$wl;
	d=ncol(udata);
	n=nrow(udata);
	fam=dstruct$fam
	nq=length(xl)
	param=param
	npar=length(param)

  out= .Fortran("ft",as.integer(npar), 
		as.double(param),
		as.integer(d),
		as.integer(n),
		as.double(udata),
		as.integer(fam),
		as.integer(nq), 
		as.double(wl),
		as.double(xl),
		as.integer(edg1), 
		as.integer(edg2),
		nllk=as.double(0.),
		grad=as.double(rep(0,npar)),
		hess=as.double(rep(0,npar*npar)),PACKAGE = "FactorModels")
	
		nllk=out$nllk; hess=matrix(out$hess,npar,npar); grad=out$grad;
	
		list(fnval=nllk, grad=grad, hess=hess)
}






#' approximate nllk for 1-factor copula model with weak res dep
#' @description This function computes the approxiamte nllk
#' for 1-factor copula model with the proxies plugged in
#' @param param parameters in the 1-factor copula model with
#' weak res dep, vector of length 2*d-1
#' @param dstruct list with  data set $data on (0,1)^d
#         $edg1, $edg2 (d-1)x1 vectors for node labels of edges 1,...,d-1
#          for the tree for residual dependence,
#          $fam for copula family labels for the first tree
#' @param iprfn print gradient and hessian if true
#' @return nllk: approximate nllk; grad: gradient of nllk;
#' hess: hessian matrix of nllk
#' @export
f901f1tproxynllk=function(param,dstruct,iprfn=FALSE){
	udata=dstruct$data
	edg1=dstruct$edg1
	edg2=dstruct$edg2
	d=ncol(udata);
	n=nrow(udata);
	fam=dstruct$fam
	param=param
	npar=length(param)
	lat=dstruct$lat

	out= .Fortran("ft2",as.integer(npar), 
		      as.double(param), 
		      as.integer(d), 
		      as.integer(n),
		      as.double(udata),
		      as.double(lat),
		      as.integer(fam),
		      as.integer(edg1), 
		      as.integer(edg2),
		      nllk=as.double(0.),
		      grad=as.double(rep(0,npar)),
		      hess=as.double(rep(0,npar*npar)),PACKAGE ="FactorModels")
	
  nllk=out$nllk; hess=matrix(out$hess,npar,npar); grad=out$grad;
	
  list(fnval=nllk, grad=grad, hess=hess)
}







#'parameter estimation using proxy method
#'@description parameter estimation using mean version of proxies in 1-factor
#'copula model
#'@param udata N*D data in u-scale
#'@param fam copula families of linking copulas
#'@param start starting values for paramters
#'@param LB upper bound for parameters, length d
#'@param UB lower bound for parameters, length d
#'@param xl  Gaussian-Legendre points
#'@param wl  Gaussian-Legendre weights
#'@param iprint print the gradient and hessian if T
#'@param ifixed paramters are set to be fixed
#'@return mlpx1 estimated parameters from mean-proxy method
#'        mlpx2 estimated parameters from new proxy method
#'        proxyMean mean verison of proxies
#'        proxyNew  new proposed proxies
#'@export
proxy1fctweak<-function(udata,fam,start,LB,UB,xl,wl,iprint,ifixed){
	d=ncol(udata)
	proxyMean=uscore(apply(udata,1,mean)) #mean proxy
	#fam=rep(family,d)
	dstruct=list(data=udata,lat=proxyMean,edg1=c(1:(d-1)),edg2=c(2:d),fam=fam)
	
	out1=pdhessminb(param=start,
			objfn=f901f1tproxynllk,
			dstruct=dstruct,
			ifixed=ifixed,eps=1e-05,
			LB=LB,UB=UB,iprint = F)$parmin[!ifixed]
	
	mlpx1=out1
	
	#new proxy
	lat_update=latUpdateOnefct(th=mlpx1[1:d],udata=udata,nq=25,xl=xl,wl=wl,fam[1])
	proxyNew=uscore(lat_update)
	dstruct=list(data=udata,lat=proxyNew,edg1=c(1:(d-1)),edg2=c(2:d),
							 fam=fam)
	out2=pdhessminb(param=start,objfn=f901f1tproxynllk,dstruct=dstruct,ifixed=ifixed,eps=1e-05,
									 LB=LB,UB=UB,iprint = F)$parmin[!ifixed]
	mlpx2=out2
	
	return(list(mlpx1=mlpx1,mlpx2=mlpx2,proxyMean=proxyMean,proxyNew=proxyNew))
}

