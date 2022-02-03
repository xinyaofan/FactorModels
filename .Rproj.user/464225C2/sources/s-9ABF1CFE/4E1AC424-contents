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
  th=matrix(th,2,d)
	out=.Fortran("latupdate2",as.double(th),as.integer(n),as.integer(d),
							 as.double(udata),as.integer(nq),as.double(xl),as.double(wl),
							 as.integer(family),
							 lat=rep(0,n),PACKAGE = "FactorModels")
	return(out)
}










#'parameter estimation using proxy method in 1-factor copula model
#'with all copula families in the same family
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
	lat_update=latUpdateOnefct(th=mlpx1,udata=udata,nq=25,xl=xl,wl=wl,
														 family=family)
	proxyNew=uscore(lat_update)
	dstruct=list(udata=udata,lat=proxyNew,fam=fam)
	out2=pdhessminb(start,f901fproxynllk,ifixed = rep(F,d),dstruct=dstruct,
									LB=LB,UB=UB,mxiter=30,eps=1.e-4,iprint=F)
	mlpx2=out2$parmin
	return(list(mlpx1=mlpx1,mlpx2=mlpx2,proxyMean=proxyMean,proxyNew=proxyNew))
}








#'return the optimization upper and lower bound for
#'the parameters in 1-factor copula model
#'@description upper and lower bound for parameters in 1-factor copula
#'model
#'@param familyvec a vector of length d for copula families
#'@return UB1 upper bound for par1
#'UB2 upper bound for par2
#'LB1 lower bound for par1
#'LB2 lower bound for par2
#'@export
#'
return_bound<-function(familyvec){
	n=length(familyvec)
	UB1=rep(NA,n)
	LB1=rep(NA,n)
	UB2=rep(NA,n)
	LB2=rep(NA,n)
	for(i in 1:n){
		tem=familyvec[i]
		if(tem==2){
			LB1[i]=-0.99
			UB1[i]=0.99
			LB2[i]=2.01
			UB2[i]=5000
		} else if (tem==17|tem==7) {
			LB1[i]=0
			UB1[i]=7
			LB2[i]=1
			UB2[i]=7
		} else if (tem==40){
			LB1[i]=-8
			UB1[i]=-1
			LB2[i]=-1
			UB2[i]=-1e-4
		} else if (tem==1) {
			LB1[i]=-0.99
			UB1[i]=0.99
			LB2[i]=0
			UB2[i]=0
		} else if (tem==5){
			LB1[i]=-99.99
			UB1[i]=99.99
			LB2[i]=0
			UB2[i]=0
		} else if(tem==20){
			UB1[i]=8
			UB2[i]=1
			LB1[i]=1
			LB2[i]=0.00001
		} else if(tem==10){
			UB1[i]=8
			LB1[i]=1
			UB2[i]=1
			LB2[i]=0.00001
		} else if(tem==14){
			UB1[i]=100
			LB1[i]=1
			UB2[i]=0
			LB2[i]=0
		}else if(tem==4){
			UB1[i]=50
			LB1[i]=1
			UB2[i]=0
			LB2[i]=0
		}
	}
	return(list(UB1=UB1,UB2=UB2,LB1=LB1,LB2=LB2))
}












#'find a good starting point for 1-factor copula model
#'@description find a good starting point for optimization in
#'1-factor copula model with diffferent copula families
#'@param fam a vector of length d indicating the copula families
#'@param tau a vector of length d indicating the kendall's tau
#'@param df the fixed averaged degree of freedom for student-t copula
#'@return st a vector of length 2*d for starting point
#'@export
#'

find_start<-function(fam,tau,df){
	d=length(fam)

	st=rep(0,2*d)
	for(i in 1:d){
		if(fam[i]==4){
			st[2*i-1]=gum.tau2cpar(tau=tau[i])
		}else if(fam[i]==5){
			st[2*i-1]=BiCopTau2Par(family=5,tau=tau[i])
		}else if(fam[i]==2){
			st[2*i-1]=BiCopTau2Par(family=2,tau=tau[i])
			st[2*i]=df
		}else if(fam[i]==17){
			st[c(2*i-1,2*i)]=bb1.tau2eqlm(tau[i])[c(1,2)]
		}else if(fam[i]==7){
			st[c(2*i-1,2*i)]=bb1.tau2eqlm(tau[i])[c(1,2)]
		}else if(fam[i]==10|fam[i]==20)
		{
			st[c(2*i-1,2*i)]=c(4.5,0.5)
		}
		#print(i)
	}
	return(st)
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
#'@param ifixed if ifixed=T, the parameters at j-th position will be
#'fixed
#'@param iprint print the gradient and hessian if T
#'@return mlpx1 estimated parameters from mean-proxy method
#'        mlpx2 estimated parameters from new proxy method
#'        proxyMean mean verison of proxies
#'        proxyNew  new proposed proxies
#'@export
proxy1fct2<-function(udata,fam,start,LB,UB,ifixed,xl,wl,iprint){
	d=ncol(udata)
	proxyMean=uscore(apply(udata,1,mean)) #mean proxy

	dstruct=list(udata=udata,lat=proxyMean,fam=fam)
	out1=pdhessminb(start,f901fproxynllk2,ifixed =ifixed,
									dstruct=dstruct,
									LB=LB,UB=UB,mxiter=30,eps=1.e-4,iprint=T)
	mlpx1=out1$parmin
	#new proxy
	# lat_update=latUpdateOnefct2(th=mlpx1,udata=udata,
	# 														nq=25,xl=xl,wl=wl,family=fam)
	# proxyNew=uscore(lat_update)
	# dstruct=list(udata=udata,lat=proxyNew,fam=fam)
	# out2=pdhessminb(mlpx1,f901fproxynllk2,ifixed = ifixed,dstruct=dstruct,
	# 								LB=LB,UB=UB,mxiter=30,eps=1.e-4,iprint=F)
	# mlpx2=out2$parmin
	return(list(mlpx1=mlpx1,proxyMean=proxyMean))
	#return(list(mlpx1=mlpx1,mlpx2=mlpx2,proxyMean=proxyMean,proxyNew=proxyNew))
}
