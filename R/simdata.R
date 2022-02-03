#' simulate data from 1-factor copula model
#' @description simulate data from 1-factor model
#' @param n sample size
#' @param parobj1 parameter vector of dimension d
#' @param qcond1 function for copula conditional inverse cdf
#' @param copname1 copula family in the copula model
#' @param ivect flag that is T if qcond1 have vectorized forms
#' @return yy N*D simulate data in the u-scale
#' @return vv simulated latent variables
#' @export

sim1fact2<-function(n,parobj1, qcond1, copname1,ivect)
{
	vv = runif(n)
	if (is.vector(parobj1))
		d = length(parobj1)
	else d = nrow(parobj1)
	copname = tolower(copname1)
	yy = matrix(0, n, d)
	if (copname == "frank" | copname == "frk" | copname == "mtcj" |
			copname == "mtcjr" | copname == "fgm" | ivect) {
		for (j in 1:d) {
			qq = runif(n)
			if (is.vector(parobj1))
				th = parobj1[j]
			else th = parobj1[j, ]
			yy[, j] = qcond1(qq, vv, th)
		}
	}
	else {
		for (j in 1:d) {
			if (is.vector(parobj1))
				th = parobj1[j]
			else th = parobj1[j, ]
			for (i in 1:n) {
				v1 = vv[i]
				qq = runif(1)
				yy[i, j] = qcond1(qq, v1, th)
			}
		}
	}
	return(list("dat"=yy,"vlat"=vv))
}










#' simulate data from 1-factor copula model with different
#' copula familes
#' @description simulate data from 1-factor copula model
#' @param n sample size
#' @param d dimension
#' @param param parameter vector of dimension 2*d (par1,par2) for
#' d linking copulas, one-parameter set par2=0
#' @param fam family vector for d linking copulas
#' same index with the vinecopula package, 4:gumbel
#' 5:frank copula;  2: student-t copula
#' @return yy N*D simulate data in the u-scale
#' @return vv simulated latent variables
#' @export

sim1fact3<-function(n,d,param,fam){
	vv = runif(n)
	yy = matrix(0, n, d)

	for (j in 1:d) {
		jj=c(2*j-1,2*j)
		th = param[jj]
		for (i in 1:n) {
			v1 = vv[i]
			qq = runif(1)
			bicop=BiCop(family=fam[j],par=th[1],par2=th[2])
			yy[i,j]=BiCopHinv(u1=qq,u2=v1,bicop)$hinv2
		}
	}

	return(list("dat"=yy,"vlat"=vv))
}







#'simulate data from bi-factor copula model
#'@description simulate data from bi-factor copula model
#'@param nn sample size
#'@param grsize group size in the model
#'@param cop  copula familes in the model
#'@param param paramters on the linking copulas
#'@return zdata simulate data in the u-scale
#'@return v0 simulated glocal latent variables
#'@return vg simulate local latent variables
#'@export
#'
simbifact2<-function(nn, grsize, cop = 5, param)
{
	ivect = F
	d = sum(grsize)
	mgrp = length(grsize)
	th1 = param[1:d]
	th2 = param[(d + 1):(2 * d)]
	if (cop == 2) {
		df0 = param[2 * d + 1]
	}
	if (cop == 3) {
		qcond = qcondgum
		ivect = T
	}
	if (cop == 5) {
		qcond = qcondfrk
		ivect = T
	}
	if (cop == 9) {
		qcond1 = qcondbb1
		qcond2 = qcondfrk
		th3 = param[(2 * d + 1):(3 * d)]
	}
	if (cop == 1 || cop == 2)
		lm2 = th2 * sqrt(1 - th1^2)
	zdata = matrix(0, nrow = nn, ncol = d)
	if (cop == 1 || cop == 2) {
		z0 = rnorm(nn)
		z = matrix(rnorm(nn * mgrp), ncol = mgrp)
	}
	else {
		z0 = runif(nn)
		z = matrix(runif(nn * mgrp), ncol = mgrp)
	}
	ind = 0
	for (jg in 1:mgrp) {
		ind1 = ind + 1
		ind2 = ind + grsize[jg]
		ind = ind + grsize[jg]
		for (ij in ind1:ind2) {
			if (cop == 1 || cop == 2) {
				zdata[, ij] = th1[ij] * z0 + z[, jg] * lm2[ij] +
					sqrt(1 - th1[ij]^2 - lm2[ij]^2) * rnorm(nn)
			}
			else if (ivect) {
				q1 = qcond(runif(nn), z[, jg], th2[ij])
				zdata[, ij] = qcond(q1, z0, th1[ij])
			}
			else if (cop == 9) {
				for (i in 1:nn) {
					q1 = qcond2(runif(1), z[i, jg], th3[ij])
					zdata[i, ij] = qcond1(q1, z0[i], c(th1[ij],
																						 th2[ij]))
				}
			}
			else {
				for (i in 1:nn) {
					q1 = qcond(runif(1), z[i, jg], th2[ij])
					zdata[i, ij] = qcond(q1, z0[i], th1[ij])
				}
			}
		}
	}
	if (cop == 2) {
		for (i in 1:nn) {
			zdata[i, ] = zdata[i, ]/sqrt(rchisq(1, df = df0)/df0)
		}
	}
	return(list(zdata=zdata,v0=z0,vg=z))
}







#'simulate data from nested oblique factor model
#'@description simulate data from nested oblique factor model
#'@param nn sample size
#'@param grsize group size in the model
#'@param cop1 linking copula familes in T1
#'@param cop2 linking copula families in each group
#'@param param paramters on the linking copulas
#'@return zdata data in the u-scale
#'@return lat0 global latent variables
#'@return lat local latent variables
#'@export
simnestfact2=function(nn,grsize,cop1,cop2,param){
	d=sum(grsize);
	mgrp=length(grsize);
	if(cop1==2) { df0=param[mgrp+d+1]; }
	if(cop1==1 || cop1==2) { z0=rnorm(nn); } else { z0=runif(nn); }
	z=matrix(0,nrow=nn,ncol=mgrp);
	zdata=matrix(0,nrow=nn,ncol=d)
	if(cop1==1 || cop1==2)
	{ for (jg in 1:mgrp) { z[,jg]= z0*param[jg]+ sqrt(1-param[jg]^2)*rnorm(nn); }}
	else if(cop1==5)
	{ for (jg in 1:mgrp) { z[,jg]= qcondfrk(runif(nn),z0,param[jg]); } }
	else if(cop1==3 || cop1==10)
	{ for (jg in 1:mgrp)
	{ for (i in 1:nn) { z[i,jg]= qcondgum(runif(1),z0[i],param[jg]); } }
	}
	ind = 0;
	for (jg in 1:mgrp)
	{ ind1=ind+1;
	ind2=ind+grsize[jg];
	ind =ind+grsize[jg];
	for(ij in ind1:ind2)
	{ ijm=ij+mgrp
	if(cop2==1 || cop2==2)
	{ zdata[,ij]=z[,jg]*param[ijm]+sqrt(1-param[ijm]^2)*rnorm(nn); }
	else if(cop2==5) { zdata[,ij]=qcondfrk(runif(nn),z[,jg],param[ijm]); }
	else if(cop2==3)
	{ for(i in 1:nn) { zdata[i,ij]=qcondgum(runif(1),z[i,jg],param[ijm]); } }
	else if(cop2==10)
	{ for(i in 1:nn)
	{ zdata[i,ij]=qcondbb1(runif(1),z[i,jg],param[c(ijm,ijm+d)]); }
	}
	else print("this copula family is not implemented");
	}
	}
	if(cop2==2)
	{ for(i in 1:nn)
	{ zdata[i,]=zdata[i,]/sqrt(rchisq(1,df=df0)/df0);
	}
	}
	return(list(zdata=zdata,lat0=z0,lat=z))
}




# sim1fact4<-function(n,d,param,fam){
# 	vv = runif(n)
# 	yy = matrix(0, n, d)
#
# 	for (j in 1:d) {
# 		jj=c(2*j-1,2*j)
# 		 th = param[jj]
# 			qq = runif(n)
# 			if(fam[j]==4){
# 				yy[,j]=qcondgum(qq, vv, th[1])
# 			}else if(fam[j]==5){
# 					yy[,j]=qcondfrk(qq,vv,th[1])
# 			}else if(fam[j]==2){
# 					yy[,j]=qcondbvtcop(qq,vv,cpar=th[1],df=5)
# 				}
# 		}
#
#
# 	return(list("dat"=yy,"vlat"=vv))
# }
#

#' #'simulate data from 1-factor with weak res dependence
#' #'@description simulate data from 1-factor with weak residual
#' #'depdence, the weak dependence is modeled by 1-truncated
#' #'D-vine
#' #'@param n sample size
#' #'@param d dimension
#' #'@param thlow1 lower bound for param in tree1
#' #'@param thupp1 upper bound for param in tree1
#' #'@param thlow2 lower bound for param in res tree
#' #'@param thupp2 upper bound for param in res tree
#' #'@param fam1 copula family in T1
#' #'@param fam2 copula family in residual tree
#' #'@return simdata simulate data in u-scale
#' #'@return thv paramters in the model
#' #'@export
#' #'
#' simWeakRes<-function(n,d,thlow1=1.67,thupp1=5.0,
#' 										 thlow2=0.95,thupp2=2.95,
#' 										 fam1=4,fam2=5){
#' 	d1=d+1
#' 	Matrix <-matrix(NA,d1,d1)
#' 	Matrix[1,]<-rep(d1,d1)
#' 	m=Dvinearray(d1-1)
#' 	flag=lower.tri(m,F)
#' 	m[flag]<-0
#' 	Matrix[2:d1,2:d1]<-m
#' 	Matrix[is.na(Matrix)]<-0
#' 	family<-c(0,rep(fam1,d1-1),0,0,rep(fam2,d1-2),rep(0,(d1-2)*d1))
#' 	family <- matrix(family, d1, d1,T)
#' 	thv1=runif(d1-1,thlow1,thupp1)
#' 	thv2=runif(d1-2,thlow2,thupp2)
#' 	thv=c(thv1,thv2)
#' 	par<-c(0,thv1,0,0,thv2,rep(0,(d1-2)*d1))
#' 	par<-matrix(par,d1,d1,T)
#' 	par2 <- matrix(0, d1, d1)
#' 	RVM <- RVineMatrix(Matrix = Matrix, family = family,
#' 										 par = par, par2 = par2,
#' 										 names = paste0("V",c(1:d1)))
#' 	simdata <- RVineSim(n, RVM)
#' 	return(list(simdata=simdata,thv=thv))
#' }
#'
#'
#'
#'
#'
#'

