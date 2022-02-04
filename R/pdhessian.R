#' optimize the negative log-likelihood
#' @description Version with ifixed as argument
#' @param param starting point on entry
#' @param objfn function to be minimized with gradient and Hessian
#' @param ifixed = vector of length(param) of True/False, such that
#      ifixed[i]=T iff param[i] is fixed at the given value
#' @param dstruct = list with  data set and other variables used by objfn
#' @param LB = lower bound vector
#' @param UB = upper bound vector
#' @param mxiter = max number of iterations
#' @param eps = tolerance for Newton-Raphson iterations
#' @param iprint = control on amount of printing
#           F for no printing of iterations
#           T for printing x^{(k)} on each iteration
#' @param bdd = bound on difference of 2 consecutive iterations
#' @return fnval = function value at minimum
#' parmin = param for minimum
#' invh = inverse Hessian
#' iconv = 1 if converged, -1 for a boundary point, 0 otherwise
#'
#' @export
#'
pdhessminb=function(param,objfn,ifixed,dstruct,LB,UB,mxiter=30,eps=1.e-6,
										bdd=5,iprint=F)
{ np=length(param)
iter=0; mxdif=1; brk=0;
while(mxdif>eps & iter<mxiter)
{ dstruct$pdf=0;
fn=objfn(param,dstruct,iprint);
fnold=fn$fnval;
# get subset of grad, hess
sgrad=fn$grad[!ifixed]; shess=fn$hess[!ifixed,!ifixed]
evd= eigen(shess,F);
eigval= evd$values;
eigvec= evd$vectors;
ind= (eigval<=0);
#replace negative eigenvalues with a positive number, keeping eigenvectors
if(sum(ind)>0 & iprint) cat("replacing eigenvalues in iteration ", iter,"\n")
eigval[ind]=1e-1;
newhess= eigvec%*%diag(eigval)%*%t(eigvec);
sdiff=try(solve(newhess,sgrad),T);
if(is.vector(sdiff) == FALSE) { brk = 1; break; }
diff=rep(0,np)
diff[!ifixed]=sdiff
mxdif=max(abs(diff))
iter=iter+1
param=param-diff
outb=0;
if(any(param<=LB) | any(param>=UB)) outb=1;
if(is.infinite(fn$fnval) | is.nan(fn$fnval)) outb=1;
if(any(is.infinite(fn$grad)) || any(is.nan(fn$grad))) outb=1;
# modified NR to limit size of step to within bounds
while(outb==1 | mxdif>bdd)
{ mxdif=mxdif/2.;
outb=0;
diff=diff/2; param=param+diff;
if(any(param<=LB) | any(param>=UB)) outb=1;
}
# find a step size to ensure the objective fn decreases
LB.ind= param<=LB +0.008
UB.ind= param>=UB -0.008
ifixed[LB.ind]=T
ifixed[UB.ind]=T
dstruct$pdf=1; # positive definite
fnnew=objfn(param,dstruct,F)$fnval;
brkvr=0;
while(fnold<=fnnew)
{ mxdif=mxdif/2.;
diff=diff/2; param=param+diff;
fnnew= objfn(param,dstruct,F)$fnval;
brkvr=brkvr+1;
if(brkvr>10)
{ print("break: step size is very small");
	fnnew=fnold-1e-5; mxdif=eps-1e-5;
}
}
if(iprint) { cat(iter, fn$fnval, mxdif,"\n") }
}
if(brk==0)
{ iconv=1;
if(iter>=mxiter) { cat("did not converge\n"); iconv=0; }
else { if(max(abs(fn$grad))>eps) iconv=-1 }
# check for positive definite
iposdef=isposdef(shess)
# inverse Hessian
invh=solve(shess)
}
if(brk==1)
{ iconv=2; #algorithm failed
param=rep(1e5,np);
fnnew=1e5;
invh=NA;
isposdef=FALSE;
}
list(parmin=param,fnval=fnnew,invh=invh,iconv=iconv,iposdef=iposdef)
}

