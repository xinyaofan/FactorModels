! Code written by Pavel Krupskii and revised by HJ
! The subroutine lgum1derivs is in file gum12fact.f90
! The subroutine lfrkderivs is in file frk12fact-lpdf.f90

! nested-factor model with Frank linking copulas for group to common latent
!         and Gumbel linking copulas from observed to group latent.
! inputs
!   npar = mgrp+dvar = number of parameters
!   th = vector of parameters; 
!        parameters for links of group latent to common latent go first
!   mgrp = number of groups 
!   n = sample size
!   dvar = sum(grsize) = dimension of data
!   grsize = vector with group sizes
!   udata = nxd matrix of uniform scores
!   lat=n * mgrp estiamted proxies 
!   nq = number of quadrature points
!   wl = vector of quadrature weights
!   xl = vector of quadrature nodes
! outputs 
!   nllk = negative log-likelihood, 
!   grad = gradient of nllk, 
!   hess = hessian of nllk
subroutine strfrkgum1(npar,th,mgrp,n,dvar,grsize,udata,lat,nq,wl,xl,nllk,grad,hess)
  implicit none
  integer npar,mgrp,dvar,n,nq,ip,jp, ind
  double precision th(npar),udata(n,dvar),wl(nq),xl(nq),uvec(dvar),lat(n,mgrp)
  double precision lpdf(npar),der1(npar),der2(npar),llk
  double precision nllk,liki,lk,grad(npar),hess(npar,npar)
  integer i,iq,jg,jg2,mj,mj2,nj,nj2,ind1,ind2,grsize(mgrp)
  double precision intj(mgrp), intjj(npar), grdj(npar), der2j(npar), int0, grd0(npar), grdi(npar)
  double precision hss0(npar,npar), hssi(npar,npar), hessi(npar,npar), hssaux(npar,npar)

  ! npar = dvar+mgrp; dvar = sum(grsize)
  nllk=0.d0; grad = 0.d0; hess=0.d0 
  do i =1,n 
    uvec = udata(i,:)  
    liki = 0.d0; grdi = 0.d0; hssi = 0.d0; 
    do iq =1,nq
      int0 = 1.d0; grd0 = 1.d0; hssaux = 0.d0 
      ind = 0; intj = 0.d0;  grdj = 0.d0; der2j = 0.d0 
      do jg =1,mgrp   !jth group                         
        !jg2=jg+dvar
        jg2=jg
        ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg);      
        !do iq2 =1,nq   
          lk = 1.d0      
          do mj = ind1,ind2 ! index for within subgroup
            nj=mj+mgrp
            call lgum1derivs(uvec(mj),lat(i,jg),th(nj),lpdf(nj),der1(nj),der2(nj))   !c_{ij,V_j}   
            lk = lk*exp(lpdf(nj))
          end do
          call lfrk1derivs(lat(i,jg),xl(iq),th(jg2),lpdf(jg2),der1(jg2),der2(jg2))   !c_{V_j,V_0}  
          llk = exp(lpdf(jg2))*lk  !llk: values of the jth integrand
          intj(jg) = llk !intj value of the jth integral
          !intjj(ind1:ind2) = intj(jg) !intjj: repeat intj(j) grsize[j] times
          intjj(ind1+mgrp:ind2+mgrp) = intj(jg) !intjj: repeat intj(j) grsize[j] times
          do mj = ind1,ind2
            nj=mj+mgrp
            grdj(nj) = llk*der1(nj) !grdj: der. of the jth int. wrt th(mj)
            der2j(nj)= llk*(der2(nj)+der1(nj)*der1(nj)) !der2j: second order der. wrt th(mj)
            hssaux(jg2,nj) = llk*der1(nj)*der1(jg2)
            hssaux(nj,jg2) = hssaux(jg2,nj)
            if (mj < ind2) then
              do mj2 = (mj+1),ind2  !2nd order der. of the jth int. wrt th(mj), th(mj2)
                nj2=mj2+mgrp
                hssaux(nj,nj2) = llk*der1(nj)*der1(nj2)
                hssaux(nj2,nj) = hssaux(nj,nj2)
              end do 
            endif         
          end do
          grdj(jg2) = llk*der1(jg2) 
          der2j(jg2)= llk*(der2(jg2)+der1(jg2)*der1(jg2))
        end do
      !end do
      !intjj(dvar+1:dvar+mgrp)=intj
      intjj(1:mgrp)=intj
      int0 = product(intj)   !product of j inner integrals  
      grd0 = int0*(grdj/intjj) !gradient of the product
      do ip=2,npar !hessian of the product
        do jp=1,ip-1
          hss0(ip,jp) = int0*grdj(ip)*grdj(jp)/intjj(ip)/intjj(jp)
          hss0(jp,ip) = hss0(ip,jp)
        end do
      end do
      !ind = 0
      ind = mgrp
      do jg=1,mgrp   !2nd order derivatives within groups 
        !jg2=jg+dvar
        jg2=jg
        ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg)
        hss0(ind1:ind2,ind1:ind2)=int0*hssaux(ind1:ind2,ind1:ind2)/intj(jg)
        hss0(ind1:ind2,jg2)=int0*hssaux(ind1:ind2,jg2)/intj(jg)
        hss0(jg2,ind1:ind2)=int0*hssaux(jg2,ind1:ind2)/intj(jg)
      end do
      do ip=1,npar
        hss0(ip,ip) = int0*der2j(ip)/intjj(ip)
      end do
      liki = liki + wl(iq)*int0   !integrating over V0
      grdi = grdi + wl(iq)*grd0
      hssi = hssi + wl(iq)*hss0
    end do
    nllk = nllk - log(liki)       !updating loglikelihood
    grad = grad - grdi/liki       !updating gradient
    do ip=1,npar                  !updating hessian
      do jp=1,ip
        hessi(ip,jp)=hssi(ip,jp)/liki-grdi(ip)*grdi(jp)/liki/liki 
        hessi(jp,ip)=hessi(ip,jp)
      end do
    end do 
    hess = hess - hessi 
  end do
  return
  end

! Function with Frank lpdf, lder1, lder2 (partial wrt cpar and 2nd order)
!  for factor 1.
! inputs
!   u1,u2 = values in (0,1)
!   cpar = scalar in (-oo,oo) but could be underflow problems for cpar>35
! outputs 
!   lpdf = log pdf, 
!   lder1 = \p lpdf/\p cpar
!   lder2 = \p^2 lpdf/\p cpar^2
subroutine lfrk1derivs(u1,u2,cpar,lpdf,lder1,lder2)
  implicit none
  double precision u1,u2,cpar,lpdf,lder1,lder2
  double precision den,den1,den2,t0,t1,t2
  t0 = exp(-cpar);
  t1 = exp(-cpar*u1);
  t2 = exp(-cpar*u2);
  den = t1+t2-t0-t1*t2;
  den1 = -u1*t1-u2*t2+t0+(u1+u2)*t1*t2;
  den2 = u1*u1*t1+u2*u2*t2-t0-(u1+u2)*(u1+u2)*t1*t2;
  lpdf = log(abs(cpar))+log(abs(1.d0-t0))-cpar*(u1+u2)-2.d0*log(abs(den));
  lder1 = 1.d0/cpar+t0/(1.d0-t0)-(u1+u2)-2.d0*den1/den;
  lder2 = -1.d0/(cpar*cpar)-t0/((1.d0-t0)*(1.d0-t0))-2.d0*den2/den+2.d0*den1*den1/(den*den); 
  return
  end
!--------------------------------------------------------------------
! Remaining code below is translated from R code written by Pavel Krupskii

! Function with Gumbel lpdf, lder1, lder2 (partial wrt cpar and 2nd order)
!  for factor 1.
! For derivatives, lpdf is function of cpar and m=(x^cpar+y^cpar)^(1/cpar)
! inputs
!   u1,u2 = values in (0,1)
!   cpar = scalar in (1,oo)
! outputs 
!   lpdf = log pdf, 
!   lder1 = \p lpdf/\p cpar
!   lder2 = \p^2 lpdf/\p cpar^2
subroutine lgum1derivs(u1,u2,cpar,lpdf,lder1,lder2)
  implicit none
  double precision u1,u2,cpar,lpdf,lder1,lder2
  double precision x,y,tx,ty,xd,yd,s,m,logs,logm,msq,dlsq,dlcu
  double precision sder1,sder2,mder1,mder2,den,den2
  x = -log(u1); y = -log(u2);
  tx = log(x); ty = log(y); 
  xd = x**cpar; yd = y**cpar; s = xd+yd; m = s**(1.d0/cpar);
  logs=log(s);
  dlsq=cpar*cpar; dlcu=dlsq*cpar
  sder1 = xd*tx+yd*ty;
  sder2 = xd*tx*tx+yd*ty*ty;
  mder1 = m*sder1/(s*cpar)-m*logs/dlsq;
  mder2 = -mder1*logs/dlsq-2.d0*m*sder1/(s*dlsq)+2.d0*m*logs/dlcu;
  mder2 = mder2+sder2*m/(s*cpar)+(mder1/s-m*sder1/(s*s))*sder1/cpar;
  den = m+cpar-1.d0; den2=den*den
  logm=log(m); msq=m*m
  lpdf = -m+log(den)+(1.d0-2.d0*cpar)*logm+(cpar-1.d0)*(tx+ty)+x+y;
  lder1 = -mder1+(mder1+1)/den-2*logm+(1-2*cpar)*mder1/m+tx+ty;
  lder2 = -mder2+mder2/den-(mder1+1.d0)**2/den2-4.d0*mder1/m+(1.d0-2.d0*cpar)*(mder2/m-mder1*mder1/msq);
  return
  end