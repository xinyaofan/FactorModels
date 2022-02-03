! program test
!   implicit none
!   double precision u1,u2,cparv(2)
!   integer icop
!   !double precision m,mder1th,mder1dl,mder2th,mder2dl,mderthd,mderu(4),mderv(4),mder2uv
!   !double precision m,mder1th,mder1dl,mder2th,mder2dl,mderthd,mderu(4),mderv(4),mder2uv
!   double precision lpdf,lder11(2),lder22(3),ldermixvv(2),ldermixuu(2)
!   double precision lde2u,lder1u,lder1v,lder2u,lder2v,lder2uv
!    
!    	read*,u1
!   	read*,u2
!   	read*,cparv
!    	read*,icop
!    	
!    	print*,u1
!   	print*,u2
!    	print*,cparv
!    	print*,icop
!    	
!    
!    	call lcop2derivt(u1,u2,cparv,icop,lpdf,lder11,lder22,lder1u,lder2u,&
 !ldermixuu,lder1v,lder2v,lder2uv,ldermixvv)
!    	print*,"lpdf"
!    	print*,lpdf
!    	print*, "der for cparv"
!     print*,lder11
!     print*,"second der for cparv"
!     print*,lder22
!     print*,"der for u"
!     print*,lder1u
!     print*,"second der for u"
!     print*,lder2u
!     print*,"mix der for u and cpar"
!     print*,ldermixuu
!     print*,"der for v"
!     print*,lder1v
!     print*,"second der for v"
!     print*,lder2v
!     print*,"mix der for uv"
!     print*,lder2uv
!     print*,"mix der for v and cpar"
!     print*,ldermixvv
! 
! stop 
! end

!this function computes the derivatives for copula log pdf=c(u1,u2,cparv)
!log density =lpdf
!lder11= \p logpdf/\cparv 2d-vector
!lder22=\p^2 logpdf/\p cpar1^2,\p^2 logpdf/\p cpar1*cpar2,\p^2 logpdf/\p cpar2^2 3d-vector
!lder1u=\p lpdf/\p u1
!lder2u=\p^2 lpdf/\p u1^2
!ldermixuu=\p^2 lpdf\p^2 u1*cpar1, \p^2 lpdf\p^2 u1*cpar2 2d-vector
!lder1v=\p lpdf/\p u2
!lder2v=\p^2 lpdf/\p u2^2
!ldermixvv=\p^2 lpdf\p^2 u2*cpar1, \p^2 lpdf\p^2 u2*cpar2 2d-vector
!lder2uv=\p^2 lpdf\p^2 u1*u2
subroutine lcop2derivt(u1,u2,cparv,icop,lpdf,lder11,lder22,lder1u,lder2u,&
ldermixuu,lder1v,lder2v,lder2uv,ldermixvv)
  implicit none
  integer icop
  double precision u1,u2,t1,t2,ut,vt,cparv(2),lder1u,lder2u,lder1v,lder2v,lder2uv
  double precision ltder1u,ltder2u,ltdermix
  double precision lder11(2),lder22(3),lder1,lder2,ldermixv,lderuv,ldermixu
  double precision lderu(4),lderv(4),ldermixuu(2),ldermixvv(2)
  double precision rho,nu,cpar,lpdf
  double precision qt
  cpar=cparv(1)
  select case(icop)
    case (1) ! bivariate normal/Gaussian !check
    call lgau2derivt(u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,ldermixu,lder1v,&
    lder2v,lder2uv,ldermixv)
    lder11(1)=lder1;
    lder11(2)=0.0d0;
    lder22(1)=lder2;
    lder22(2)=0.0d0;
    lder22(3)=0.0d0;
    ldermixuu(1)=ldermixu;
    ldermixuu(2)=0.0d0;
    ldermixvv(1)=ldermixv;
    ldermixvv(2)=0.0d0;
    case (2) !student t !check
    rho=cparv(1);nu=cparv(2);
    t1=qt(u1,nu);t2=qt(u2,nu);
    call lt2derivs(t1,t2,rho,nu,lpdf,lder1,lder2,ltder1u,ltder2u,ltdermix)
    lder11(1)=lder1; 
    lder11(2)=0.0d0;
    lder22(1)=lder2; 
    lder22(2)=0.0d0;
    lder22(3)=0.0d0;
    lder1u=ltder1u; 
    lder2u=ltder2u; 
    ldermixuu(1)=ltdermix; 
    ldermixuu(2)=0.0d0;
    ldermixvv(1)=0.0d0;
    ldermixvv(2)=0.0d0;
    lder2uv=0.0d0;
    lder1v=0.0d0;lder2v=0.0d0;
    case (4) ! gumbel !check
    call lgum2derivt(u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,ldermixu,&
    lder1v,lder2v,lder2uv,ldermixv)
    lder11(1)=lder1; lder11(2)=0.0d0;
    lder22(1)=lder2; lder22(2)=0.0d0;lder22(3)=0.0d0;
    ldermixuu(1)=ldermixu; 
    ldermixuu(2)=0.0d0;
    ldermixvv(1)=ldermixv;
    ldermixvv(2)=0.0d0;
    case (5) ! frank !check
    call lfrk2derivt(u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,ldermixu, &
    lder1v,lder2v,lder2uv,ldermixv)
    lder11(1)=lder1; 
    lder11(2)=0.0d0;
    lder22(1)=lder2; 
    lder22(2)=0.0d0;
    lder22(3)=0.0d0;
    ldermixuu(1)=ldermixu;
    ldermixuu(2)=0.0d0;
    ldermixvv(1)=ldermixv; 
    ldermixvv(2)=0.0d0;
    case (7) !bb1 copula
    call lbb1derivs(u1,u2,cparv,lpdf,lder11,lder22,lderu,lderv,lderuv)
    lder11=lder11;
    lder22=lder22;
    lder1u=lderu(1)
    lder2u=lderu(2)
    ldermixuu(1)=lderu(3); 
    ldermixuu(2)=lderu(4);
    lder1v=lderv(1); 
    lder2v=lderv(2);
    lder2uv=lderuv;
    ldermixvv(1)=lderv(3);
    ldermixvv(2)=lderv(4);
    case (10) !bb8 copula
    call lbb8derivs(u1,u2,cparv,lpdf,lder11,lder22,lderu,lderv,lderuv)
    lder11=lder11;
    lder22=lder22;
    lder1u=lderu(1)
    lder2u=lderu(2)
    ldermixuu(1)=lderu(3); 
    ldermixuu(2)=lderu(4);
    lder1v=lderv(1); 
    lder2v=lderv(2);
    lder2uv=lderuv;
    ldermixvv(1)=lderv(3);
    ldermixvv(2)=lderv(4);
    case(14) !survival gumbel !check
    ut=1.0d0-u1
    vt=1.0d0-u2
    call lgum2derivt(ut,vt,cpar,lpdf,lder1,lder2,lder1u,lder2u,ldermixu,&
    lder1v,lder2v,lder2uv,ldermixv)
    lder11(1)=lder1; 
    lder11(2)=0.0d0;
    lder22(1)=lder2; 
    lder22(2)=0.0d0;
    lder22(3)=0.0d0;
    !lder1u=-lderu(1);
    lder1u=-lder1u;
    !lder2u=lderu(2);
    !lder1v=-lderv(1);
    lder1v=-lder1v;
    !lder2v=lderv(2);
    lder2uv=lder2uv;
    ldermixuu(1)=-ldermixu;
    ldermixuu(2)=0.0d0;
    ldermixvv(1)=-ldermixv;
    ldermixvv(2)=0.0d0;
    case (17) !survival bb1 copula
    ut=1.0d0-u1
    vt=1.0d0-u2
    call lbb1derivs(ut,vt,cparv,lpdf,lder11,lder22,lderu,lderv,lderuv)
    lder11=lder11;
    lder22=lder22;
    lder1u=-lderu(1)
    lder2u=lderu(2)
    ldermixuu(1)=-lderu(3);
    ldermixuu(2)=-lderu(4);
    lder1v=-lderv(1);
    lder2v=lderv(2);
    lder2uv=lderuv;
    ldermixvv(1)=-lderv(3);
    ldermixvv(2)=-lderv(4);
    case (20) !survival bb8
    ut=1.0d0-u1
    vt=1.0d0-u2
    call lbb8derivs(ut,vt,cparv,lpdf,lder11,lder22,lderu,lderv,lderuv)
    lder11=lder11;
    lder22=lder22;
    lder1u=-lderu(1)
    lder2u=lderu(2)
    ldermixuu(1)=-lderu(3);
    ldermixuu(2)=-lderu(4);
    lder1v=-lderv(1);
    lder2v=lderv(2);
    lder2uv=lderuv;
    ldermixvv(1)=-lderv(3);
    ldermixvv(2)=-lderv(4);
  end select
  return
end

  
! 1st and 2nd order derivatives of m = s^(1/delta) wrt theta and delta
! inputs
!   u1,u2 = values in (0,1)
!   theta >0
!   delta >1
! outputs
!   m = s^(1/delta); s =  [u1^(-theta)-1]^delta + [u2^(-theta)-1]^delta
!   mder1th = \p m/\p theta
!   mder1dl = \p m/\p delta
!   mder2th = \p^2 m/\p theta^2
!   mder2dl = \p^2 m/\p delta^2
!   mderthd = \p^2 m/\p theta \p delta
subroutine mderivs(u1,u2,cparv,m,mder1th,mder1dl,mder2th,mder2dl,mderthd,&
mderu,mderv,mder2uv)
  implicit none
  double precision u1,u2,theta,delta,m,mder1th,mder1dl,mder2th,mder2dl,mderthd
  double precision t1,t2,tu1,tu2,ttu1,ttu2,td01,td02,td11,td12,td21,td22
  double precision s,sder1th,sder1dl,sder2th,sder2dl,sderthd,m1,ts,m1der1dl 
  double precision t1a,t2a,dlsq,dlcu
  double precision xder1u,yder1v,xder2u,mder1u,mder1v,mder2u,mder2uv,mder2v,del1,th1
  double precision yder2v,xder2uth,yder2vth,xder1th,xder2udl,yder2vdl,yder1th
  double precision mder2vdl,mder2vth,mder2udl,mder2uth
  double precision mderu(4),mderv(4),cparv(2)
  theta=cparv(1);delta=cparv(2)
  del1=1.0d0/delta-1.0d0
  th1=-theta-1.0d0
  
  t1 = u1**(-theta); t2 = u2**(-theta)
  t1a=t1-1.d0; t2a=t2-1.d0
  tu1 = -log(u1); tu2 = -log(u2)
  ttu1 = log(t1a); ttu2 = log(t2a)
  td01 = (t1a)**delta; td02 = (t2a)**delta
  td11=td01/t1a; td12=td02/t2a;
  td21=td11/t1a; td22=td12/t2a;
 
  s = td01+td02
  sder1th = delta*(td11*t1*tu1+td12*t2*tu2)
  sder1dl = td01*ttu1+td02*ttu2
  sder2th = delta*(delta-1.d0)*(td21*t1*t1*tu1*tu1+td22*t2*t2*tu2*tu2)
  sder2th = sder2th+delta*(td11*t1*tu1*tu1+td12*t2*tu2*tu2)
  sder2dl = td01*ttu1*ttu1+td02*ttu2*ttu2
  sderthd = sder1th/delta+delta*(td11*ttu1*tu1*t1+td12*ttu2*tu2*t2)

  m = s**(1.d0/delta); m1 = m/s
  ts = log(s)
  dlsq=delta*delta; dlcu=delta*dlsq
  mder1th = m1*sder1th/delta
  mder1dl = m1*sder1dl/delta - m*ts/dlsq
  m1der1dl = mder1dl/s - m*sder1dl/s**2
  mder2th = (1.d0-delta)*m1*sder1th**2/(dlsq*s)+m1*sder2th/delta
  mder2dl = 2.d0*m*ts/dlcu-mder1dl*ts/dlsq-2.d0*m1*sder1dl/dlsq
  mder2dl = mder2dl+sder2dl*m1/delta+sder1dl*m1der1dl/delta
  mderthd = -m1*sder1th/dlsq+sder1th*m1der1dl/delta+m1*sderthd/delta
  
  xder1u=-delta*theta*td11*u1**(th1) !correct
  yder1v=-delta*theta*td12*u2**(th1) !correct
  xder2u=-delta*theta*(-theta*(delta-1.0d0)*td21*(u1**(2*th1))+&
  th1*td11*u1**(th1-1.0d0))!correct
  yder2v=-delta*theta*(-theta*(delta-1.0d0)*td22*(u2**(2*th1))+&
  th1*td12*u2**(th1-1.0d0))!correct
  xder2uth=-delta*td11*u1**th1&
  -delta*theta*(-(delta-1.0d0)*td21*u1**(-theta)*log(u1)*u1**th1&
  -u1**th1*log(u1)*td11)!correct
  yder2vth=-delta*td12*u2**th1-&
  delta*theta*(-(delta-1.0d0)*td22*u2**(-theta)*log(u2)*u2**th1&
  -u2**th1*log(u2)*td12)!correct
  xder1th=delta*(td11*t1*tu1)
  yder1th=delta*(td12*t2*tu2)
  xder2udl=-theta*u1**th1*(td11+delta*td11*ttu1)
  yder2vdl=-theta*u2**th1*(td12+delta*td12*ttu2)
  
  mder1u=xder1u*s**del1/delta !correct
  mder1v=yder1v*s**del1/delta !correct
  mder2uv=del1*(s**(del1-1))*yder1v*xder1u/delta!correct
  mder2u=del1*s**(del1-1.0d0)*xder1u**2/delta+xder2u*s**del1/delta!correct
  mder2v=del1*(s**(del1-1.0d0))*(yder1v)**2/delta+yder2v*s**(del1)/delta!correct
  mder2uth=(xder2uth*s**del1+del1*s**(del1-1.0d0)*sder1th*xder1u)/delta !correct
  mder2udl=xder2udl*s**del1/delta+(-s**del1/delta**2+&
  (1.0d0/delta)*s**(del1)*(-log(s)/delta**2+sder1dl*del1/s))*xder1u
  mder2vth=(yder2vth*s**del1+del1*s**(del1-1.0d0)*sder1th*yder1v)/delta!correct
  mder2vdl=yder2vdl*s**del1/delta+(-s**del1/delta**2+&
  (1.0d0/delta)*s**(del1)*(-log(s)/delta**2+sder1dl*del1/s))*yder1v
  
  mderu(1)=mder1u;mderu(2)=mder2u;mderu(3)=mder2uth;mderu(4)=mder2udl;
  mderv(1)=mder1v;mderv(2)=mder2v;mderv(3)=mder2vth;mderv(4)=mder2vdl;

  return
  
  end

! this routine is in bb1facts.f90
!Function with bb8 logpdf =log f(u1,u2,th,delta) 
!lder1, lder2 (partial wrt cparv, 1st and 2nd order) lder1 is 2-dimenisonal vector
!lder2 is 
 !also lderu (partial wrt u1, 1st and 2nd order and also wrt parmaters)
! ***  also lderv (partial wrt u2, 1st order and 2nd order and also wrt parameters) and others
! and lderuv  (partial wrt u and v) 
! outputs 
!lpdf = log pdf, 
!lder11 = \p lpdf/\p nu  \p lpdf/\p delta !2-dimensional vector
!lder22 = \p^2 lpdf/\p^2 th^2, \p^2 lpdf/\p^2 th delta, lpdf/\p^2 delta^2!  3d-vector
!lderu = \p lpdf/\p u, \p^2 lpdf/\p u^2, \p^2 lpdf/\p u th, \p^2 lpdf/\p u delta  4d-vector
!lderv = \p lpdf/\p v, \p^2 lpdf/\p v^2, \p^2 lpdf/\p v th, \p^2 lpdf/\p v delta  4d-vector
!lderuv= \p^2 lpdf/\p uv
subroutine lbb1derivs(u1,u2,cparv,lpdf,lder11,lder22,lderu,lderv,lderuv)
  implicit none
  double precision u1,u2,cparv(2),lpdf,lder11(2),lder22(3)
  double precision theta,delta,der1th,der1dl,der2th,der2dl,derthd
  double precision t10,t1der1th,t1der1dl,t1der2th,t1der2dl,t1derthd
  double precision t20,t2der1th,t2der1dl,t2der2th,t2der2dl,t2derthd
  double precision t30,t3der1th,t3der1dl,t3der2th,t3der2dl,t3derthd
  double precision t40,t4der1th,t4der1dl,t4der2th,t4der2dl,t4derthd
  double precision m,mder1th,mder1dl,mder2th,mder2dl,mderthd,mderu(4),mderv(4)
  double precision den,coef,t1,t2,tu1,tu2
  double precision mp1,msq,thsq,den2,thtem,dltem,dl1,t1a,t2a,smlog
  double precision mder1u,mder2u,mder2uth,mder2udl,mder1v,mder2v,mder2vth,mder2vdl,mder2uv
  double precision lderu(4),lderv(4),lderuv,tem,tem1u,tem1v,lpdf1u,lpdf2u,lpdf2uth,lpdf2udl
  double precision lpdf2uv,lpdf1v,lpdf2v,lpdf2vth,lpdf2vdl,dl2,th1
  double precision tem1,tem2,tem3,th2,tem1th,tem1dl,tem2uth,tem2vth
  theta=cparv(1); delta=cparv(2);

  
  call  mderivs(u1,u2,cparv,m,mder1th,mder1dl,mder2th,mder2dl,mderthd,mderu,mderv,mder2uv)
  !print*,mderu
  !print*,mderv
  mder1u=mderu(1);mder2u=mderu(2);mder2uth=mderu(3);mder2udl=mderu(4)
  mder1v=mderv(1);mder2v=mderv(2);mder2vth=mderv(3);mder2vdl=mderv(4);
  mp1=1.d0+m; msq=m*m; thsq=theta*theta
  thtem=2.d0+1.d0/theta
  t10 = -(thtem)*log(mp1)
  t1der1th = log(mp1)/thsq - (thtem)*mder1th/mp1
  t1der1dl = -(thtem)*mder1dl/mp1
  t1der2th = -2.d0*log(mp1)/theta**3+(2.d0/thsq)*mder1th/mp1
  t1der2th = t1der2th-(thtem)*(mder2th/mp1-mder1th*mder1th/(mp1*mp1))
  t1der2dl = -(thtem)*(mder2dl/mp1-mder1dl*mder1dl/(mp1*mp1))
  t1derthd = mder1dl/(thsq*(mp1))-(thtem)*(mderthd/(mp1)-mder1th*mder1dl/(mp1*mp1))
  
  tem=theta*(delta-1.0d0)+(theta*delta+1.d0)*m
  tem1u=(theta*delta+1.0d0)*mder1u
  tem1v=(theta*delta+1.0d0)*mder1v
  tem1th=delta-1.0d0+mder1th+delta*(m+mder1th*theta)
  tem1dl=theta+mder1dl+theta*(m+mder1dl*delta)
  tem2uth=delta*mder1u+mder2uth*(theta*delta+1.0d0)
  tem2vth=delta*mder1v+mder2vth*(theta*delta+1.0d0)
 !  print*,tem
!   print*,tem1th
!   print*,tem2uth
  dltem=1.d0-2.d0*delta; dl1=delta-1.d0
  t20 = (dltem)*log(m)
  t2der1th = (dltem)*mder1th/m
  t2der1dl = -2.d0*log(m)+(dltem)*mder1dl/m
  t2der2th = (dltem)*(mder2th/m-mder1th*mder1th/msq)
  t2der2dl = -4.d0*mder1dl/m+(dltem)*(mder2dl/m-mder1dl*mder1dl/msq)
  t2derthd = -2.d0*mder1th/m+(dltem)*(mderthd/m-mder1th*mder1dl/msq)

  coef = theta*delta+1.d0
  den = theta*(dl1)+(coef)*m
  den2=den*den
  t30 = log(theta*(dl1)+coef*m)
  t3der1th = (dl1+delta*m+coef*mder1th)/den
  t3der1dl = (theta+theta*m+coef*mder1dl)/den
  t3der2th = (2.d0*delta*mder1th+coef*mder2th)/den
  t3der2th = t3der2th-(dl1+delta*m+coef*mder1th)**2/den2
  t3der2dl = (2.d0*theta*mder1dl+coef*mder2dl)/den
  t3der2dl = t3der2dl-(theta+theta*m+coef*mder1dl)**2/den2;
  t3derthd = (1.d0+m+theta*mder1th+delta*mder1dl+coef*mderthd)/den;
  t3derthd = t3derthd - (theta+theta*m+coef*mder1dl)*(dl1+delta*m+coef*mder1th)/den2

  t1 = u1**(-theta); t2 = u2**(-theta)
  t1a=t1-1.d0; t2a=t2-1.d0; smlog=log(t1a)+log(t2a)
  tu1 = -log(u1); tu2 = -log(u2)
  t40 = (dl1)*(smlog)+(theta+1.d0)*(tu1+tu2)
  t4der1th = (dl1)*(t1*tu1/t1a+t2*tu2/t2a)+tu1+tu2  
  t4der1dl = smlog
  t4der2th = -(dl1)*(t1*tu1*tu1/t1a**2+t2*tu2*tu2/t2a**2)
  t4der2dl = 0.d0
  t4derthd = t1*tu1/t1a+t2*tu2/t2a

  lpdf = t10+t20+t30+t40
  der1th = t1der1th+t2der1th+t3der1th+t4der1th
  der1dl = t1der1dl+t2der1dl+t3der1dl+t4der1dl
  der2th = t1der2th+t2der2th+t3der2th+t4der2th
  der2dl = t1der2dl+t2der2dl+t3der2dl+t4der2dl
  derthd = t1derthd+t2derthd+t3derthd+t4derthd
  lder11(1)=der1th; lder11(2)=der1dl  ! gradient of log pdf
  lder22(1)=der2th; lder22(2)=derthd; lder22(3)=der2dl; ! Hessian terms
  
  th2=(-1.0d0/theta)-2.0d0
  dl2=1.0d0-2.0d0*delta
  th1=-theta-1.0d0
  
  tem1=th2*mder1u/(m+1.0d0)+(dl2)*mder1u/m
  tem2=coef*mder1u/tem+dl1*(-theta*u1**(-theta-1.0d0)/(u1**(-theta)-1.0d0))&
  -(theta+1.0d0)/u1
  lpdf1u=tem1+tem2 !correct
  
  tem1=th2*mder1v/(m+1.0d0)+(dl2)*mder1v/m
  tem2=coef*mder1v/tem+dl1*(-theta*u2**(-theta-1.0d0)/(u2**(-theta)-1.0d0))&
  -(theta+1.0d0)/u2
  lpdf1v=tem1+tem2 !correct
  
  tem1=th2*(-mder1u**2.0d0/(m+1.0d0)**2.0d0+mder2u/(m+1.0d0))&
  +dl2*(-mder1u**2.0d0/m**2+mder2u/m)
  tem2=coef*(-tem1u*mder1u/tem**2.0d0+mder2u/tem)+(theta+1.0d0)/u1**2.0d0
  tem3=(dl1/(u1**(-theta)-1)**2.d0)*(-theta*th1*u1**(th1-1.0d0)*(u1**(-theta)-1.0d0)&
  -(theta*u1**th1)**2.0d0)
  lpdf2u=tem1+tem2+tem3 !correct
  
  
  tem1=th2*(-mder1v**2.0d0/mp1**2.0d0+mder2v/mp1)+dl2*(-mder1v**2.0d0/m**2+mder2v/m)
  tem2=coef*(-tem1v*mder1v/tem**2.0d0+mder2v/tem)+(theta+1.0d0)/u2**2.0d0
  tem3=(dl1/(u2**(-theta)-1)**2.d0)*(-theta*th1*u2**(th1-1.0d0)*(u2**(-theta)-1.0d0)&
  -(theta*u2**th1)**2.0d0)
  lpdf2v=tem1+tem2+tem3 !correct
  
  tem1=th2*(-mder1u*mder1v/mp1**2+mder2uv/mp1)+dl2*(-mder1u*mder1v/m**2+mder2uv/m)
  tem2=coef*(-tem1v*mder1u/tem**2+mder2uv/tem)
  lpdf2uv=tem1+tem2!correct
  
  tem1=mder1u/(mp1*thsq)+(-1.0d0/theta-2.0d0)*(-mder1th*mder1u/mp1**2+mder2uth/mp1)&
  +dltem*(-mder1th*mder1u/m**2+mder2uth/m)
  tem2=(-tem1th*tem1u/tem/tem)+tem2uth/tem
  tem3=(-t1/u1+t1*log(u1)*theta/u1)*t1a-t1*log(u1)*theta*t1/u1
  lpdf2uth=tem1+tem2+dl1*tem3/t1a**2-1.0d0/u1
  
  tem1=mder1v/(mp1*thsq)+(-1.0d0/theta-2.0d0)*(-mder1th*mder1v/mp1**2+mder2vth/mp1)&
  +dltem*(-mder1th*mder1v/m**2+mder2vth/m)
  tem2=(-tem1th*tem1v/tem/tem)+tem2vth/tem
  tem3=(-t2/u2+t2*log(u2)*theta/u2)*t2a-t2*log(u2)*theta*t2/u2
  lpdf2vth=tem1+tem2+dl1*tem3/t2a**2-1.0d0/u2
  
  tem1=th2*(-mder1dl*mder1u/mp1**2+mder2udl/mp1)-2.0d0*mder1u/m
  tem2=dl2*(-mder1dl*mder1u/m**2+mder2udl/m)-tem1dl*coef*mder1u/tem**2.0d0
  tem3=(theta*mder1u+coef*mder2udl)/tem
  lpdf2udl=tem1+tem2+tem3-theta*(t1/u1)/(t1-1.0d0)!correct
  
  tem1=th2*(-mder1dl*mder1v/mp1**2+mder2vdl/mp1)-2.0d0*mder1v/m
  tem2=dl2*(-mder1dl*mder1v/m**2+mder2vdl/m)-tem1dl*coef*mder1v/tem**2.0d0
  tem3=(theta*mder1v+coef*mder2vdl)/tem
  lpdf2vdl=tem1+tem2+tem3-theta*(t2/u2)/(t2-1.0d0)!correct
  
  lderu(1)=lpdf1u;lderu(2)=lpdf2u;lderu(3)=lpdf2uth;lderu(4)=lpdf2udl;
  lderv(1)=lpdf1v;lderv(2)=lpdf2v;lderv(3)=lpdf2vth;lderv(4)=lpdf2vdl;
  
  lderuv=lpdf2uv
  return
end

! Function with Gumbel log pdf c(u1,u2,cpar) for tree 2, and 
!   lder1, lder2 (partial wrt cpar, 1st and 2nd order)
!   also lder1u lder2u (partial wrt u1, 1st and 2nd order)
! ***  also lder1v (partial wrt u2, 1st order) and others
!   and ldermixu  (partial wrt cpar and u1) 
! inputs
!   u1,u2 = values in (0,1)
!   cpar = scalar in (1,100)
! outputs 
!   lpdf = log pdf, 
!   lder1 = \p lpdf/\p cpar
!   lder2 = \p^2 lpdf/\p cpar^2
!   lder1u = \p lpdf/\p u1
!   lder2u = \p^2 lpdf/\p u1^2
!   ldermixu = \p^2 lpdf/\p u1 \p cpar
!   lder1v = \p lpdf/\p u2
!   lder2v = \p^2 lpdf/\p u2^2
!   lder2uv = \p^2 lpdf/\p u1 \p u2
!   ldermixv = \p^2 lpdf/\p u2 \p cpar
!  some derivs added; function name with ending of 't'
subroutine  lgum2derivt(u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,&
ldermixu,lder1v,lder2v,lder2uv,ldermixv)
  implicit none
  double precision u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,lder2uv
  double precision x,y,tx,ty,xd,yd,s,m,logs,logm,msq,dlsq,dlcu
  double precision sder1,sder2,mder1,mder2,den,den2
  double precision mu,m2u,u1sq,muder1,u2sq,lder1v,lder2v,ldermixv,ldermixu
  double precision mv,m2v,mvder1,m2uv
  x = -log(u1); y = -log(u2);
  tx = log(x); ty = log(y); 
  xd = x**cpar; yd = y**cpar; s = xd+yd; m = s**(1.d0/cpar);
  logs=log(s);
  dlsq=cpar*cpar; dlcu=dlsq*cpar
  ! for 1-factor and 2-factor models
  sder1 = xd*tx+yd*ty;
  sder2 = xd*tx*tx+yd*ty*ty;
  mder1 = m*sder1/(s*cpar)-m*logs/dlsq;
  mder2 = -mder1*logs/dlsq-2.d0*m*sder1/(s*dlsq)+2.d0*m*logs/dlcu;
  mder2 = mder2+sder2*m/(s*cpar)+(mder1/s-m*sder1/(s*s))*sder1/cpar;
  den = m+cpar-1.d0; den2=den*den
  logm=log(m); msq=m*m
  lpdf = -m+log(den)+(1.d0-2.d0*cpar)*logm+(cpar-1.d0)*(tx+ty)+x+y;
  lder1 = -mder1+(mder1+1.d0)/den-2.d0*logm+(1.d0-2.d0*cpar)*mder1/m+tx+ty;
  lder2 = -mder2+mder2/den-(mder1+1.d0)**2/den2-4.d0*mder1/m&
  +(1.d0-2.d0*cpar)*(mder2/m-mder1*mder1/msq);
  ! for 2-factor model
  u1sq=u1*u1;
  mu = -m*xd/(u1*s*x); 
  m2u = (1.d0-cpar)*m*xd*xd/(u1*s*x)**2+(cpar-1.d0)*m*xd/(s*x*x*u1sq)+m*xd/(s*x*u1sq);
  muder1 = -(mder1/s-m*sder1/(s*s))*xd/(x*u1)-m*xd*tx/(s*u1*x);
  lder1u = -mu+mu/den+(1.d0-2.d0*cpar)*mu/m-(cpar-1.d0)/(u1*x)-1.d0/u1;
  lder2u = -m2u+m2u/den-mu*mu/den2+(1.d0-2.d0*cpar)*(m2u/m-mu*mu/msq)+(cpar-1.d0)/(u1sq*x); 
  lder2u = lder2u-(cpar-1.d0)/(x*x*u1sq)+1.d0/u1sq;
  ldermixu = -muder1+muder1/den-mu*(mder1+1.d0)/den2-2.d0*mu/m&
  +(1.d0-2.d0*cpar)*(muder1/m-mu*mder1/msq)-1.d0/(x*u1);
  
  u2sq=u2*u2;
  mv = -m*yd/(u2*s*y); 
  m2v = (1.d0-cpar)*m*yd*yd/(u2*s*y)**2+(cpar-1.d0)*m*yd/(s*y*y*u2sq)+m*yd/(s*y*u2sq);
  mvder1 = -(mder1/s-m*sder1/(s*s))*yd/(y*u2)-m*yd*ty/(s*u2*y);
  lder1v = -mu+mu/den+(1.d0-2.d0*cpar)*mu/m-(cpar-1.d0)/(u2*y)-1.d0/u2;
  lder2v = -m2u+m2u/den-mu*mu/den2+(1.d0-2.d0*cpar)*(m2u/m-mu*mu/msq)+(cpar-1.d0)/(u2sq*y); 
  lder2v = lder2v-(cpar-1.d0)/(y*y*u2sq)+1.d0/u2sq;
  ldermixv = -muder1+muder1/den-mu*(mder1+1.d0)/den2-2.d0*mu/m&
  +(1.d0-2.d0*cpar)*(muder1/m-mu*mder1/msq)-1.d0/(y*u2);
  
  m2uv=(x**(cpar-1.0d0)/u1/u2)*(cpar*(1.0d0/cpar-1.0d0)*m*y**(cpar-1.0d0)/s/s)
  lder2uv=-m2uv+(-mu*mv/den**2+m2uv/den)+(1.0d0-2*cpar)*(m2uv/m-mu*mv/m**2)
  return
  end
  


! Function with t log pdf = log f(t1,t2,rho,nu) for factor 2, and
!   lder1, lder2 (partial wrt rho, 1st and 2nd order)
!   also lder1u lder2u (partial wrt t1, 1st and 2nd order)
!   and ldermix  (partial wrt rho and t1 )
! inputs
!   t1 = qt(u1,nu) transformed to t scale
!   t2 = qt(u2,nu) transformed to t scale
!   rho = scalar in (-1,1)
!   nu = degree of freedom >0
! outputs 
!   lpdf = log pdf, 
!   lder1 = \p lpdf/\p rho
!   lder2 = \p^2 lpdf/\p rho^2
!   ltder1u = \p lpdf/\p t1
!   ltder2u = \p^2 lpdf/\p t1^2
!   ltdermix = \p^2 lpdf/\p u1 \p rho
subroutine lt2derivs(t1,t2,rho,nu,lpdf,lder1,lder2,ltder1u,ltder2u,ltdermix)
  implicit none
  double precision t1,t1sq,t2,t2sq,rho,nu,lpdf,lder1,lder2,ltder1u,ltder2u,ltdermix
  double precision pi,con,lgdif,coef,coef2,coef3,dt1,dt2,den,den1,den2
  double precision dd1,dd12,ltder1t1,ltder2t1,ltdert1mx, reg
  double precision lgamma   ! function in C code lgamma.c, in gnu fortran
  pi=3.14159265358979323846d0; ! pi
  con=0.5d0*log(pi*nu);
  !lgdif=log(gamma(0.5*(nu+1)))-log(gamma(0.5*nu));
  lgdif=lgamma(0.5d0*(nu+1.d0))-lgamma(0.5d0*nu);
  coef = 1.d0-rho*rho; coef2 = coef*coef; 
  t1sq = t1*t1; t2sq = t2*t2;
  den = 1.d0+(t1sq-2.d0*rho*t1*t2+t2sq)/(nu*coef);
  den1 = 2.d0*rho*(den-1.d0)/coef-2.d0*t1*t2/(nu*coef);
  den2 = 2.d0*rho*den1/coef+2.d0*(den-1.d0)*(2.d0-coef)/coef2-4.d0*rho*t1*t2/(nu*coef2);
  dt1 = lgdif-con-0.5d0*(nu+1.d0)*log(1.d0+t1sq/nu);
  dt2 = lgdif-con-0.5d0*(nu+1.d0)*log(1.d0+t2sq/nu);
  dd1 = exp(dt1); dd12 = dd1*dd1; 
  reg = t1-rho*t2;  ! added temporary variable
  coef3 = den*coef*nu; ! moved
  lpdf = -log(2.d0*pi)-0.5d0*log(coef)-0.5d0*(nu+2.d0)*log(den)-dt1-dt2;
  lder1 = rho/coef-0.5d0*(nu+2.d0)*den1/den;
  lder2 = (2.d0-coef)/coef2-0.5d0*(nu+2.d0)*(den2/den-den1*den1/(den*den));
  !ltder1t1 = -(nu+2)*(t1-rho*t2)/coef3 + (nu+1)*t1/(nu+t1sq);
  !ltder2t1 = -(nu+2)/coef3+2*(nu+2)*(t1-rho*t2)**2/coef3**2+(nu+1)*(nu-t1sq)/(nu+t1sq)**2;
  !ltdert1mx = (nu+2)*(t2/coef3-2*rho*(t1-rho*t2)/coef3/coef+(t1-rho*t2)*den1/den/coef3);
  ltder1t1 = -(nu+2.d0)*reg/coef3 + (nu+1.d0)*t1/(nu+t1sq);
  ltder2t1 = -(nu+2.d0)/coef3+2.d0*(nu+2.d0)*(reg*reg)/coef3**2+&
  (nu+1.d0)*(nu-t1sq)/(nu+t1sq)**2;
  ltdert1mx = (nu+2.d0)*(t2/coef3-2.d0*rho*reg/coef3/coef+reg*den1/den/coef3);
  ltder1u =ltder1t1/dd1;
  ltder2u =ltder2t1/dd12+(nu+1.d0)*t1*ltder1t1/dd12/(nu+t1sq);
  ltdermix=ltdert1mx/dd1;
  return
  end
  
! Function with Frank log pdf c(u1,u2,cpar) for tree 2, and 
!   lder1, lder2 (partial wrt cpar, 1st and 2nd order)
!   also lder1u lder2u (partial wrt u1, 1st and 2nd order)
! ***  also lder1v (partial wrt u2, 1st order) and others
!   and ldermixu  (partial wrt cpar and u1) 
! inputs
!   u1,u2 = values in (0,1)
!   cpar = scalar in (-oo,oo)
! outputs 
!   lpdf = log pdf, 
!   lder1 = \p lpdf/\p cpar
!   lder2 = \p^2 lpdf/\p cpar^2
!   lder1u = \p lpdf/\p u1
!   lder2u = \p^2 lpdf/\p u1^2
!   ldermixu = \p^2 lpdf/\p u1 \p cpar
!   lder1v = \p lpdf/\p u2
!   lder2v = \p^2 lpdf/\p u2^2
!   lder2uv = \p^2 lpdf/\p u1 \p u2
!   ldermixv = \p^2 lpdf/\p u2 \p cpar
!  some derivs added; function name with ending of 't'
subroutine lfrk2derivt(u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,ldermixu,&
 lder1v,lder2v,lder2uv,ldermixv)
  implicit none
  double precision u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,ldermixu
  double precision lder1v,lder2v,lder2uv,ldermixv
  double precision den,den1,den2,den1u,den2u,denmixu,t0,t1,t2
  double precision den1v,den2v,denmixv

  t0 = exp(-cpar);
  t1 = exp(-cpar*u1);
  t2 = exp(-cpar*u2);
  den = t1+t2-t0-t1*t2;
  den1 = -u1*t1-u2*t2+t0+(u1+u2)*t1*t2;
  den2 = u1*u1*t1+u2*u2*t2-t0-(u1+u2)*(u1+u2)*t1*t2;
  den1u = -cpar*t1*(1.d0-t2);  
  den1v = -cpar*t2*(1.d0-t1);   ! added
  den2u = cpar*cpar*t1*(1.d0-t2);
  den2v = cpar*cpar*t2*(1.d0-t1); ! added
  denmixu = t1*(-1.d0+cpar*u1+t2-(u1+u2)*cpar*t2); 
  denmixv = t2*(-1.d0+cpar*u2+t1-(u1+u2)*cpar*t1); ! added 
   
  !lpdf = log(abs(cpar))+log(abs(1-t0))-cpar*(u1+u2)-2*log(abs(den));
  ! 121023 maybe later add the limits as cpar->0
  ! pdf = cpar*(1-t0)/den^2 where
  !    1-t0 has same sign as cpar,  den has same sign as cpar
  lpdf = log(cpar*(1.d0-t0))-cpar*(u1+u2)-2.d0*log(abs(den));
  lder1 = 1.d0/cpar+t0/(1.d0-t0)-(u1+u2)-2.d0*den1/den;
  lder2 = -1.d0/(cpar*cpar)-t0/((1.d0-t0)*(1.d0-t0))-2.d0*den2/den&
  +2.d0*den1*den1/(den*den); 
  lder1u = -cpar-2.d0*den1u/den;
  lder1v = -cpar-2.d0*den1v/den;  ! added
  lder2u = -2.d0*den2u/den+2.d0*den1u*den1u/(den*den);
  lder2v = -2.d0*den2v/den+2.d0*den1v*den1v/(den*den); ! added
  lder2uv = 2.d0*cpar*cpar*t1*t2/den+2.d0*den1u*den1v/(den*den); ! added
  ldermixu = -1.d0-2.d0*denmixu/den+2.d0*den1u*den1/(den*den);
  ldermixv = -1.d0-2.d0*denmixv/den+2.d0*den1v*den1/(den*den); ! added
  return
  end
  
! Function with bivariate Gaussian copula log pdf c(u1,u2,cpar) for tree 2, and 
!   lder1, lder2 (partial wrt cpar, 1st and 2nd order)
!   also lder1u lder2u (partial wrt u1, 1st and 2nd order)
! ***  also lder1v (partial wrt u2, 1st order) and others
!   and ldermixu  (partial wrt cpar and u1) 
! inputs
!   u1,u2 = values in (0,1)
!   cpar = scalar in (-oo,oo)
! outputs 
!   lpdf = log pdf, 
!   lder1 = \p lpdf/\p cpar
!   lder2 = \p^2 lpdf/\p cpar^2
!   lder1u = \p lpdf/\p u1
!   lder2u = \p^2 lpdf/\p u1^2
!   ldermixu = \p^2 lpdf/\p u1 \p cpar
!   lder1v = \p lpdf/\p u2
!   lder2v = \p^2 lpdf/\p u2^2
!   lder2uv = \p^2 lpdf/\p u1 \p u2
!   ldermixv = \p^2 lpdf/\p u2 \p cpar
!  some derivs added; function name with ending of 't'
subroutine lgau2derivt(u1,u2,rho,lpdf,lder1,lder2,lder1u,lder2u,ldermixu,&
lder1v,lder2v,lder2uv,ldermixv)
  implicit none
  double precision u1,u2,rho,lpdf,lder1,lder2,lder1u,lder2u,ldermixu
  double precision lder1v,lder2v,lder2uv,ldermixv
  double precision qnorms, dnorms
  double precision x,y,x2,y2,rh2,den,den2,con,qf1,qf2,dx,dy,dx2,dy2

  x=qnorms(u1); x2=x*x
  y=qnorms(u2); y2=y*y
  rh2=rho**2; den=1.d0-rh2; den2=den*den
  con= -0.5d0*log(den)
  qf1= -(x2+y2-2.d0*rho*x*y)/(2*den)
  qf2=(x2+y2)/2.d0
  lpdf=con+qf1+qf2
  lder1 = rho/den+2.d0*qf1*rho/den + x*y/den
  dx=dnorms(x); dy=dnorms(y); dx2=dx*dx; dy2=dy*dy
  lder1u = (x-(x-rho*y)/den)/dx
  lder1v = (y-(y-rho*x)/den)/dy
  lder2u = (rho*(y-rho*x)/den)*x/dx2 - rh2/den/dx2
  lder2v = (rho*(x-rho*y)/den)*y/dy2 - rh2/den/dy2
  lder2uv = rho/den/dx/dy
  lder2=(1.d0+rh2)/den2 + (qf2*2.d0*(-1.d0-3.d0*rh2)+2.d0*rho*(3.d0+rh2)*x*y)/den2/den
  ldermixu=((1.d0+rh2)*y-2.d0*rho*x)/den2/dx
  ldermixv=((1.d0+rh2)*x-2.d0*rho*y)/den2/dy
  return
  end

!Function with bb8 logpdf =log f(u1,u2,nu,delta) 
!lder1, lder2 (partial wrt rho, 1st and 2nd order)
 !also lderu (partial wrt u1, 1st and 2nd order and also wrt parmaters)
! ***  also lderv (partial wrt u2, 1st order and 2nd order and also wrt parameters) and others
! and lderuv  (partial wrt u and v) 
! outputs 
!lpdf = log pdf, 
!lder11 = \p lpdf/\p nu  \p lpdf/\p delta !2-dimensional vector
!lder22 = \p^2 lpdf/\p^2 nu^2 \p^2 \p^2 lpdf/\p^2 nu delta lpdf/\p^2 delta^2  !3-dimensional vector
!lderu = \p lpdf/\p u, \p^2 lpdf/\p u^2, \p^2 lpdf/\p u nu, \p^2 lpdf/\p u delta 
!lderv = \p lpdf/\p v, \p^2 lpdf/\p v^2, \p^2 lpdf/\p v nu, \p^2 lpdf/\p v delta 
!lderuv= \p^2 lpdf/\p uv
subroutine lbb8derivs(u1,u2,cparv,lpdf,lder11,lder22,lderu,lderv,lderuv)
 implicit none
 double precision u1,u2,cparv(2),lpdf,lder11(2),lder22(3)
 double precision eta,x,y,delta,dl1,dlu1,dlu2,nu,nu1,nu2
 double precision eta1nu,eta1del,x1nu,x1del,y1nu,y1del,eta2nu,eta2del,eta2nudel
 double precision tem,temm,tem1,tem2,tem3,tem4,lpdf1nu,lpdf1del,tem1nu,tem1dl,eta1dl
 double precision x2nu,y2nu,y2nudel,y2del,x2nudel,x2del,lpdf2nu,lpdf2del,lpdf2nudel
 double precision xy1mixdel,tem2nu,tem2del,tem2nudel
 double precision lderu(4),lderv(4),lderuv,lpdf1u,lpdf2u,lpdf1v,lpdf2v
 double precision lpdf2udl,lpdf2unu,lpdf2vdl,lpdf2vnu
 double precision tem1u,tem1v,x1u,y1v,tem2u,tem2v,x2u,y2v,x2unu,y2vnu
 double precision x2udl,y2vdl,tem2unu,tem2udl,tem2vnu,tem2vdl
 double precision tem2uv
 nu=cparv(1)
 delta=cparv(2)
 dl1=1.0d0-delta
 nu1=nu-1.0d0
 nu2=1.0d0/nu-2.0d0
 dlu1=delta*u1
 dlu2=delta*u2
 eta=1.0d0-dl1**nu
 x=1.0d0-(1.0d0-delta*u1)**nu
 y=1.0d0-(1.0d0-delta*u2)**nu
 tem=x*y/eta !x*y/eta
 temm=1.0d0-tem
 lpdf=-log(eta)+log(delta)+(nu2)*log(temm)+log(nu-tem)+(nu1)*log(1.0d0-dlu1)+&
 (nu1)*log(1-dlu2)
 
 eta1nu=-log(dl1)*((dl1)**nu)!\p eta/\p nu
 eta1del=nu*((dl1)**(nu1)) !\p eta/\p delta
 x1nu=-log(1.0d0-dlu1)*((1-dlu1)**nu) !\p x/\p nu
 x1del=(nu*(1.0d0-dlu1)**(nu1))*u1 !\p x/\p delta
 y1nu=-log(1.0d0-dlu2)*((1-delta*u2)**nu)  !\p y/\p nu
 y1del=nu*((1.0d0-dlu2)**(nu1))*u2 !\p y/\p delta
 tem1nu=-(x*y/eta**2)*eta1nu+(x1nu*y+y1nu*x)/eta !\p tem/\p nu
 tem1dl=(x1del*y+y1del*x)/eta-eta1del*x*y/eta**2!\p tem/\p delta
 
 tem1=-eta1nu/eta-log(temm)/nu**2+nu2*(1.0d0/(temm))*(-tem1nu)
 tem2=(1.0d0/(nu-tem))*(1.0d0-tem1nu)+log(1.0d0-dlu1)+log(1.0d0-dlu2)
 lpdf1nu=tem1+tem2!\p lpdf/\p nu
 
 tem1=-eta1del/eta+1.0d0/delta+nu2*(1.0d0/temm)*(-tem1dl)
 tem2=-(1.0d0/(nu-tem))*(tem1dl)-(nu1)*u1/(1.0d0-delta*u1)-nu1*u2/(1.0d0-delta*u2)
 lpdf1del=tem1+tem2 !\p lpdf/\p delta
 
 lder11(1)=lpdf1nu;lder11(2)=lpdf1del;
 
 eta2nu=-((log(dl1))**2)*(dl1**nu) !\p^2 eta/\p nu^2
 eta2del=-nu*(nu1)*(dl1)**(nu-2.0d0) !\p^2 eta/\p delta^2
 eta2nudel=nu*log(dl1)*(dl1**nu1)+dl1**nu1!\p^2 eta/\p nu delta
 x2nu=-(log(1.0d0-dlu1)**2)*((1.0d0-dlu1)**nu) !\x^2 eta/\p nu^2
 y2nu=-(log(1.0d0-delta*u2)**2)*((1.0d0-dlu2)**nu)!\y^2 eta/\p nu^2
 x2del=-u1*u1*nu*nu1*((1-dlu1)**(nu-2.0d0))!\x^2 eta/\p delta^2
 y2del=-u2*u2*nu*nu1*((1-dlu2)**(nu-2.0d0))!\y^2 eta/\p delta^2
 x2nudel=u1*nu*(1-dlu1)**(nu1)*log(1-dlu1)+u1*((1-dlu1)**nu)/(1-dlu1) !\x^2 eta/\p nu delta
 y2nudel=u2*nu*(1-dlu2)**(nu1)*log(1-dlu2)+u2*((1-dlu2)**nu)/(1-dlu2)!\y^2 eta/\p nu delta
 
 tem1=(2.0d0/eta**3)*x*y*((eta1nu)**2)-(eta2nu*x*y+(x1nu*y+y1nu*x)*eta1nu)/eta**2
 tem2=-eta1nu*x1nu*y/eta**2+(x2nu*y+x1nu*y1nu)/eta
 tem3=-eta1nu*y1nu*x/eta**2+(y2nu*x+x1nu*y1nu)/eta
 tem2nu=tem1+tem2+tem3 !\p^2 tem/\p nu^2
 
 xy1mixdel=x1del*y1del
 tem1=(2.0d0*x*y*((eta1del)**2))/eta**3-(eta2del*x*y+(x1del*y+y1del*x)*eta1del)/eta**2
 tem2=-eta1del*y*x1del/eta**2+(xy1mixdel+x2del*y)/eta
 tem3=-eta1del*x*y1del/eta**2+(xy1mixdel+y2del*x)/eta
 tem2del=tem1+tem2+tem3!\p^2 tem/\p delta^2
 
 
 tem1=2.0d0*eta1nu*eta1del*x*y/eta**3.0d0
 tem2=-(eta2nudel*x*y+(x1del*y+y1del*x)*eta1nu)/eta**2.0d0
 tem3=-eta1del*x1nu*y/eta**2+(x2nudel*y+y1del*x1nu)/eta
 tem4=-eta1del*y1nu*x/eta**2+(y2nudel*x+x1del*y1nu)/eta
 tem2nudel=tem1+tem2+tem3+tem4 !\p^2 tem/\p nu delta
 
 tem1=eta1nu**2/eta**2-eta2nu/eta+2.0d0*log(temm)/nu**3.0d0
 tem2=tem1nu/(temm)/nu**2.0d0+tem1nu/(temm)/nu**2.0d0
 tem3=(1.0d0/nu-2.0d0)*(-(tem1nu**2.0d0/(temm)**2.0d0)-tem2nu/(temm))
 tem4=-(1.0d0-tem1nu)**2.0d0/(nu-tem)**2.0d0-tem2nu/(nu-tem)
 lpdf2nu=tem1+tem2+tem3+tem4!\p^2 lpdf/\p nu^2 
 
 tem1=eta1del**2/eta**2-eta2del/eta-1.0d0/delta**2
 !tem2=(-(1.0d0/(1.0d0-tem)**2)*(tem1del)**2+(-tem2del)/(1.0d0-tem)))*(1.0d0/nu-2.0d0)
 tem2=nu2*(-tem1dl**2/(temm)**2-tem2del/temm)
 tem3=-(tem1dl**2/(nu-tem)**2+tem2del/(nu-tem))
 tem4=-(nu1*u1**2.0d0/(1-dlu1)**2.0d0+nu1*u2**2.0d0/(1-dlu2)**2.0d0)
 lpdf2del=tem1+tem2+tem3+tem4!\p^2 lpdf/\p delta^2 
 
 tem1=eta1nu*eta1del/eta**2-eta2nudel/eta+tem1dl/nu**2/(1.0d0-tem)
 tem2=nu2*(-(tem1nu*tem1dl/(1-tem)**2)-tem2nudel/(1.0d0-tem))
 tem3=(1.0d0-tem1nu)*tem1dl/(nu-tem)**2-tem2nudel/(nu-tem)
 tem4=-u1/(1.0d0-dlu1)-u2/(1.0d0-dlu2)
 lpdf2nudel=tem1+tem2+tem3+tem4
 
 lder22(1)=lpdf2nu
 lder22(2)=lpdf2nudel
 lder22(3)=lpdf2del
 
 x1u=delta*nu*(1-dlu1)**nu1
 y1v=delta*nu*(1-dlu2)**nu1
 tem1u=y*x1u/eta
 tem1v=x*y1v/eta
 
 lpdf1u=-nu2*tem1u/temm-tem1u/(nu-tem)-nu1*delta/(1-dlu1)
 lpdf1v=-nu2*tem1v/temm-tem1v/(nu-tem)-nu1*delta/(1-dlu2)
 
 x2u=-(delta**2)*nu*nu1*(1-dlu1)**(nu-2.0d0)
 y2v=-(delta**2)*nu*nu1*(1-dlu2)**(nu-2.0d0)
 tem2u=y*x2u/eta
 tem2v=x*y2v/eta
 
 lpdf2u=-nu2*tem1u**2/temm**2-nu2*(tem2u/temm)-tem1u**2/(nu-tem)**2&
 -tem2u/(nu-tem)-nu1*delta**2/(1-dlu1)**2
 lpdf2v=-nu2*tem1v**2/temm**2-nu2*(tem2v/temm)-tem1v**2/(nu-tem)**2&
 -tem2v/(nu-tem)-nu1*delta**2/(1-dlu2)**2
 
 x2unu=delta*(1-dlu1)**nu1+delta*nu*log(1-dlu1)*(1-dlu1)**nu1
 tem2unu=-eta1nu*y*x1u/eta**2+(y1nu*x1u+x2unu*y)/eta 
 tem1=tem1u/nu**2/temm+(-tem1nu*tem1u/temm**2-tem2unu/temm)*nu2
 tem2=tem1u*(1-tem1nu)/(nu-tem)**2-tem2unu/(nu-tem)-delta/(1-dlu1)
 lpdf2unu=tem1+tem2
 
 eta1dl=nu*dl1**nu1
 x2udl=nu*(1.0d0-dlu1)**nu1-delta*nu*nu1*u1*(1.0d0-dlu1)**(nu1-1.0d0)
 tem2udl=-eta1dl*y*x1u/eta**2+(y1del*x1u+x2udl*y)/eta
 tem1=-nu2*tem1dl*tem1u/temm**2-nu2*tem2udl/temm
 tem2=-tem1dl*tem1u/(nu-tem)**2-tem2udl/(nu-tem)-nu1/(1-dlu1)**2
 lpdf2udl=tem1+tem2
 
 y2vnu=delta*(1-dlu2)**nu1+delta*nu*log(1-dlu2)*(1-dlu2)**nu1
 tem2vnu=-eta1nu*x*y1v/eta**2+(x1nu*y1v+y2vnu*x)/eta 
 tem1=tem1v/nu**2/temm+(-tem1nu*tem1v/temm**2-tem2vnu/temm)*nu2
 tem2=tem1v*(1-tem1nu)/(nu-tem)**2-tem2vnu/(nu-tem)-delta/(1-dlu2)
 lpdf2vnu=tem1+tem2
 
 y2vdl=nu*(1.0d0-dlu2)**nu1-delta*nu*nu1*u2*(1.0d0-dlu2)**(nu1-1.0d0)
 tem2vdl=-eta1dl*x*y1v/eta**2+(x1del*y1v+y2vdl*x)/eta
 tem1=-nu2*tem1dl*tem1v/temm**2-nu2*tem2vdl/temm
 tem2=-tem1dl*tem1v/(nu-tem)**2-tem2vdl/(nu-tem)-nu1/(1-dlu2)**2
 lpdf2vdl=tem1+tem2
 
 tem2uv=y1v*x1u/eta
 tem1=-nu2*tem1u*tem1v/temm**2-tem2uv*nu2/temm
 tem2=-tem1u*tem1v/(nu-tem)**2-tem2uv/(nu-tem)
 lderuv=tem1+tem2
 
 
 lderu(1)=lpdf1u;lderu(2)=lpdf2u;lderu(3)=lpdf2unu;lderu(4)=lpdf2udl;
 lderv(1)=lpdf1v;lderv(2)=lpdf2v;lderv(3)=lpdf2vnu;lderv(4)=lpdf2vdl;
return
end
 
