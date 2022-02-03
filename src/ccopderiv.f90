 ! program test
!    implicit none
!    double precision u1,u2,cparv(2)
!    double precision cder1vec(2),cder2vec(3)
!   double precision ccdf,cder1,cder2
!   integer icop
!    
!    	read*,u1
!   	read*,u2
!    	read*,cparv
!    	read*,icop
!    	
!    	print*,u1
!   	print*,u2
!    	print*,cparv
!    	print*,icop
!    	
!    call ccopderiv(u1,u2,cparv,icop,ccdf,cder1vec,cder2vec)
!    !call cgumderivs(u1,u2,cpar,ccdf,cder1,cder2)
!    print *, "conditional pdf"
!    print*,ccdf
!    print*,"der1 for par"
!    print*,cder1vec
!    print*,"der2 for par"
!    print*,cder2vec
!   
! stop 
! end



!This function computes the conditional cdf of bivariate copula pcond(u1|u2)
!and also the derivatives of conditional pdf wrt parameters(1st and 2nd order)
!cder1vec is 2d-vector while cder2vec is 3d-vector
!For one-parameter copula family, only the first element in cder1vec,cder2vec is present.
subroutine ccopderiv(u1,u2,cparv,icop,ccdf,cder1vec,cder2vec)
  implicit none
  integer icop
  double precision u1,u2,t1,t2,ut,vt,cparv(2),ccdf,cder1vec(2),cder2vec(3),cder11(2),cder22(3)
  double precision rho,nu,cpar,cder1,cder2
  double precision qt
  cpar=cparv(1)
  select case(icop)
    case (1) ! bivariate normal/Gaussian !check
    rho=cparv(1);nu=5000.0d0;
    t1=qt(u1,nu);t2=qt(u2,nu);
    call ctderivs(t1,t2,rho,nu,ccdf,cder1,cder2)
    cder1vec(1)=cder1;cder1vec(2)=0.0d0
    cder2vec(1)=cder2;cder2vec(2)=0.0d0;cder2vec(3)=0.0d0
    case (2) !student t !check
    rho=cparv(1);nu=cparv(2);
    t1=qt(u1,nu);t2=qt(u2,nu);
    call ctderivs(t1,t2,rho,nu,ccdf,cder1,cder2)
    cder1vec(1)=cder1;cder1vec(2)=0.0d0
    cder2vec(1)=cder2;cder2vec(2)=0.0d0;cder2vec(3)=0.0d0
    case (4) ! gumbel !check
    call cgumderivs(u1,u2,cpar,ccdf,cder1,cder2)
    cder1vec(1)=cder1;cder1vec(2)=0.0d0
    cder2vec(1)=cder2;cder2vec(2)=0.0d0;cder2vec(3)=0.0d0
    case (5) ! frank !check
    call clfrkderivs(u1,u2,cpar,ccdf,cder1,cder2)
    cder1vec(1)=cder1;cder1vec(2)=0.0d0
    cder2vec(1)=cder2;cder2vec(2)=0.0d0;cder2vec(3)=0.0d0
    case(7)!bb1 
    call cbb1derivs(u1,u2,cparv,ccdf,cder11,cder22)
    cder1vec=cder11
    cder2vec=cder22
    case(10)!bb8
    call cbb8derivs(u1,u2,cparv,ccdf,cder11,cder22)
    cder1vec=cder11
    cder2vec=cder22
    case(14) !survival gumbel !check
    ut=1.0d0-u1
    vt=1.0d0-u2
    call cgumderivs(ut,vt,cpar,ccdf,cder1,cder2)
    ccdf=1.0d0-ccdf
    cder1vec(1)=-cder1;cder1vec(2)=0.0d0
    cder2vec(1)=-cder2;cder2vec(2)=0.0d0;cder2vec(3)=0.0d0

    case(17)!survival bb1
    ut=1.0d0-u1
    vt=1.0d0-u2
    call cbb1derivs(ut,vt,cparv,ccdf,cder11,cder22)
    ccdf=1.0d0-ccdf
    cder1vec=-cder11
    cder2vec=-cder22
    
    case(20)!survival bb8
    ut=1.0d0-u1
    vt=1.0d0-u2
    call cbb8derivs(ut,vt,cparv,ccdf,cder11,cder22)
    ccdf=1.0d0-ccdf
    cder1vec=-cder11
    cder2vec=-cder22
  end select
  return
end



! function with t condcdf T_{1|2}(t1|t2;rho,nu+1) for factor 1, and 
! ccdf cder1, cder2 (partial wrt rho, 1st and 2nd order)
! inputs
!   t1 = qt(u1,nu) transformed to t scale
!   t2 = qt(u2,nu) transformed to t scale
!   rho = scalar in (-1,1)
!   nu = degree of freedom >0
! outputs 
!   ccdf = conditional cdf, 
!   cder1 = \p ccdf/\p rho
!   cder2 = \p^2 ccdf/\p rho^2
subroutine ctderivs(t1,t2,rho,nu,ccdf,cder1,cder2)
  implicit none
  double precision t1,t2,rho,rho2,nu,ccdf,cder1,cder2
  double precision logpi,r2,cr,xt,tem1,tem2,const,gtem1,gtem2
  double precision pt   ! function in C code ptbeta.c
  double precision lgamma   ! function in C code lgamma.c, in gnu fortran
  logpi=1.1447298858494d0
  rho2 = rho*rho;
  r2 = sqrt((nu+t2*t2)/(nu+1.d0));
  cr = sqrt(1.d0-rho2)
  xt = (t1-rho*t2)/cr/r2;
  ccdf = pt(xt,nu+1.d0); 
  tem1 = (t1*rho-t2)/cr**3/r2;
  tem2 = xt*xt/(nu+1.d0);
  const = exp(lgamma(0.5d0*nu+1.d0)-lgamma(0.5d0*nu+0.5d0)-0.5d0*logpi-0.5d0*log(nu+1.d0))
  cder1 = const*(1.d0+tem2)**(-0.5d0*nu-1.d0)*tem1;
  gtem1 = (t1+2.d0*t1*rho2-3.d0*t2*rho)/cr**5/r2;
  gtem2 = 2.d0*(rho*(t1*t1+t2*t2)-t1*t2*(rho2+1.d0))/(nu+1.d0)/cr**4/r2**2;
  cder2 = -const*(0.5d0*nu+1.d0)*(1.d0+tem2)**(-0.5d0*nu-2.d0)*gtem2*tem1;
  cder2 = cder2+const*(1.d0+tem2)**(-0.5d0*nu-1.d0)*gtem1;
  return
  end

! Function with Gumbel condcdf C_{1|2}(u1|u2;cpar) for factor 1, and 
! ccdf cder1, cder2 (partial wrt cpar, 1st and 2nd order)
! inputs
!   u1,u2 = values in (0,1)
!   cpar = scalar in (1,oo)
! outputs
!   ccdf = conditional cdf, 
!   cder1 = \p ccdf/\p cpar
!   cder2 = \p^2 ccdf/\p cpar^2
subroutine cgumderivs(u1,u2,cpar,ccdf,cder1,cder2)
  implicit none
  double precision u1,u2,cpar,ccdf,cder1,cder2
  double precision x,y,tx,ty,xd,yd,s,m,logs,logm,msq,dlsq,dlcu
  double precision sder1,sder2,mder1,mder2
  double precision lccdf,lcder1,lcder2
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
  logm=log(m); msq=m*m
  lccdf = y-m+(1.d0-cpar)*(logm-ty)
  !lcder1 = -mder1+(1.d0-cpar)*mder1/m-logm+tx;
  lcder1 = -mder1+(1.d0-cpar)*mder1/m-logm+ty;  ! given y
  lcder2 = -mder2-2.d0*mder1/m+(1.d0-cpar)*(mder2/m-mder1*mder1/msq);
  ccdf = exp(lccdf)
  cder1 = ccdf*lcder1
  cder2 = ccdf*(lcder1*lcder1+lcder2)
  return
  end



! Function with Frank condcdf C_{1|2}(u1|u2;cpar) for factor 1, and
! ccdf cder1, cder2 (partial wrt cpar, 1st and 2nd order)
! inputs
!   u1,u2 = values in (0,1)
!   cpar = scalar in (-oo,oo)
! outputs
!   ccdf = conditional cdf, 
!   cder1 = \p ccdf/\p cpar
!   cder2 = \p^2 ccdf/\p cpar^2
subroutine clfrkderivs(u1,u2,cpar,ccdf,cder1,cder2)
  implicit none
  double precision u1,u2,cpar,ccdf,cder1,cder2
  double precision t0,t1,t2,den,den1,den2,lccdf,lcder1,lcder2
  
  t0 = exp(-cpar);
  t1 = exp(-cpar*u1);
  t2 = exp(-cpar*u2);
  den = t1+t2-t0-t1*t2;
  den1 = -u1*t1-u2*t2+t0+(u1+u2)*t1*t2;
  den2 = u1*u1*t1+u2*u2*t2-t0-(u1+u2)*(u1+u2)*t1*t2;
  !lccdf = -cpar*u2+log(1-t1)-log(den);
  lccdf = -cpar*u2+log((1.d0-t1)/den);
  lcder1 = -u2+u1*t1/(1.d0-t1)-den1/den;
  lcder2 = -u1*u1*t1/((1.d0-t1)*(1.d0-t1))-den2/den+den1*den1/(den*den);
  ccdf = exp(lccdf)
  cder1 = ccdf*lcder1
  cder2 = ccdf*(lcder1*lcder1+lcder2)
  return
  end


! this routine is in bb1facts.f90
! Function with BB1 condcdf C_{1|2}(u1|u2;theta,delta) for factor 1, and
! ccdf cder11(2), cder22(3) (partial wrt param, 1st and 2nd order),
! second deriv has theta^2, mix, delta^2
! inputs
!   u1,u2 = values in (0,1)
!   param = (theta, delta), theta>0, delta>1
! outputs
!   ccdf = conditional cdf, 
!   cder11 = \p ccdf/\p param
!   cder22= \p^2 ccdf/\p param \p param^T
!
subroutine cbb1derivs(u1,u2,cparv,ccdf,cder11,cder22)
  implicit none
  double precision u1,u2,cparv(2),ccdf,cder11(2),cder22(3)
  double precision theta,delta
  double precision m,mder1th,mder1dl,mder2th,mder2dl,mderthd
  double precision t2,tu2,cf1,logm
  double precision mp1,msq,thsq,dl1n,t2a,lt2a,lmp1
  double precision lcdf,lder1th,lder1dl,lder2th,lder2dl,lderthd
  double precision cder1th,cder1dl,cder2th,cder2dl,cderthd

  theta=cparv(1); delta=cparv(2)

  call mderivs2(u1,u2,theta,delta,m,mder1th,mder1dl,mder2th,mder2dl,mderthd)
  mp1=1.d0+m; msq=m*m; thsq=theta*theta
  lmp1=log(mp1); logm=log(m);

  t2 = u2**(-theta); tu2 = -log(u2);
  t2a=t2-1.d0; lt2a=log(t2a)
  cf1 = 1.d0+1.d0/theta;
  dl1n=1.d0-delta

  lcdf = -cf1*lmp1+(dl1n)*(logm-lt2a)+(theta+1.d0)*tu2;
  lder1th = lmp1/thsq-cf1*mder1th/mp1;
  lder1th = lder1th+(dl1n)*(mder1th/m-t2*tu2/t2a)+tu2;
  lder1dl = -cf1*mder1dl/mp1+(dl1n)*mder1dl/m-logm+lt2a;
  lder2th = -2.d0*lmp1/(thsq*theta)+(2.d0/thsq)*mder1th/mp1-cf1*(mder2th/mp1&
  -mder1th**2/(mp1*mp1));
  lder2th = lder2th+(dl1n)*(mder2th/m-mder1th**2/msq+tu2*tu2*t2/(t2a*t2a));
  lder2dl = -cf1*(mder2dl/mp1-mder1dl**2/(mp1*mp1))-2.d0*mder1dl/m;
  lder2dl = lder2dl+(dl1n)*(mder2dl/m-mder1dl**2/msq);
  lderthd = (1.d0/thsq)*mder1dl/mp1-cf1*(mderthd/mp1-mder1th*mder1dl/(mp1*mp1))-mder1th/m;
  lderthd = lderthd+(dl1n)*(mderthd/m-mder1th*mder1dl/msq)+t2*tu2/t2a;

  ccdf = exp(lcdf);
  cder1th=ccdf*lder1th;
  cder1dl=ccdf*lder1dl;
  cder2dl=ccdf*(lder2dl+lder1dl**2);
  cder2th=ccdf*(lder2th+lder1th**2);
  cderthd=ccdf*(lderthd+lder1dl*lder1th);
  cder11(1)=cder1th; cder11(2)=cder1dl  ! gradient of log pdf
  cder22(1)=cder2th; cder22(2)=cderthd; cder22(3)=cder2dl; ! Hessian terms
  return
end



! this routine is in bb1facts.f90
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
subroutine mderivs2(u1,u2,theta,delta,m,mder1th,mder1dl,mder2th,mder2dl,mderthd)
  implicit none
  double precision u1,u2,theta,delta,m,mder1th,mder1dl,mder2th,mder2dl,mderthd
  double precision t1,t2,tu1,tu2,ttu1,ttu2,td01,td02,td11,td12,td21,td22
  double precision s,sder1th,sder1dl,sder2th,sder2dl,sderthd,m1,ts,m1der1dl 
  double precision t1a,t2a,dlsq,dlcu

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
  return
end

! Function with BB8 condcdf C_{1|2}(u1|u2;theta,delta) for factor 1, and
! ccdf cder1(2), cder2(3) (partial wrt param, 1st and 2nd order),
! second deriv has theta^2, mix, delta^2
! inputs
!   u1,u2 = values in (0,1)
!   param = (nu, delta), theta\in(1,8), delta\in (0,1)
! outputs
!   ccdf = conditional cdf, 
!   cder11 = \p ccdf/\p param
!   cder22 = \p^2 ccdf/\p param \p param^T
!
subroutine cbb8derivs(u1,u2,cparv,ccdf,cder11,cder22)
 implicit none
  double precision u1,u2,cparv(2),ccdf,cder11(2),cder22(3)
  double precision nu,delta,x,y,tem,eta,eta1nu,eta1del
  double precision x1nu,x1del,y1nu,y1del,tem1nu,tem1del,mm
  double precision tem1,tem2,tem3,tem4,eta2nu,eta2del,eta2nudel
  double precision x2nu,y2nu,x2del,y2del,x2nudel,y2nudel
  double precision tem2nu,tem2del, tem2nudel
  double precision lpdf,lpdf1nu,lpdf1del,lpdf2del,lpdf2nu,lpdf2nudel
  double precision cder1nu,cder1del,cder2nu,cder2del,xy1mixdel,cder2nudel
  
  nu=cparv(1)
  delta=cparv(2)
  eta=1.0d0-(1.0d0-delta)**nu
  x=1.0d0-(1.0d0-delta*u1)**nu
  y=1.0d0-(1.0d0-delta*u2)**nu
  tem=x*y/eta
  mm=1.0d0/nu-1.0d0

  eta1nu=-log(1.0d0-delta)*((1.0d0-delta)**nu) !\p eta/\p nu
  eta1del=nu*((1.0d0-delta)**(nu-1.0d0))!\p eta/\p delta
  x1nu=-log(1.0d0-delta*u1)*((1-delta*u1)**nu) !\p x/\p nu
  x1del=(nu*(1.0d0-delta*u1)**(nu-1))*u1!\p x/\p delta
  y1nu=-log(1.0d0-delta*u2)*((1-delta*u2)**nu)!\p y/\p nu
  y1del=nu*((1.0d0-delta*u2)**(nu-1))*u2!\p y/\p delta
  tem1nu=-(x*y/eta**2)*eta1nu+(x1nu*y+y1nu*x)/eta  !\p tem/\p nu
  tem1del=(x1del*y+y1del*x)/eta-eta1del*x*y/eta**2!\p tem/\p delta
 
  eta2nu=-((log(1.0d0-delta))**2)*((1.0d0-delta)**nu) !\p^2 eta/\p nu^2
  eta2del=-nu*(nu-1.0d0)*(1.0d0-delta)**(nu-2.0d0) !\p^2 eta/\p delta^2
 !\p^2 eta/\p delta*nu
  eta2nudel=nu*log(1-delta)*((1.0d0-delta)**(nu-1.0d0))+(1.0d0-delta)**(nu-1.0d0)
  x2nu=-(log(1.0d0-delta*u1)**2)*((1.0d0-delta*u1)**nu)!\p^2 x/\p nu^2
  y2nu=-(log(1.0d0-delta*u2)**2)*((1.0d0-delta*u2)**nu)!\p^2 y/\p nu^2
  x2del=-u1*u1*nu*(nu-1.0d0)*((1-delta*u1)**(nu-2.0d0))!\p^2 x/\p delta^2
  y2del=-u2*u2*nu*(nu-1.0d0)*((1-delta*u2)**(nu-2.0d0))!\p^2 y/\p delta^2
 !\p^2 x/\p delta*nu
  x2nudel=u1*nu*(1-delta*u1)**(nu-1.0d0)*log(1-delta*u1)+u1*((1-delta*u1)**nu)/(1-delta*u1)
!\p^2 y/\p delta*nu
  y2nudel=u2*nu*(1-delta*u2)**(nu-1.0d0)*log(1-delta*u2)+u2*((1-delta*u2)**nu)/(1-delta*u2)
 
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
  tem2nudel=tem1+tem2+tem3+tem4 !\p^2 tem/\p nu*delta
 
  lpdf=-log(eta)+log(x)+(1.0d0/nu-1.0d0)*log(1.0d0-tem)+(1.0d0-1.0d0/nu)*log(1.0d0-y)
  ccdf=exp(lpdf)
 !ccdf=x*(1-tem)**(1.0d0/nu-1.0d0)*(1.0d0-y)**(1.0d0-1.0d0/nu)/eta
 
  tem1=-eta1nu/eta+x1nu/x-log(1.0d0-tem)/nu**2-tem1nu*mm/(1.0d0-tem)
  tem2=log(1.0d0-y)/nu**2+y1nu*mm/(1.0d0-y)
  lpdf1nu=tem1+tem2
  cder1nu=ccdf*lpdf1nu !\p ccdf/\p nu
 
  lpdf1del=-eta1del/eta+x1del/x+(1.0d0/nu-1.0d0)*(-tem1del)/(1.0d0-tem)&
  +(1.0d0-1.0d0/nu)*(-y1del)/(1.0d0-y)
  cder1del=ccdf*lpdf1del!\p ccdf/\p delta
 
  tem1=(eta1nu)**2/eta**2-eta2nu/eta-(x1nu)**2/x**2+x2nu/x
  tem2=2*log(1.0d0-tem)/nu**3+tem1nu/(1.0d0-tem)/nu**2
  tem3=-mm*tem1nu**2/(1.0d0-tem)**2-(tem2nu*mm-tem1nu/nu**2)/(1.0d0-tem)
  tem4=-2.0d0*log(1.0d0-y)/nu**3-y1nu/(1.0d0-y)/nu**2-(-y2nu*mm/(1.0d0-y)&
  +(1.0d0/(1.0d0-y)/nu**2-y1nu*mm/(1.0d0-y)**2)*y1nu)
  lpdf2nu=tem1+tem2+tem3+tem4 !\p^2 lpdf/\p nu^2
  cder2nu=ccdf*(lpdf1nu)**2+lpdf2nu*ccdf!\p^2 ccdf/\p nu^2
 
  tem1=eta1del**2/eta**2-eta2del/eta-x1del**2/x**2+x2del/x
  tem2=(-tem1del**2/(1.0d0-tem)**2-tem2del/(1.0d0-tem))*(1.0d0/nu-1.0d0)
  tem3=-mm*(-y2del/(1.0d0-y)-y1del**2/(1.0d0-y)**2)
  lpdf2del=tem1+tem2+tem3!\p^2 lpdf/\p delta^2
  cder2del=ccdf*(lpdf1del)**2+lpdf2del*ccdf!\p^2 ccdf/\p delta^2
 
  tem1=eta1nu*eta1del/eta**2-eta2nudel/eta
  tem2=-x1del*x1nu/x**2+x2nudel/x+tem1del/(1.0d0-tem)/nu**2
  tem3=-(mm*tem1nu*tem1del/(1.0d0-tem)**2+mm*tem2nudel/(1.0d0-tem))
  tem4=-y1del/(1.0d0-y)/nu**2+mm*(y2nudel/(1.0d0-y)+y1del*y1nu/(1.0d0-y)**2)
  lpdf2nudel=tem1+tem2+tem3+tem4!\p^2 lpdf/\p delta*nu
  cder2nudel=ccdf*lpdf1nu*lpdf1del+lpdf2nudel*ccdf!\p^2 ccdf/\p nu*delta
  cder11(1)=cder1nu;cder11(2)=cder1del
  cder22(1)=cder2nu;cder22(2)=cder2nudel;cder22(3)=cder2del
return
end