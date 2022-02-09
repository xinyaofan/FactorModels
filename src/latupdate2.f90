
!main function 
! PROGRAM test
!   implicit none
!   integer nq,n,d,i
!   double precision,dimension(:,:),allocatable::udata
!   integer,dimension(:),allocatable::family
!   double precision,dimension(:),allocatable::xl,wl,lat,th
! 
!   read*,nq
!   print*,nq
!   allocate(xl(nq),wl(nq))
!   read*,xl
!   read*,wl
!   
!   print*,xl
!   print*,wl
!   
!   read*,n
!   read*,d
!   allocate(udata(n,d),th(2*d),lat(n),family(d))
!   read *, (udata(i,:), i=1,n) 
!   
!   read*,th
!   print*,th
!   read*,family
!   print*,family
!   
!   call latupdate2(th,n,d,udata,nq,xl,wl,family,lat)
!   
!   print*,lat
!   
!   deallocate(udata,xl,wl,family)
! end PROGRAM



!latent update in one-factor copula model; different copula families are allowed

subroutine latupdate3 (th,n,d,udata,nq,xl,wl,family,lat)
  implicit none
  integer n,nq,i,iq,j,d,jj(2),jj2
  double precision num,dem,lpdf,rho,nu,t1,t2
  integer family(d)
  double precision lder11,lder22,lderu,lderv,lderuv,lder1,lder2,ltder1u,ltder2u,ltdermix
  double precision udata(n,d),th(2*d),wl(nq),xl(nq),uvec(d),val(nq),tem(nq),lat(n)
  double precision ut,vt,qt
  do i=1,n
    uvec=udata(i,:)
    val=0.0d0;
    do j=1,d
      jj=(/2*j-1,2*j/)
      jj2=2*j-1
     do iq=1,nq
      if (family(j)==4) then
         call lgum(uvec(j),xl(iq),th(jj2),lpdf)
      elseif(family(j)==5) then
         call lfrk(uvec(j),xl(iq),th(jj2),lpdf)
      elseif(family(j)==1)then
         call lgau(uvec(j),xl(iq),th(jj2),lpdf)
      elseif(family(j)==2) then
        rho=th(2*j-1);nu=th(2*j);
        t1=qt(uvec(j),nu);t2=qt(xl(iq),nu);
        call lt2derivs(t1,t2,rho,nu,lpdf,lder1,lder2,ltder1u,ltder2u,ltdermix)
      elseif(family(j)==7) then
        call lbb1derivs(uvec(j),xl(iq),th(jj),lpdf,lder11,lder22,lderu,lderv,lderuv)
      elseif(family(j)==17) then
        ut=1.0d0-uvec(j)
        vt=1.0d0-xl(iq)
        call lbb1derivs(ut,vt,th(jj),lpdf,lder11,lder22,lderu,lderv,lderuv)
      elseif(family(j)==10) then
        call lbb8derivs(uvec(j),xl(iq),th(jj),lpdf,lder11,lder22,lderu,lderv,lderuv)
      elseif(family(j)==20) then
        ut=1.0d0-uvec(j)
        vt=1.0d0-xl(iq)
        call lbb8derivs(ut,vt,th(jj),lpdf,lder11,lder22,lderu,lderv,lderuv)
     end if
      val(iq)=val(iq)+lpdf
      end do
    end do
    tem=exp(val)
    num=dot_product(wl,tem*xl)
    dem=dot_product(wl,tem)
    lat(i)=num/dem  
  end do   
return
end subroutine





!this functions returns the log-density of Gumbel copula
!inputs: u1,u2,cpar:paramter in the copula 
!output: lpdf: log density of Gumbel copula 
! subroutine lgum(u1,u2,cpar,lpdf)
!   implicit none
!   double precision u1,u2,cpar,lpdf
!   double precision x,y,tx,ty,logm,den,m,s,xd,yd
!   
!   x = -log(u1); y = -log(u2);
!   tx = log(x); ty = log(y);
!   xd = x**cpar; yd = y**cpar; s = xd+yd; m = s**(1.d0/cpar); 
!   den = m+cpar-1.d0
!   
!   logm=log(m);
!   lpdf = -m+log(den)+(1.d0-2.d0*cpar)*logm+(cpar-1.d0)*(tx+ty)+x+y;
!   
! end subroutine


!This functions returns the log-density of Frank copula
!inputs: u1,u2,cpar:paramter in the copula 
!output: lpdf: log density of Frank copula 
! subroutine lfrk(u1,u2,cpar,lpdf)
!   implicit none
!   double precision u1,u2,cpar,lpdf
!   double precision den,t0,t1,t2
! 
!   t0 = exp(-cpar);
!   t1 = exp(-cpar*u1);
!   t2 = exp(-cpar*u2);
!   den = t1+t2-t0-t1*t2;
!   lpdf = log(cpar*(1.d0-t0))-cpar*(u1+u2)-2.d0*log(abs(den));
!   
! end subroutine
! 
! subroutine lgau(u1,u2,rho,lpdf)
!   implicit none
!   double precision u1,u2,rho,lpdf
!   double precision qnorms
!   double precision x,y,x2,y2,rh2,con,den,den2,qf1,qf2
! 
!   x=qnorms(u1); x2=x*x
!   y=qnorms(u2); y2=y*y
!   rh2=rho**2; den=1.d0-rh2; den2=den*den
!   con= -0.5d0*log(den)
!   qf1= -(x2+y2-2.d0*rho*x*y)/(2*den)
!   qf2=(x2+y2)/2.d0
!   lpdf=con+qf1+qf2
!   return
!   end