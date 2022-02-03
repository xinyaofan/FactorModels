! PROGRAM test   ! main function
! implicit none
! 	integer npar,mgrp,dvar,n,i,nq,iq
! 	double precision v0,int11,int22
! 	double precision,dimension(:,:),allocatable::udata
! 	double precision,dimension(:),allocatable::th
! 	integer,dimension(:),allocatable::grsize
! 	double precision,dimension(:),allocatable::xl,wl,v0mat
! 	double precision,dimension(:,:),allocatable::vgmat
! 	
! 	read *,nq
!  	allocate(xl(nq),wl(nq))
!  	do iq =1,nq   ! equidistant and equally weighted for testing
!           xl(iq)=iq/(nq+1.d0)
!           wl(iq)=1.d0/nq
!   	end do
!     read *,npar
!     read *,dvar
!     read *,mgrp
!     allocate(grsize(mgrp))
!     read *,grsize
!     read*,n
!     allocate(th(npar))
!  	read *,th
!  	allocate(udata(n,dvar))
!  	read *, (udata(i,:), i=1,n) 
!     allocate(v0mat(n),vgmat(n,mgrp))
!  	call latupdatebifact (npar,th,mgrp,dvar,n,grsize,udata,nq,xl,wl,v0mat,vgmat)
!  	print *, "the estimats for global latent is "
!  	print *,v0mat
!  	print *,"the estimates for local latent is"
!  	do i=1,n
!         print *, vgmat(i,:)
!     end do
!  	deallocate(grsize,th,udata,xl,wl,v0mat,vgmat)
! END PROGRAM

! program test
!   implicit none
!   integer family,nq,n,d,i
!   double precision,dimension(:,:),allocatable::udata
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
!   allocate(udata(n,d),th(d),lat(n))
!   read *, (udata(i,:), i=1,n) 
!   
!   read*,th
!   print*,th
!   read*,family
!   print*,family
!   
!   call latupdate (th,n,d,udata,nq,xl,wl,family,lat)
!   
!   print*,lat
!   
!   deallocate(udata,xl,wl)
! end program 

! subroutine innerint(npar,th,mgrp,dvar,grsize,uvec,nq,wl,xl,v0,vg,int11,int22)
! 	implicit none
! 	integer npar,mgrp,dvar,n,nq,ind,i
! 	double precision th(npar),wl(nq),xl(nq)
!   	double precision pdf(npar),ccdf(dvar),uvec(dvar)
!   	double precision nllk,llk,lk1,lk2,v0,vg(mgrp),int11,int22
!   	integer iq2,jg,mj,ind1,ind2,grsize(mgrp)
!   	double precision intj(mgrp), intj2(mgrp),int1(mgrp),int2(mgrp),int0
! 	
! 	 int0 = 1.d0; ind = 0; intj = 0.d0; intj2=0.d0; int1=0.d0; int2=0.d0;
! 	  do jg =1,mgrp   !jth group                         
!         ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg);      
!         do iq2 =1,nq   
!           lk1 = 1.d0; lk2 = 1.d0   
!           do mj = ind1,ind2 ! within group index
!             call pcondfrk(uvec(mj),v0,th(mj),ccdf(mj)) !C_{ij|V_0}
!             if(ccdf(mj) < 0.00001) ccdf(mj) = 0.00001
!             if(ccdf(mj) > 0.99999) ccdf(mj) = 0.99999    
!             call dfrk(uvec(mj),v0,th(mj),pdf(mj)) !c_{ij,V_0}  
!             call dfrk(ccdf(mj),xl(iq2),th(dvar+mj),pdf(dvar+mj)) !c_{i,V_j;V_0}  
!             lk1 = lk1*pdf(mj+dvar)
!             lk2 = lk2*pdf(mj)
!           end do
!           int1(jg)=int1(jg)+wl(iq2)*lk1
!           int2(jg)=lk2
!           llk=lk1*lk2
!           intj(jg)=intj(jg)+wl(iq2)*llk*xl(iq2)
!           intj2(jg)=intj2(jg)+wl(iq2)*llk
!           
!         end do   
!       end do
!       do i=1,mgrp
!       	vg(i)=intj(i)/intj2(i)
!       end do
!       int11=product(int1)
!       int22=product(int2)
!    return
!    end 
!    
!    
! 
! subroutine latupdatebifact (npar,th,mgrp,dvar,n,grsize,udata,nq,xl,wl,v0mat,vgmat)
! 	implicit none
! 	integer npar,mgrp,n,nq,ind,iq,i,dvar
! 	integer int1vec(nq),int2vec(nq),grsize(mgrp)
! 	double precision th(npar),udata(n,dvar),wl(nq),xl(nq),uvec(dvar)
! 	double precision v0,vg(mgrp),int11,int22
! 	double precision num,dem,tem(nq),tem1(nq),v0mat(n),vgmat(n,mgrp)
! 	do i=1,n
! 	uvec = udata(i,:)   
! 		int1vec=0.d0;int2vec=0.d0;
! 		do iq =1,nq
! 		v0=xl(iq)
! 		call innerint(npar,th,mgrp,dvar,grsize,uvec,nq,wl,xl,v0,vg,int11,int22)
! 		int1vec(iq)=int11;int2vec(iq)=int22
! 		tem1(iq)=int11*int22
! 		tem(iq)=xl(iq)*wl(iq)
!      	end do
!      num=dot_product(tem,tem1)
! 	 dem=dot_product(tem1,wl)
! 	 v0mat(i)=num/dem
! 	 v0=v0mat(i)
! 	 call innerint(npar,th,mgrp,dvar,grsize,uvec,nq,wl,xl,v0,vg,int11,int22)
! 	 vgmat(i,:)=vg
! 	end do
! return
! end subroutine

!latent update in one-factor copula model; only the same family are allowed
subroutine latupdate (th,n,d,udata,nq,xl,wl,family,lat)
  implicit none
  integer n,nq,i,iq,j,d
  integer family
  double precision num,dem,lpdf
  double precision udata(n,d),th(d),wl(nq),xl(nq),uvec(d),val(nq),tem(nq),lat(n)
  do i=1,n
    uvec=udata(i,:)
    val=0.0d0;
    do j=1,d
     do iq=1,nq
     if (family==4) then
     call lgum(uvec(j),xl(iq),th(j),lpdf)
     elseif(family==5) then
      call lfrk(uvec(j),xl(iq),th(j),lpdf)
      elseif(family==2) then
      call lgau(uvec(j),xl(iq),th(j),lpdf)
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
subroutine lgum(u1,u2,cpar,lpdf)
  implicit none
  double precision u1,u2,cpar,lpdf
  double precision x,y,tx,ty,logm,den,m,s,xd,yd
  
  x = -log(u1); y = -log(u2);
  tx = log(x); ty = log(y);
  xd = x**cpar; yd = y**cpar; s = xd+yd; m = s**(1.d0/cpar); 
  den = m+cpar-1.d0
  
  logm=log(m);
  lpdf = -m+log(den)+(1.d0-2.d0*cpar)*logm+(cpar-1.d0)*(tx+ty)+x+y;
  
end subroutine


!This functions returns the log-density of Frank copula
!inputs: u1,u2,cpar:paramter in the copula 
!output: lpdf: log density of Frank copula 
subroutine lfrk(u1,u2,cpar,lpdf)
  implicit none
  double precision u1,u2,cpar,lpdf
  double precision den,t0,t1,t2

  t0 = exp(-cpar);
  t1 = exp(-cpar*u1);
  t2 = exp(-cpar*u2);
  den = t1+t2-t0-t1*t2;
  lpdf = log(cpar*(1.d0-t0))-cpar*(u1+u2)-2.d0*log(abs(den));
  
end subroutine

subroutine lgau(u1,u2,rho,lpdf)
  implicit none
  double precision u1,u2,rho,lpdf
  double precision qnorms
  double precision x,y,x2,y2,rh2,con,den,den2,qf1,qf2

  x=qnorms(u1); x2=x*x
  y=qnorms(u2); y2=y*y
  rh2=rho**2; den=1.d0-rh2; den2=den*den
  con= -0.5d0*log(den)
  qf1= -(x2+y2-2.d0*rho*x*y)/(2*den)
  qf2=(x2+y2)/2.d0
  lpdf=con+qf1+qf2
  return
  end
!compute the density of frank copula	
! subroutine dfrk(x1,x2,cpar,pdf)
! 	implicit none
! 	double precision:: x1,x2
! 	double precision :: cpar,tem1,tem2,tem,t1,tem3,pdf
! 	
! 	t1=1.d0-exp(-cpar);
! 	tem1=exp(-cpar*x1); tem2=exp(-cpar*x2);
! 	tem=t1-(1.d0-tem1)*(1.d0-tem2);
! 	tem3=cpar*tem1*tem2*t1
! 	pdf=tem3/(tem*tem);
! 	return
! end subroutine dfrk
!compute the conditional pdf of frank copula

! subroutine pcondfrk(x2,x1,cpar,ccdf)
! 	implicit none
! 	double precision:: x2,x1
! 	double precision :: cpar,cpar1,tem,ccdf
! 	
! 	cpar1=1.d0-exp(-cpar)  ! 1 in double precision is 1.d0 in fortran
! 	tem=1.d0-exp(-cpar*x1) ! 1.d0
! 	ccdf=(1.d0-tem)/(cpar1/(1.d0-exp(-cpar*x2))-tem) ! 1.d0
! end subroutine pcondfrk
