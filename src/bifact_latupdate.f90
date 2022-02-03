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

 !PROGRAM test   ! main function
 !	implicit none
 !   double precision::x1,x2
 !   double precision::pdf,ccdf,cpar
 !   read*,x1
 !   read*,x2
  !  read*,cpar
  !  call pcondgum(x2,x1,cpar,ccdf)
 !   print*,ccdf
!END program
!compute the inner integral for frank copula 
 subroutine innerint(npar,th,mgrp,dvar,grsize,uvec,nq,wl,xl,v0,vg,int11,int22)
  implicit none
  integer npar,mgrp,dvar,nq,ind,i
  double precision th(npar),wl(nq),xl(nq)
  double precision pdf(npar),ccdf(dvar),uvec(dvar)
  double precision llk,lk1,lk2,v0,vg(mgrp),int11,int22
  integer iq2,jg,mj,ind1,ind2,grsize(mgrp)
  double precision intj(mgrp), intj2(mgrp),int1(mgrp),int2(mgrp),int0

   int0 = 1.d0; ind = 0; intj = 0.d0; intj2=0.d0; int1=0.d0; int2=0.d0;
   do jg =1,mgrp   !jth group                         
         ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg);      
         do iq2 =1,nq   
           lk1 = 1.d0; lk2 = 1.d0   
           do mj = ind1,ind2 ! within group index
             call pcondfrk(uvec(mj),v0,th(mj),ccdf(mj)) !C_{ij|V_0}
             if(ccdf(mj) < 0.00001) ccdf(mj) = 0.00001
             if(ccdf(mj) > 0.99999) ccdf(mj) = 0.99999    
             call dfrk(uvec(mj),v0,th(mj),pdf(mj)) !c_{ij,V_0}  
             call dfrk(ccdf(mj),xl(iq2),th(dvar+mj),pdf(dvar+mj)) !c_{i,V_j;V_0}  
             lk1 = lk1*pdf(mj+dvar)
             lk2 = lk2*pdf(mj)
           end do
           int1(jg)=int1(jg)+wl(iq2)*lk1
           int2(jg)=lk2
           llk=lk1*lk2
           intj(jg)=intj(jg)+wl(iq2)*llk*xl(iq2)
           intj2(jg)=intj2(jg)+wl(iq2)*llk
           
         end do   
       end do
       do i=1,mgrp
        vg(i)=intj(i)/intj2(i)
       end do
       int11=product(int1)
       int22=product(int2)
    return
   end 
    
!compute the inner integral for Gumbel copula 
 subroutine innerint2(npar,th,mgrp,dvar,grsize,uvec,nq,wl,xl,v0,vg,int11,int22)
  implicit none
  integer npar,mgrp,dvar,nq,i,ind
  double precision th(npar),wl(nq),xl(nq)
  double precision pdf(npar),ccdf(dvar),uvec(dvar)
  double precision llk,lk1,lk2,v0,vg(mgrp),int11,int22
  integer iq2,jg,mj,ind1,ind2,grsize(mgrp)
  double precision intj(mgrp), intj2(mgrp),int1(mgrp),int2(mgrp),int0
 
  int0 = 1.d0; ind = 0; intj = 0.d0; intj2=0.d0; int1=0.d0; int2=0.d0;
    do jg =1,mgrp   !jth group                         
         ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg);      
         do iq2 =1,nq   
           lk1 = 1.d0; lk2 = 1.d0   
           do mj = ind1,ind2 ! within group index
             call pcondgum(uvec(mj),v0,th(mj),ccdf(mj)) !C_{ij|V_0}
             if(ccdf(mj) < 0.00001) ccdf(mj) = 0.00001
             if(ccdf(mj) > 0.99999) ccdf(mj) = 0.99999    
             call dgum(uvec(mj),v0,th(mj),pdf(mj)) !c_{ij,V_0}  
             call dgum(ccdf(mj),xl(iq2),th(dvar+mj),pdf(dvar+mj)) !c_{i,V_j;V_0}  
             lk1 = lk1*pdf(mj+dvar)
             lk2 = lk2*pdf(mj)
           end do
           int1(jg)=int1(jg)+wl(iq2)*lk1
           int2(jg)=lk2
           llk=lk1*lk2
           intj(jg)=intj(jg)+wl(iq2)*llk*xl(iq2)
           intj2(jg)=intj2(jg)+wl(iq2)*llk
           
         end do   
       end do
       do i=1,mgrp
        vg(i)=intj(i)/intj2(i)
       end do
       int11=product(int1)
       int22=product(int2)
    return
   end  
! 
subroutine latupdatebifact (npar,th,mgrp,dvar,n,grsize,udata,nq,xl,wl,v0mat,vgmat)
  implicit none
  integer npar,mgrp,n,nq,iq,i,dvar
  integer grsize(mgrp)
  double precision th(npar),udata(n,dvar),wl(nq),xl(nq),uvec(dvar)
  double precision v0,vg(mgrp),int11,int22,int1vec(nq),int2vec(nq)
  double precision num,dem,tem(nq),tem1(nq),v0mat(n),vgmat(n,mgrp)
  do i=1,n
     uvec = udata(i,:)   
     int1vec=0.d0;int2vec=0.d0;
     do iq =1,nq
      v0=xl(iq)
      call innerint(npar,th,mgrp,dvar,grsize,uvec,nq,wl,xl,v0,vg,int11,int22)
      int1vec(iq)=int11;int2vec(iq)=int22
      tem1(iq)=int11*int22
      tem(iq)=xl(iq)*wl(iq)
      end do
      num=dot_product(tem,tem1)
      dem=dot_product(tem1,wl)
      v0mat(i)=num/dem
      v0=v0mat(i)
      call innerint(npar,th,mgrp,dvar,grsize,uvec,nq,wl,xl,v0,vg,int11,int22)
      vgmat(i,:)=vg
    end do
 return
end subroutine


subroutine latupdatebifact2 (npar,th,mgrp,dvar,n,grsize,udata,nq,xl,wl,v0mat,vgmat)
  implicit none
  integer npar,mgrp,n,nq,iq,i,dvar
  integer grsize(mgrp)
  double precision th(npar),udata(n,dvar),wl(nq),xl(nq),uvec(dvar)
  double precision v0,vg(mgrp),int11,int22,int1vec(nq),int2vec(nq)
  double precision num,dem,tem(nq),tem1(nq),v0mat(n),vgmat(n,mgrp)
  do i=1,n
     uvec = udata(i,:)   
     int1vec=0.d0;int2vec=0.d0;
     do iq =1,nq
      v0=xl(iq)
      call innerint2(npar,th,mgrp,dvar,grsize,uvec,nq,wl,xl,v0,vg,int11,int22)
      int1vec(iq)=int11;int2vec(iq)=int22
      tem1(iq)=int11*int22
      tem(iq)=xl(iq)*wl(iq)
      end do
     num=dot_product(tem,tem1)
     dem=dot_product(tem1,wl)
     v0mat(i)=num/dem
     v0=v0mat(i)
     call innerint2(npar,th,mgrp,dvar,grsize,uvec,nq,wl,xl,v0,vg,int11,int22)
     vgmat(i,:)=vg
   end do
 return
end subroutine
! subroutine latupdate (npar,th,n,d,udata,nq,xl,wl,family,lat)
!   implicit none
!   integer npar,n,nq,i,iq
!   integer family(d),th(d)
!   double precision num,dem
!   double precision udata(n,d),wl(nq),xl(nq),uvec(d),val(nq),tem(nq),lat(n)
!   do i=1,n
!     uvec=udata(i,:)
!     do j=1,d
!     if (fam==4) !gumbel
!       do iq=1,nq
!       call lgum(uvec(j),xl(iq),cpar=th(j),lpdf)
!       val(iq)=val(iq)+lpdf
!     else
!      do iq=1,nq
!       call lfrk(uvec(j),xl(iq),cpar=th(j),lpdf)
!       val(iq)=val(iq)+lpdf
!       end do
!     end do
!     tem=exp(val)
!     num=dot_product(wl,dot_product(tem,xl))
!     dem=dot_product(wl,tem)
!     lat(i)=num/dem  
!   end do   
! return
! end subroutine


!compute the density of frank copula	
!  subroutine dfrk(x1,x2,cpar,pdf)
!  	implicit none
!  	double precision:: x1,x2
!  	double precision :: cpar,tem1,tem2,tem,t1,tem3,pdf
 
!  	t1=1.d0-exp(-cpar);
!  	tem1=exp(-cpar*x1); tem2=exp(-cpar*x2);
!  	tem=t1-(1.d0-tem1)*(1.d0-tem2);
!  	tem3=cpar*tem1*tem2*t1
!  	pdf=tem3/(tem*tem);
!  	return
!  end subroutine dfrk
! !compute the conditional pdf of frank copula

!  subroutine pcondfrk(x2,x1,cpar,ccdf)
!  	implicit none
!  	double precision:: x2,x1
!  	double precision :: cpar,cpar1,tem,ccdf
 
!  	cpar1=1.d0-exp(-cpar)  ! 1 in double precision is 1.d0 in fortran
!  	tem=1.d0-exp(-cpar*x1) ! 1.d0
!  	ccdf=(1.d0-tem)/(cpar1/(1.d0-exp(-cpar*x2))-tem) ! 1.d0
!  end subroutine pcondfrk


! !compute the density of Gumbel copula
! subroutine dgum(x1,x2,cpar,pdf)
! 	implicit none
! 	double precision::x1,x2
! 	double precision::cpar,pdf,l1,l2,tem1,tem2,sm,tem,cdf
! 	l1=-log(x1)
! 	l2=-log(x2)
! 	tem1=l1**cpar
! 	tem2=l2**cpar
! 	sm=tem1+tem2
! 	tem=sm**(1.0d0/cpar)
! 	cdf=exp(-tem)
! 	pdf=cdf*tem*tem1*tem2*(tem+cpar-1.0d0)
! 	pdf=pdf/(sm*sm*l1*l2*x1*x2)
! 	return
! end subroutine dgum

! subroutine pcondgum(x2,x1,cpar,ccdf)
! 	implicit none
! 	double precision::x2,x1,x,y
! 	double precision::cpar,tem1,tem2,sum,tem,ccdf
! 	x=-log(x1)
! 	y=-log(x2)
! 	tem1=x**cpar
! 	tem2=y**cpar
! 	sum=tem1+tem2
! 	tem=sum**(1.0d0/cpar)
! 	ccdf=exp(-tem)
! 	ccdf=ccdf*(1.0d0+tem2/tem1)**(-1.0d0+1.0d0/cpar)
! 	ccdf=ccdf/x1	
! end subroutine pcondgum