! PROGRAM test   ! main function
! implicit none
! 	integer npar,mgrp,dvar,n,i
! 	double precision,dimension(:,:),allocatable::udata,vlat
! 	double precision,dimension(:),allocatable::th
! 	integer,dimension(:),allocatable::grsize
! 	double precision nllk
!     read *,npar
!     read *,dvar
!     read *,mgrp
!     allocate(grsize(mgrp))
!     read *,grsize
!     read*,n
!     allocate(th(npar),udata(n,dvar),vlat(n,mgrp+1))
!  	read *,th
!  	read *, (udata(i,:), i=1,n) 
!  	! do i=1,n
!     !      print *, udata(i,:)
!     !    end do
!  	read *, (vlat(i,:), i=1,n) 
!  	! do i=1,n
!      !     print *, vlat(i,:)
!      !   end do
!     
!  	call frkproxynllk(npar,th,mgrp,n,dvar,grsize,udata,vlat,nllk)
!  	
!  	print *, "the negative log-likelihood is"
!  	print *,nllk
!  	deallocate(grsize,udata,vlat,th)
! END PROGRAM


subroutine frkproxynllk(npar,th,mgrp,n,dvar,grsize,udata,vlat,nllk)
  implicit none
  integer npar,mgrp,dvar,n
  double precision th(npar),udata(n,dvar),uvec(dvar),vlat(n,mgrp+1)
  integer grsize(mgrp)
  double precision nllk,lk1,lk2,liki,ccdf(dvar),pdf(npar)
  integer i,jg,mj,ind1,ind2,ind,int0
  double precision intj(mgrp)

   nllk=0.d0; 
   do i =1,n 
       uvec = udata(i,:) 
       int0 = 1.d0; ind = 0.d0; intj = 1.d0;  
      do jg =1,mgrp   !jth group                         
        ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg);lk1 = 1.d0; lk2 = 1.d0;  
          do mj = ind1,ind2 ! within group index
            call pcondfrk(uvec(mj),vlat(i,1),th(mj),ccdf(mj)) !C_{ij|V_0}
            call dfrk(uvec(mj),vlat(i,1),th(mj),pdf(mj)) !c_{ij,V_0} 
            call dfrk(ccdf(mj),vlat(i,jg+1),th(dvar+mj),pdf(dvar+mj)) !c_{i,V_j;V_0}  
    
            lk1 = lk1*pdf(mj)
            lk2 = lk2*pdf(dvar+mj)
          end do
          intj(jg) = intj(jg)*lk1*lk2 !intj value of the jth integral
        end do   
      liki = product(intj)   !product of j inner integrals  
      nllk = nllk - log(liki)       !updating loglikelihood
  end do
  return
  end

subroutine gumproxynllk(npar,th,mgrp,n,dvar,grsize,udata,vlat,nllk)
  implicit none
  integer npar,mgrp,dvar,n
  double precision th(npar),udata(n,dvar),uvec(dvar),vlat(n,mgrp+1)
  integer grsize(mgrp)
  double precision nllk,lk1,lk2,liki,ccdf(dvar),pdf(npar)
  integer i,jg,mj,ind1,ind2,ind,int0
  double precision intj(mgrp)

  nllk=0.d0; 
  do i =1,n 
       uvec = udata(i,:) 
       int0 = 1.d0; ind = 0.d0; intj = 1.d0;  
      do jg =1,mgrp   !jth group                         
        ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg);lk1 = 1.d0; lk2 = 1.d0;  
          do mj = ind1,ind2 ! within group index
            call pcondgum(uvec(mj),vlat(i,1),th(mj),ccdf(mj)) !C_{ij|V_0}
            call dgum(uvec(mj),vlat(i,1),th(mj),pdf(mj)) !c_{ij,V_0} 
            call dgum(ccdf(mj),vlat(i,jg+1),th(dvar+mj),pdf(dvar+mj)) !c_{i,V_j;V_0}  
            lk1 = lk1*pdf(mj)
            lk2 = lk2*pdf(dvar+mj)
          end do
          intj(jg) = intj(jg)*lk1*lk2 !intj value of the jth integral
        end do   
      liki = product(intj)   !product of j inner integrals  
      nllk = nllk - log(liki)       !updating loglikelihood
  end do
  return
  end
  
  
!compute the density of Gumbel copula
subroutine dgum(x1,x2,cpar,pdf)
  implicit none
  double precision::x1,x2
  double precision::cpar,pdf,l1,l2,tem1,tem2,sm,tem,cdf
  l1=-log(x1)
  l2=-log(x2)
  tem1=l1**cpar
  tem2=l2**cpar
  sm=tem1+tem2
  tem=sm**(1.0d0/cpar)
  cdf=exp(-tem)
  pdf=cdf*tem*tem1*tem2*(tem+cpar-1.0d0)
  pdf=pdf/(sm*sm*l1*l2*x1*x2)
  return
end subroutine dgum
  
!compute the density of frank copula	
subroutine dfrk(x1,x2,cpar,pdf)
  implicit none
  double precision:: x1,x2
  double precision :: cpar,tem1,tem2,tem,t1,tem3,pdf
  t1=1.d0-exp(-cpar);
  tem1=exp(-cpar*x1); tem2=exp(-cpar*x2);
  tem=t1-(1.d0-tem1)*(1.d0-tem2);
  tem3=cpar*tem1*tem2*t1
  pdf=tem3/(tem*tem);
  return
end subroutine dfrk


!compute the conditional pdf of Gumbel copula
subroutine pcondgum(x2,x1,cpar,ccdf)
  implicit none
  double precision::x2,x1,x,y,tem,sum
  double precision::cpar,tem1,tem2,ccdf
  x=-log(x1)
  y=-log(x2)
  tem1=x**cpar
  tem2=y**cpar
  sum=tem1+tem2
  tem=sum**(1.0d0/cpar)
  ccdf=exp(-tem)
  ccdf=ccdf*(1+tem2/tem1)**(-1.0d0+1.0d0/cpar)
  ccdf=ccdf/x1
end subroutine pcondgum


!compute the conditional pdf of frank copula
subroutine pcondfrk(x2,x1,cpar,ccdf) 
  implicit none
  double precision:: x2,x1
  double precision :: cpar,cpar1,tem,ccdf

  cpar1=1.d0-exp(-cpar)  ! 1 in double precision is 1.d0 in fortran
  tem=1.d0-exp(-cpar*x1) ! 1.d0
  ccdf=(1.d0-tem)/(cpar1/(1.d0-exp(-cpar*x2))-tem) ! 1.d0
end subroutine pcondfrk
