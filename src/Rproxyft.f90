! program test 
!  implicit none
!  integer n,d,nq,i,iq,j,npar,ip
!  integer, dimension(:), allocatable :: edg1,edg2
!  integer,dimension(:), allocatable :: fam
!  double precision, dimension(:,:), allocatable :: udata,hess
!  double precision, dimension(:), allocatable :: grad,param0,lat
!  double precision nllk
!  !read *,nq
!  !allocate ( xl(nq),wl(nq) )
!  
!  !do iq=1,nq
!  !xl(iq)=iq/(nq+1.0d0)
!  !wl(iq)=1.0d0/nq
!  !end do
!  
! ! read *, xl(:)
!  !read *, wl(:)
!  read *,n
!  read*,d
!  npar=d*4-2
!  allocate (udata(n,d),lat(n), param0(npar), grad(npar), hess(npar,npar) )
!  allocate (edg1(d-1), edg2(d-1) )
!  do i=1,n
!    read *, udata(i,:)
!  end do
!  
!  do i=1,n
!  read*,lat(i)
!  end do
!  
!  allocate(fam(2*d-1))
!  read*,fam
!  
!  read *, (param0(ip),ip=1,npar)
!  read *, edg1(:)
!  read *, edg2(:)
! 
! !  print*,n
! !  print*,d
! !   print*,edg1
! !  print*,edg2
! !  print*,udata(1,:)
! !  print*,fam
! !  print*,param0
! !print*,wl
! !print*,xl
!  call ft2(npar,param0,d,n,udata,fam,edg1,edg2, nllk,grad,hess)
!  print *, nllk
!  print*,grad
!  !print "(8f10.5)", grad
!  print *, "grad above, hess below"
!  !print *," "
!  !print*,hess
!  print "(10f10.5)", hess
!  deallocate (param0,fam,udata,lat,grad,hess,edg1,edg2)
!  stop
!  end

! 1-factor 1-truncated model
! inputs 
!   fam = copula families in the first tree (in an natural order)
!   npar = #parameters = 4*d-2
!   param0 = parameter vector (dimension d=npar), 
!       param0 has th1,dl1,th2,dl2, ... as vector theta>0, delta>1
!   d = #variables
!   n = sample size
!   udata = nxd matrix of uniform scores
!   nq = #quadrature points
!   wl = vector of quadrature weights
!   xl = vector of quadrature nodes
!   edg1 = node 1 vector of edges on residual tree  (length d-1)
!   edg2 = node 2 vector of edges on residual tree 
! outputs 
!   nllk = negative log-likelihood, 
!   grad = gradient of nllk, 
!   hess = hessian of nllk
subroutine ft2(npar,param0,d,n,udata,lat,fam,edg1,edg2,nllk,grad,hess)
  implicit none
  integer npar,d,n,edg1(d-1),edg2(d-1),fam(2*d-1)
  double precision param0(npar),udata(n,d),lat(n)
  double precision nllk,grad(npar),hess(npar,npar)
  integer i,j,j2,jj(2),jj2(2),jj3(2)
  integer d2,d1,ie,e1,e2,jj1(2)
  integer, dimension(:), allocatable ::famvec1,famvec2
  double precision, dimension(:), allocatable :: uvec,lpdf,fval1,grd
  double precision, dimension(:,:), allocatable :: cder1,cder2
  double precision, dimension(:), allocatable :: ccdf
  double precision, dimension(:,:), allocatable :: fval2,hss,der1,der2
  double precision, dimension(:,:),allocatable :: partr,param
  double precision, dimension(:), allocatable :: lgrad
  double precision, dimension(:,:), allocatable :: lhess,gmat
  double precision fval,ljpdf,tmat(2,2)
  double precision lpdfr,lder1u,lder1v,lder2u,lder2v,lder2uv
  double precision ldermixuu(2),ldermixvv(2),ldercpar(2),ldercpar2(3)
  ! npar=4*d-2
  allocate ( uvec(d),lpdf(d), der1(2,d), der2(3,d),fval1(npar), grd(npar))
  allocate ( fval2(npar,npar), hss(npar,npar), param(2,d) )
  allocate ( ccdf(d), cder1(2,d), cder2(3,d) )
  allocate (lgrad(npar), lhess(npar,npar), gmat(npar,npar), partr(2,d-1) )
  allocate(famvec1(d),famvec2(d-1))
  d2=2*d; d1=d-1
  famvec1=fam(1:d)
  famvec2=fam((d+1):(d2-1))
  param=reshape(param0(1:d2),(/2,d/)) !param for the frist tree
  partr=reshape(param0((d2+1):npar),(/2,d-1/))
  
  nllk=0.d0; grad=0.d0; hess=0.d0; ! initialize
  do i=1,n ! loop over rows of data set
    uvec=udata(i,:)
    !integl=0.d0; integl1=0.d0; integl2=0.d0; ! initialize for integrals
    !do iq=1,nq ! loop over quadrature points
        lgrad=0.d0; lhess=0.d0; 
      ! get \p log copden(u[j],v,th[j]) / \p th[j] , j=1,...,d
      do j=1,d
        call lcop2derivt(lat(i),uvec(j),param(:,j),famvec1(j),lpdf(j),der1(:,j),&
          der2(:,j),lder1u,lder2u,ldermixuu,lder1v,lder2v,lder2uv,ldermixvv)
        call ccopderiv(uvec(j),lat(i),param(:,j),famvec1(j),ccdf(j),cder1(:,j),cder2(:,j))
      end do
      ! update integrand, and derivs wrt copula parameters 
      ! fval1: vector of partial deriv wrt param[j], j=1,...,np
      ! fval2: matrix of 2nd partial deriv wrt param[j], param[j2], 
      ljpdf=sum(lpdf)
      !print*,ljpdf
      do j=1,d
        jj=(/2*j-1,2*j/)
        lgrad(jj)=lgrad(jj)+der1(:,j)
        !print*,lgrad(jj)
        tmat(1,1)=der2(1,j); tmat(2,2)=der2(3,j);
        tmat(1,2)=der2(2,j); tmat(2,1)=der2(2,j);
        lhess(jj,jj)=lhess(jj,jj)+tmat
      end do
     !print*,lgrad(1:2)
      do ie=1,d1
        e1=edg1(ie); e2=edg2(ie)
        jj1=(/2*e1-1,2*e1/); jj2=(/2*e2-1,2*e2/); jj3=(/2*ie-1,2*ie/);
        !tem=lgrad(jj1) 
       ! print*,"tem now"
        !print*,tem
        call lcop2derivt(ccdf(e1),ccdf(e2),partr(:,ie),famvec2(ie),&
        lpdfr,ldercpar,ldercpar2,lder1u,lder2u,ldermixuu,&
        lder1v,lder2v,lder2uv,ldermixvv)
        !print*,"tem"
       ! print*,tem
        lgrad(jj1)=lgrad(jj1) +lder1u*cder1(:,e1)
        lgrad(jj2)=lgrad(jj2)+lder1v*cder1(:,e2)
        lgrad(d2+jj3)=lgrad(d2+jj3)+ldercpar!correct
        ljpdf=ljpdf+lpdfr
    
       !compute hessian
          call outer(2,2,cder1(:,e1),cder1(:,e1),tmat)
          lhess(jj1,jj1)=lhess(jj1,jj1)+lder2u*tmat
     
          call outer(2,2,cder1(:,e2),cder1(:,e2),tmat)
          lhess(jj2,jj2)=lhess(jj2,jj2)+lder2v*tmat
        
          call outer(2,2,cder1(:,e1),cder1(:,e2),tmat)
          lhess(jj1,jj2)=lhess(jj1,jj2)+lder2uv*tmat
          
          call outer(2,2,cder1(:,e2),cder1(:,e1),tmat)
          lhess(jj2,jj1)=lhess(jj2,jj1)+lder2uv*tmat
          
          tmat(1,1)=cder2(1,e1); tmat(2,2)=cder2(3,e1);
          tmat(1,2)=cder2(2,e1); tmat(2,1)=cder2(2,e1);
          lhess(jj1,jj1)=lhess(jj1,jj1)+lder1u*tmat
          
          tmat(1,1)=cder2(1,e2); tmat(2,2)=cder2(3,e2);
          tmat(1,2)=cder2(2,e2); tmat(2,1)=cder2(2,e2);
          lhess(jj2,jj2)=lhess(jj2,jj2)+lder1v*tmat
          
          tmat(1,1)=ldercpar2(1);tmat(2,2)=ldercpar2(3);
          tmat(1,2)=ldercpar2(2);tmat(2,1)=ldercpar2(2);
          lhess(d2+jj3,d2+jj3)=lhess(d2+jj3,d2+jj3)+tmat
       
          call outer(2,2,ldermixuu,cder1(:,e1),tmat)
          
          lhess(jj1,d2+jj3)=lhess(jj1,d2+jj3)+transpose(tmat)
          lhess(d2+jj3,jj1)=lhess(d2+jj3,jj1)+tmat
          
          call outer(2,2,ldermixvv,cder1(:,e2),tmat)    
          lhess(jj2,d2+jj3)=lhess(jj2,d2+jj3)+transpose(tmat)
          lhess(d2+jj3,jj2)=lhess(d2+jj3,jj2)+tmat     
      ! end do
      ! update quadrature loops
       fval=exp(ljpdf)
       fval1=fval*lgrad
       call outer(npar,npar,lgrad,lgrad,gmat)
       !print*,gmat
       fval2=fval*(lhess+gmat)
       !ww=wl(iq)
       !integl=integl+fval*ww
       !integl1=integl1+fval1*ww
       !integl2=integl2+fval2*ww
    end do
    ! update contribution to negative log-likelihood
     !if(integl<=0.d0) integl=1.d-100
     nllk=nllk-log(fval);
     grd=fval1/fval;grad=grad-grd;
     !grd=integl1/integl; grad=grad-grd;
     do j=1,npar
       do j2=1,npar
         hss(j,j2)=fval2(j,j2)/fval-grd(j)*grd(j2);
         !hss(j,j2)=integl2(j,j2)/integl - grd(j)*grd(j2); 
       end do
     end do
     hess=hess-hss
  end do
  deallocate (uvec,lpdf,der1,der2,fval1,grd)
  deallocate (fval2,hss,param)
  deallocate ( ccdf,cder1,cder2)
  deallocate ( lgrad,lhess,gmat,partr )
  deallocate(famvec1,famvec2)
  return
  end

! outer product of two vectors
! input
!   na = length of avec
!   nb = length of bvec
!   avec = vector 1
!   bvec = vector 2
! output
!   abmat = na x nb matrix with outer product of avec and bvec 
subroutine outer(na,nb,avec,bvec,abmat)
  implicit none
  integer na,nb,i
  double precision avec(na),bvec(nb),abmat(na,nb)
  ! na=size(avec); nb=size(bvec)
  do i =1,na  
    abmat(i,:)=avec(i)*bvec
  end do
  return
  end