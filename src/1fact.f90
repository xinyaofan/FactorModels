!test main program
! program test 
!  implicit none
!  integer npar,d,n,ip,i
!  double precision nllk
!  integer,dimension(:),allocatable::fam
!  double precision,dimension(:), allocatable :: lgrad,lat,param
!  double precision,dimension(:,:), allocatable :: udata,lhess
!  
!  read*,n
!  read*,d
!  
!  allocate(udata(n,d),param(d),lgrad(d),lhess(d,d),lat(n))
!  do i=1,n
!  read*,udata(i,:)
!  end do
!  
!  read*,(lat(i),i=1,n)
!  
!  
!  allocate(fam(d))
!  read*,fam
!  
!  read*,(param(ip),ip=1,d)
!  
! 
!  
!  call onefact(param,d,n,udata,lat,fam,nllk,lgrad,lhess)
!  
!  print*,lgrad
!  print "(10f10.5)", lhess
!  print*,nllk
!  
!  deallocate(udata,param,fam,lat,lgrad,lhess)
!  end

! 1-factor model where bivariate linking copulas are one-parameter family
! inputs 
!   fam = copula families in the first tree 
!   npar = #parameters = d
!   param0 = parameter vector (dimension d=npar), 
!   d = #variables
!   n = sample size
!   udata = nxd matrix of uniform scores
!  lat= a data vector of length n
! outputs 
!   nllk = negative log-likelihood, 
!   lgrad = gradient of nllk, 
!   lhess = hessian of nllk
subroutine onefact(param,d,n,udata,lat,fam,nllk,lgrad,lhess)
  implicit none
  integer d,n,fam(d),i,j
  double precision nllk,lgrad(d),lhess(d,d),lsumpdf,param(d),lat(n),udata(n,d)
  double precision, dimension(:), allocatable :: uvec,grad,lpdf
  double precision, dimension(:,:),allocatable ::hess,der1,der2
  double precision ljpdf
  double precision lder1u,lder1v,lder2u,lder2v,lder2uv
  double precision ldermixuu(2),ldermixvv(2)
  
  allocate(uvec(d),grad(d),lpdf(d),hess(d,d),der1(2,d),der2(3,d))
  nllk=0.d0; lgrad=0.d0; lhess=0.d0; lsumpdf=0.0d0;! initialize
  do i=1,n ! loop over rows of data set
    uvec=udata(i,:)
    !print*,lat(i)
    grad=0.0d0;hess=0.0d0
    do j=1,d
     call lcop2derivt(lat(i),uvec(j),param(j),fam(j),lpdf(j),der1(:,j),&
      der2(:,j),lder1u,lder2u,ldermixuu,lder1v,lder2v,lder2uv,ldermixvv)

    grad(j)=grad(j)+der1(1,j)
    hess(j,j)=hess(j,j)+der2(1,j)
    end do
    ljpdf=sum(lpdf)
    lsumpdf=lsumpdf+ljpdf
    
    lgrad=lgrad-grad
    lhess=lhess-hess
  end do
  nllk=-lsumpdf
  deallocate(uvec,grad,lpdf,hess,der1,der2)
  return
  end


! 1-factor model where bivariate linking copulas allowing two-parameter family
! inputs 
!   fam = copula families in the first tree 
!   npar = #parameters = d
!   param: param matrix 2*d, j-th column is param for j-th linking copulas, 
!   d = #variables
!   n = sample size
!   udata = nxd matrix of uniform scores
!  lat= a data vector of length n
! outputs 
!   nllk = negative log-likelihood, 
!   lgrad = gradient of nllk, 
!   lhess = hessian of nllk
subroutine onefact2(param,d,n,udata,lat,fam,nllk,lgrad,lhess)
   implicit none
  integer d,n,fam(d),i,j,jj(2)
  double precision nllk,lgrad(2*d),lhess(2*d,2*d),lsumpdf,param(2,d),lat(n),udata(n,d)
  double precision, dimension(:), allocatable :: uvec,grad,lpdf
  double precision, dimension(:,:),allocatable ::hess,der1,der2
  double precision ljpdf
  double precision tmat(2,2)
  double precision lder1u,lder1v,lder2u,lder2v,lder2uv
  double precision ldermixuu(2),ldermixvv(2)


  allocate(uvec(d),grad(2*d),lpdf(d),hess(2*d,2*d),der1(2,d),der2(3,d))
  nllk=0.d0; lgrad=0.d0; lhess=0.d0; lsumpdf=0.0d0;! initialize
  do i=1,n ! loop over rows of data set
    uvec=udata(i,:)
    !print*,lat(i)
    grad=0.0d0;hess=0.0d0
    do j=1,d
      jj=(/2*j-1,2*j/)
      call lcop2derivt(lat(i),uvec(j),param(:,j),fam(j),lpdf(j),der1(:,j),&
      der2(:,j),lder1u,lder2u,ldermixuu,lder1v,lder2v,lder2uv,ldermixvv)
      grad(jj)=grad(jj)+der1(:,j)
      !lgrad(jj)=lgrad(jj)+der1(:,j)
      !print*,lgrad(jj)
      tmat(1,1)=der2(1,j); tmat(2,2)=der2(3,j);
      tmat(1,2)=der2(2,j); tmat(2,1)=der2(2,j);
      hess(jj,jj)=hess(jj,jj)+tmat
    end do
  
    ljpdf=sum(lpdf)
    lsumpdf=lsumpdf+ljpdf
    
    lgrad=lgrad-grad
    lhess=lhess-hess
  end do
  nllk=-lsumpdf
  deallocate(uvec,grad,lpdf,hess,der1,der2)
  return
  end


  