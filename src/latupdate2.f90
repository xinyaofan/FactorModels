program test
  implicit none
  integer family,nq,n,d,i
  double precision,dimension(:,:),allocatable::udata
  double precision,dimension(:),allocatable::xl,wl,lat,th

  read*,nq
  print*,nq
  allocate(xl(nq),wl(nq))
  read*,xl
  read*,wl
  
  print*,xl
  print*,wl
  
  read*,n
  read*,d
  allocate(udata(n,d),th(d),lat(n))
  read *, (udata(i,:), i=1,n) 
  
  read*,th
  print*,th
  read*,family
  print*,family
  
  call latupdate2 (th,n,d,udata,nq,xl,wl,family,lat)
  
  print*,lat
  
  deallocate(udata,xl,wl)
end program 



!latent update in one-factor copula model; different copula families are allowed

subroutine latupdate2 (th,n,d,udata,nq,xl,wl,family,lat)
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
