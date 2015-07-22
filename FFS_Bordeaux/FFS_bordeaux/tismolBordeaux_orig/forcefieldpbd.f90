!BM---------------------------------------
module forcefieldPBD

contains

 !BF------------------------------------------
  function forcePBDshift(x,syst,PBDparam)
  use types;use assign_objects;use alloc_objects
  implicit none
  type(system_type),       intent(in)::syst
  type(potPBD_type),       intent(in)::PBDparam
  double precision, intent(in)::x(syst%Npart,1)
  double precision::forcePBDshift(syst%Npart,1)
  double precision::y(syst%Npart,1)
  integer::i,imin
  double precision::xmin


  forcePBDshift=0.d0
  imin=1
  xmin=x(1,1)
  do i=2,syst%Npart
    if (x(i,1) < xmin) then
      xmin=x(i,1)
      imin=i
    endif
  enddo
  y(:,1)=x(:,1)-xmin+PBDparam%surface
  forcePBDshift=forcePBD(y,syst,PBDparam)
  forcePBDshift(imin,1)=forcePBDshift(imin,1)-&
                        SUM(forcePBDshift(1:syst%Npart,1))

  end function forcePBDshift
  !EF------------------------------------------------------

  !BF------------------------------------------
  function forcePBD(x,syst,PBDparam)
  use types;use assign_objects;use alloc_objects
  implicit none
  type(system_type),       intent(in)::syst
  type(potPBD_type),       intent(in)::PBDparam
  double precision, intent(in)::x(syst%Npart,1)
  double precision::forcePBD(syst%Npart,1)
  double precision::Y(0:syst%Npart+1),x1D(syst%Npart)
  integer::N
  double precision::a(syst%Npart),da2(syst%Npart),rho,alpha,al_2,S
  double precision::eexp, dyr, syr, dyl,syl,facr,facl
  integer::i
  logical::PBC,BIAS

  if (syst%dim/=1) stop !remove after test
  N=syst%Npart
  PBC=syst%PBC
  a=PBDparam%a
  da2=PBDparam%da2
  rho=PBDparam%rho
  alpha=PBDparam%alpha
  al_2=PBDparam%al_2
  S=PBDparam%S
  BIAS=PBDparam%BIAS
  x1d(1:N)=x(1:N,1)

  call extend_chain(x1D,N,PBC,Y)

  !Now loop over base pairs without concerning boundary problems
  do i=1,N
    eexp=exp(-a(i)*Y(i))
    dyr=Y(i+1)-Y(i)
    syr=Y(i+1)+Y(i)
    dyl=Y(i-1)-Y(i)
    syl=Y(i-1)+Y(i)
    
    facr= 1.d0 + rho*exp(-alpha*syr)*(1.d0 + al_2*dyr)
    facl= 1.d0 + rho*exp(-alpha*syl)*(1.d0 + al_2*dyl)
 
    forcePBD(i,1)=S*( dyr*facr + dyl*facl ) + da2(i)*eexp*(eexp - 1.d0)
  enddo
  if (BIAS) forcePBD=forcePBD+FBIAS(x,N,PBDparam)
 
  end function forcePBD
  !EF------------------------------------------

 !BF------------------------------------------
  function fbias(x,N,PBDparam)
  use types;use assign_objects;use alloc_objects
  implicit none
  integer,       intent(in)::N
  type(potPBD_type),       intent(in)::PBDparam
  double precision, intent(in)::x(N,1)
  double precision::fbias(N,1)
  integer::iminarr(1:1),imin
  double precision::Xmin


  fbias=0.d0
  iminarr=minloc(X(:,1));imin=iminarr(1)
  Xmin=X(imin,1)
  if (Xmin>PBDparam%BIASCUTOFF) then
    fbias(imin,1)=-PBDparam%BIASPREFAC*PBDparam%BIASEXP*(Xmin-PBDparam%BIASCUTOFF)**(PBDparam%BIASEXP-1.d0)
  endif

  end function fbias
  !EF------------------------------------------

  !BF-------------------------------------------
  function EpotPBD(x,syst,PBDparam)
  use types;use assign_objects;use alloc_objects
  implicit none
  type(system_type),       intent(in)::syst
  type(potPBD_type),       intent(in)::PBDparam
  double precision, intent(in)::x(syst%Npart,1)
  double precision::EpotPBD
  double precision::Y(0:syst%Npart+1),x1D(syst%Npart)
  integer::N
  double precision::a(syst%Npart),d(syst%Npart),rho,alpha,S_2
  double precision::ecoupl,Esite,facharm,facanh
  integer::i
  logical::PBC,BIAS

  N=syst%Npart
  PBC=syst%PBC
  a=PBDparam%a
  d=PBDparam%d
  rho=PBDparam%rho
  alpha=PBDparam%alpha
  S_2=PBDparam%S_2
  BIAS=PBDparam%BIAS
  x1d(1:N)=x(1:N,1)
  call extend_chain(x1D,N,PBC,Y)
  
 ecoupl=0.d0;Esite=0.d0
  !Now loop over base pairs without concerning boundary problems
  do i=1,N
    facharm=  (Y(i) - Y(i+1))**2
    facanh =  1.d0  + rho*exp(  -alpha*( Y(i) + Y(i+1) )  )
    ecoupl =  ecoupl+ facharm* facanh
    esite  =  esite + d(i)*(exp(-a(i)*Y(i)) - 1.d0)**2
  enddo
  ecoupl= S_2*ecoupl
  EpotPBD =  Ecoupl+ Esite

   if (BIAS) EpotPBD=EpotPBD+E_BIAS(x,N,PBDparam)

  end function EpotPBD
  !EF-------------------------------------------

!BF---------------------------------------------------------------------------
! This function gives the potential energy due to the bias-force
!-----------------------------------------------------------------------------
function E_bias(x,N,PBDparam)
use types;use assign_objects;use alloc_objects
implicit none
integer,       intent(in)::N
type(potPBD_type),       intent(in)::PBDparam
double precision, intent(in)::x(N,1)
double precision::E_bias
integer::iminarr(1:1),imin
double precision::Xmin

  E_bias=0.d0
  iminarr=minloc(X(:,1));imin=iminarr(1)
  Xmin=X(imin,1)
  if (Xmin>PBDparam%BIASCUTOFF) then
    E_bias=PBDparam%BIASPREFAC*(Xmin-PBDparam%BIASCUTOFF)**PBDparam%BIASEXP
  endif

end function E_bias
!EF---------------------------------------------------------------------------


  !BS-------------------------------------------
  subroutine extend_chain(x,N,PBC,Y)
  implicit none
  integer, intent(in)::N
  double precision, intent(in)::x(N)
  logical, intent(in)::PBC
  double precision, intent(out)::Y(0:N+1)

  Y(1:N)=x(1:N)
  if (PBC) then
    Y(0)=x(N)
    Y(N+1)=x(1)
  else
    Y(0)=x(1)
    Y(N+1)=x(N)
  endif

  end subroutine extend_chain
  !ES-------------------------------------------

end module forcefieldPBD
!EM---------------------------------------

