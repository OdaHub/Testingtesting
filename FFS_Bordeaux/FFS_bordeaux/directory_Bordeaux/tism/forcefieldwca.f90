!BM---------------------------------------
module forcefieldWCA

contains

 !BF------------------------------------------
  function forceWCA(x,syst,WCAparam)
  use types;use assign_objects;use alloc_objects
  implicit none
  type(system_type),       intent(in)::syst
  type(potWCA_type),       intent(in)::WCAparam
  double precision, intent(in)::x(syst%Npart,2)
  double precision::forceWCA(syst%Npart,2)
  double precision::A(2),B(2),LBOX,epsilon,sigma,h,w,r0
  integer::i,j,N

  forceWCA=0.d0

  LBOX=syst%BOXLENGTH
  N=syst%Npart

  epsilon=WCAparam%epsilon
  sigma=WCAparam%sigma
  h=WCAparam%h
  w=WCAparam%w
  r0=WCAparam%r0

  do i=1,N
    do j=1,N
      if (i==j) cycle
      A(:) =X(i,:); B(:)=X(j,:)
      FORCEWCA(i,:)=FORCEWCA(i,:)+ & 
                    pair_forceA(A,B,i,j,LBOX,epsilon,sigma,h,w,r0)
    enddo
  enddo

  end function forceWCA
  !EF------------------------------------------------------

  !BF------------------------------------------------------  
  function EpotWCA(x,syst,WCAparam)
  use types;use assign_objects;use alloc_objects
  use arimetric
  implicit none
  type(system_type),       intent(in)::syst
  type(potWCA_type),       intent(in)::WCAparam
  double precision, intent(in)::x(syst%Npart,2)
  double precision::EpotWCA
  double precision::A(2),B(2),LBOX,epsilon,sigma,h,w,r0,Rij
  integer::i,j,N

  LBOX=syst%BOXLENGTH
  N=syst%Npart
  epsilon=WCAparam%epsilon
  sigma=WCAparam%sigma
  h=WCAparam%h
  w=WCAparam%w
  r0=WCAparam%r0

  EpotWCA=0.d0
  do i=1,N-1
    A(:)=X(i,:)
    do j=i+1,N
      B(:)=X(j,:)
      Rij=distance(A,B,2,LBOX)
      EpotWCA=EpotWCA+Vpair(Rij,i,j,epsilon,sigma,h,w,r0)
    enddo
  enddo
  end function EpotWCA
  !EF------------------------------------------------------

  !BF--------------------------------------------------------------------
  function pair_forceA(A,B,i,j,LBOX,epsilon,sigma,h,w,r0)
  implicit none
  integer, intent(in)::i,j
  double precision, intent(in)::A(2),B(2)
  double precision, intent(in)::LBOX,epsilon,sigma,h,w,r0
  double precision::pair_forceA(2)
  double precision::RAB(2),dAB

  pair_forceA=0.d0
  rAB(:)=A(:)-B(:)
  rAB=rAB-LBOX*NINT(rAB/LBOX)
  dAB=sqrt( dot_product(rAB,rAB) )

  pair_forceA=-dVWCA(dAB,epsilon,sigma,r0)*RAB(:)/dAB

  if ( ( (i==1).AND.(j==2) ).OR.( (i==2).AND.(j==1) ) ) then
    pair_forceA=pair_forceA-dVdw(dAB,h,w,r0)*RAB(:)/dAB
  endif

  end function pair_forceA
!EF--------------------------------------------------------------------

!BF-----------------------------------------------------
function dVwca(r,epsilon,sigma,r0)
implicit none
double precision, intent(in)::r,epsilon,sigma,r0
double precision::dVwca

dVwca=0.d0
if (r< r0) then
  dVwca=4*epsilon*( -12*sigma*(r/sigma)**(-13) +6*sigma* (r/sigma)**(-7) )
endif
!!think there is an error. Check for sigma/=1 !!
end function dVwca
!EF----------------------------------------------------

!BF--------------------------------------------------------------------
function dVdw(r,h,w,r0) ! derivative dV/dr of Symmetric Double Well Potential
implicit none
double precision, intent(in)::r,h,w,r0
double precision::dVdw


dVdw=-4*h*(  1-( (r-r0-w)**2/w**2 )  )*   ( (r-r0-w)/w**2 )

end function dVdw
!EF--------------------------------------------------------------------

!BF-----------------------------------------------------
function Vpair(r,i,j,epsilon,sigma,h,w,r0)
implicit none
double precision, intent(in)::r,epsilon,sigma,h,w,r0
integer, intent(in)::i,j
double precision::Vpair

Vpair=Vwca(r,epsilon,sigma,r0)
if ( ( (i==1).AND.(j==2) ).OR.( (i==2).AND.(j==1) ) ) then
  Vpair=Vpair+Vdw(r,h,w,r0)
endif
end function Vpair
!EF----------------------------------------------------

!BF-----------------------------------------------------
function Vwca(r,epsilon,sigma,r0)
implicit none
double precision, intent(in)::r,epsilon,sigma,r0
double precision::Vwca

Vwca=0.d0
if (r< r0) then
  Vwca=4*epsilon*( (r/sigma)**(-12) - (r/sigma)**(-6) ) + epsilon
endif

end function Vwca
!EF----------------------------------------------------

!BF-----------------------------------------------------
function Vdw(r,h,w,r0)
implicit none
double precision, intent(in)::r,h,w,r0
double precision::Vdw

Vdw=h*(  1.d0 - (r-r0-w)**2/w**2  )**2
end function Vdw
!EF----------------------------------------------------




end module forcefieldWCA
