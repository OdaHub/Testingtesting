!BM---------------------------------------
module forcefieldIOT

contains

!BF------------------------------------------
function forceIOT(x,syst,IOTparam)
use types;use assign_objects;use alloc_objects
!use arimetric
implicit none
type(system_type),       intent(in)::syst
type(potiot_type),       intent(in)::iotparam
double precision, intent(in)::x(syst%Npart,2)
double precision::forceIOT(syst%Npart,2)
double precision::A(2),B(2),LBOX,eps,sig,eps2,sig2,r0,r02
integer::i,j,N
logical::cooperativity
double precision::coordination,Rcoop,ND,NN,NC,fcoop
double precision::FORCEcoop(syst%Npart,2),DFDC,Rij,RBA(2)
double precision::pf(2),fclust(3,3,2),FC,Vclust,Vshift

!!write(333,*) "force"

fclust=0.d0
forceIOT=0.d0

LBOX=syst%BOXLENGTH
N=syst%Npart

eps=IOTparam%eps
sig=IOTparam%sig
eps2=IOTparam%eps2
sig2=IOTparam%sig2
r0=IOTparam%r0
r02=IOTparam%r02


COOPERATIVITY=IOTparam%cooperativity
Rcoop=IOTparam%Rcoop
ND=IOTparam%ND
NN=IOTparam%NN
NC=IOTparam%NC
Vshift=IOTparam%Vshift
fcoop=IOTparam%fcoop



do i=1,N
  do j=1,N
    if (i==j) cycle
    A(:) =X(i,:); B(:)=X(j,:)
    pf=pair_forceA(A,B,i,j,LBOX,eps,sig,eps2,sig2,r0,r02)
    if ((i<=3).AND.(j<=3)) fclust(i,j,:)=pf(:)
    FORCEIOT(i,:)=FORCEIOT(i,:)+ pf(:) 
  enddo
enddo

!!do i=1,N
!!  write(333,*) FORCEIOT(i,:)
!!enddo 
!!write(333,*) "aaaa"

if (COOPERATIVITY) then

  FORCECOOP=0.d0
  coordination=Coordination_number(x,syst,iotparam)
  FC=FCOORD(COORDINATION,Fcoop,NC,NN)
  Vclust=0.d0

  !calculate force for i=2,3
  do i=2,3
    do j=1,3
      if (j==i) cycle
      FORCECOOP(i,:)=FORCECOOP(i,:)+fclust(i,j,:)
    enddo
  enddo 
  FORCECOOP(2:3,:)=(FC-1.d0)*FORCECOOP(2:3,:)

  !calculate Vclust
  do i=1,2
   do j=i+1,3
     A(:)=X(i,:);B(:)=X(j,:);rBA(:)=B(:)-A(:)
     rBA=rBA-LBOX*NINT(rBA/LBOX);Rij=sqrt( dot_product(rBA,rBA) )
     Vclust=Vclust+VpairIOT(Rij,i,j,eps,sig,eps2,sig2,r0,r02,Vshift)
   enddo
  enddo
!!  write(333,*) "Vclust=",Vclust
  
  !calculate force for i=4,...
  A(:)=X(1,:)
  do j=4,N
    B(:)=X(j,:)
    rBA(:)=B(:)-A(:)
    rBA=rBA-LBOX*NINT(rBA/LBOX)
    Rij=sqrt( dot_product(rBA,rBA) )
    !coordination=coordination+smoothstep(Rij,Rcoop,ND)
    FORCECOOP(j,:)=dsstep(Rij,Rcoop,ND)*RBA/Rij
    !FORCECOOP(1,:)=FORCECOOP(1,:)-FORCECOOP(j,:)
  enddo  
  DFDC=(1.d0-fcoop)*dsstep(coordination,NC,NN)
  FORCECOOP(4:N,:)=-(Vclust+eps2)*DFDC*FORCECOOP(4:N,:)

  !calculate force i=1
  FORCECOOP(1,1)= -sum(FORCECOOP(2:N,1))
  FORCECOOP(1,2)= -sum(FORCECOOP(2:N,2))

!!do i=1,N
!!  write(333,*) FORCECOOP(i,:)
!!enddo
!!write(333,*) "bbbb"


  FORCEIOT=FORCEIOT+FORCECOOP

endif
!!do i=1,N
!!  write(333,*) FORCEIOT(i,:)
!!enddo
!!write(333,*) "zzz"

end function forceIOT
!EF------------------------------------------------------

!BF------------------------------------------------------  
function EpotIOT(x,syst,IOTparam)
use types;use assign_objects;use alloc_objects
use arimetric
implicit none
type(system_type),       intent(in)::syst
type(potIOT_type),       intent(in)::IOTparam
double precision, intent(in)::x(syst%Npart,2)
double precision::EpotIOT
double precision::A(2),B(2),LBOX,eps,sig,r0,eps2,sig2,r02,Rij,Vshift
integer::i,j,N
logical::cooperativity
double precision::coordination,Rcoop,ND,NN,NC,fcoop,Vclust,Vp,Ecoop

!!write(333,*) "pot"

LBOX=syst%BOXLENGTH
N=syst%Npart
eps=IOTparam%eps
sig=IOTparam%sig
eps2=IOTparam%eps2
sig2=IOTparam%sig2
r0=IOTparam%r0
r02=IOTparam%r02
Vshift=IOTparam%Vshift
COOPERATIVITY=IOTparam%cooperativity
Rcoop=IOTparam%Rcoop
ND=IOTparam%ND
NN=IOTparam%NN
NC=IOTparam%NC
fcoop=IOTparam%fcoop

EpotIOT=0.d0
COORDINATION=0.d0
Vclust=0.d0
do i=1,N-1
  A(:)=X(i,:)
  do j=i+1,N
    B(:)=X(j,:)
    Rij=distance(A,B,2,LBOX)
    Vp=VpairIOT(Rij,i,j,eps,sig,eps2,sig2,r0,r02,Vshift)
    EpotIOT=EpotIOT+Vp
    if (COOPERATIVITY) then
      if (i==1) then
        if (j<=3) then
          Vclust=Vclust+Vp
        else
          COORDINATION=COORDINATION+smoothstep(Rij,Rcoop,ND)
        endif
      endif
      if ((i==2).AND.(j==3)) Vclust=Vclust+Vp
    endif
  enddo
enddo
!!write(333,*) EpotIOT, "aaa"

if (COOPERATIVITY) then
  Ecoop=(FCOORD(COORDINATION,Fcoop,NC,NN)-1.d0)* (Vclust+eps2)
  EpotIOT=EpotIOT + Ecoop
!!  write(333,*) Ecoop, Epotiot,"bbb:"
endif

end function EpotIOT
!EF------------------------------------------------------

!BF--------------------------------------------------------
function LJ(eps,sig,x)
implicit none
double precision, intent(in)::eps,sig,x
double precision::LJ
  LJ=4*eps*((sig/x)**12-(sig/x)**6)
end function LJ
!EF--------------------------------------------------------

!BF--------------------------------------------------------------------
function pair_forceA(A,B,i,j,LBOX,eps,sig,eps2,sig2,r0,r02)
implicit none
integer, intent(in)::i,j
double precision, intent(in)::A(2),B(2)
double precision, intent(in)::LBOX,eps,sig,eps2,sig2,r0,r02
double precision::pair_forceA(2)
double precision::RAB(2),dAB

!i=/j
pair_forceA=0.d0
rAB(:)=A(:)-B(:)
rAB=rAB-LBOX*NINT(rAB/LBOX)
dAB=sqrt( dot_product(rAB,rAB) )

if (  ( (i==1).AND.((j==2).OR.(j==3)) ) &
  .OR.( (j==1).AND.((i==2).OR.(i==3)) )  ) then
  pair_forceA=-dLJ(dAB,eps2,sig2,r0)*RAB(:)/dAB
else if  ( ((i==1).AND.(j>3)).OR.((j==1).AND.(i>3))  ) then
  pair_forceA=-dLJ(dAB,eps2,sig2,r02)*RAB(:)/dAB
else
  pair_forceA=-dLJ(dAB,eps,sig,r0)*RAB(:)/dAB
endif

end function pair_forceA
!EF--------------------------------------------------------------------

!BF-----------------------------------------------------
function dLJ(r,eps,sig,r0)
implicit none
double precision, intent(in)::r,eps,sig,r0
double precision::dLJ

dLJ=0.d0
if (r< r0) then
  dLJ=4*eps*( -12*(sig**12/r**13) +6*(sig**6/r**7) )
endif

end function dLJ
!EF----------------------------------------------------

!BF-----------------------------------------------------
function VpairIOT(r,i,j,eps,sig,eps2,sig2,r0,r02,Vs)
implicit none
double precision, intent(in)::r,eps,sig,sig2,eps2,r0,r02,Vs
integer, intent(in)::i,j
double precision::VpairIOT

!i<j 
VpairIOT=0.d0
if (i==1) then
  if ((j==2).OR.(j==3)) then
    if (r<r0) VpairIOT=LJ(eps2,sig2,r)-Vs
  else
    if (r<r02) VpairIOT=LJ(eps2,sig2,r)+eps2 
  endif
else 
  if (r<r0) VpairIOT=LJ(eps,sig,r)+eps
endif

end function VpairIOT
!EF----------------------------------------------------

!BF-----------------------------------------------------
function smoothstep(Rij,Rcoop,ND)
implicit none
double precision, intent(in)::Rij,Rcoop,ND
double precision::smoothstep

  smoothstep=1.d0/(1.d0+exp(ND*(Rij-Rcoop)))

end function smoothstep
!EF------------------------------------------------------

!BF--------------------------------------------------------
FUNCTION FCOORD(C,Fcoop,NC,NN)
implicit none
double precision, intent(in)::C,Fcoop,NC,NN
double precision::FCOORD
  FCOORD=Fcoop+(1.D0-Fcoop)*smoothstep(C,NC,NN) 
end FUNCTION FCOORD
!EF--------------------------------------------------------

!BF-----------------------------------------------------
function dsstep(Rij,Rcoop,ND)
implicit none
double precision, intent(in)::Rij,Rcoop,ND
double precision::dsstep
double precision::expf
  expf=exp(ND*(Rij-Rcoop))
  dsstep=-ND*expf/((1.d0+expf)**2)

end function dsstep
!EF------------------------------------------------------


!BF-----------------------------------------------------
function Coordination_number(x,syst,iot)
use types;use assign_objects;use alloc_objects
use arimetric
implicit none
type(system_type),       intent(in)::syst
type(potiot_type),       intent(in)::iot
double precision, intent(in)::x(syst%Npart,2)
double precision::Coordination_number
double precision::A(2),B(2),LBOX,Rij
integer::j,N
double precision::Rcoop,ND

LBOX=syst%BOXLENGTH
N=syst%Npart
Rcoop=IOT%Rcoop
ND=IOT%ND

Coordination_number=0.d0
A(:)=X(1,:)
do j=4,N
  B(:)=X(j,:)
  Rij=distance(A,B,2,LBOX)
  COORDINATION_number=COORDINATION_number+smoothstep(Rij,Rcoop,ND)
enddo

end function Coordination_number
!EF------------------------------------------------------


end module forcefieldIOT
!EM-------------------------------------------------------------
