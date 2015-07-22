!BM---------------------------------------
module forcefield

contains

  !BF------------------------------------------
  function force(x,syst,pot)
  use types;use assign_objects;use alloc_objects
  use forcefieldPBD
  use forcefieldNicolis
  use forcefieldWCA
  use forcefieldiot
  use forcefieldharm
  implicit none
  type(system_type),       intent(in)::syst
  type(potential_type),    intent(in)::pot
  double precision, intent(in)::x(syst%Npart,syst%dim)
  double precision::force(syst%Npart,syst%dim)

  select case(pot%POTENTIAL)
    case("PBD")
      force=forcePBD(x,syst,pot%potPBD)
    case("PBDshift")
      force=forcePBDshift(x,syst,pot%potPBD)
    case("HARMOSC")
      force=-pot%potHARMOSC%k*x
    case("DOUBLEWELL")
      force(1,1)=-4*pot%potDW%k4*x(1,1)**3-2*pot%potDW%k2*x(1,1)
    case("2DHARM")
      force=force2dharm(x(1,1),x(1,2),pot%pot2DHARM)
    case("NICOLIS1D")
      force=forceNICOLIS1D(x,syst,pot%potNICOLIS)
    case("NICOLIS2D")
      force=forceNICOLIS2D(x,syst,pot%potNICOLIS)
    case("WCA")
      force=forceWCA(x,syst,pot%potWCA)
    case("IONTRANS")
      force=forceIOT(x,syst,pot%potiot)
    case default
      print *,"ERROR force: POTENTIAL=",pot%POTENTIAL
      stop
  end select
end function force 
!EF------------------------------------------

  !BF------------------------------------------
  function EKIN(v,masses_2,N,dim)
  implicit none
  integer, intent(in)::N,dim
  double precision, intent(in)::v(N,dim),masses_2(N)
  double precision::EKIN
  integer::i
  double precision::vec(dim)
  
  EKIN=0.d0
  do i=1,N
    vec(1:dim)=V(i,1:dim)
    EKIN=EKIN+masses_2(i)*dot_product(vec,vec)
  enddo

end function EKIN 
!EF------------------------------------------

!BF-----------------------------------------------------------------
function E_NOSE_THERMOS(extcoord,NHparam,syst)
use types;use assign_objects;use alloc_objects
implicit none
type(extended_coordinates_type),intent(in)::extcoord
type(NOSE_HOOVER_type),intent(in)::NHparam
type(system_type), optional::syst
double precision::E_NOSE_THERMOS
integer::i,NT
double precision::Q_Nose_2(NHparam%Nthermo)
double precision::xi(NHparam%Nthermo),vxi(NHparam%Nthermo)
double precision::NFkbT,kBT

NT=NHparam%Nthermo
Q_Nose_2=NHparam%Q_Nose_2
xi=extcoord%xi
vxi=extcoord%vxi
NFkbT=syst%NFkbT
kbT=syst%kbT

E_NOSE_THERMOS=0.d0
do i=1, NT
  E_NOSE_THERMOS=E_NOSE_THERMOS+Q_Nose_2(i)*vxi(i)**2
enddo
E_NOSE_THERMOS=E_NOSE_THERMOS+NFkbT*xi(1)
do i=2, NT
  E_NOSE_THERMOS=E_NOSE_THERMOS+kbT*xi(i)
enddo
end function E_NOSE_THERMOS
!EF------------------------------------------------------------------

!BF------------------------------------------------------------------
! This function gives back the potential energy
!---------------------------------------------------------------------
function Epot(x,syst,pot)
use types;use assign_objects;use alloc_objects
use forcefieldPBD
use forcefieldNicolis
use forcefieldWCA
use forcefieldIOT
use forcefieldharm
implicit none
type(system_type),       intent(in)::syst
type(potential_type),    intent(in)::pot
double precision, intent(in)::x(syst%Npart,syst%dim)
double precision::Epot

 select case(pot%POTENTIAL)
    case("PBD")
      Epot=EpotPBD(x,syst,pot%potPBD)
    case("NICOLIS1D")
      Epot=EpotNICOLIS1D(x,syst,pot%potNICOLIS)
    case("NICOLIS2D")
      Epot=EpotNICOLIS2D(x,syst,pot%potNICOLIS)
    case("HARMOSC")
      Epot=0.5*pot%potHARMOSC%k*dot_product(x(1,:),x(1,:))
    case("DOUBLEWELL")
      Epot=pot%potDW%k4*x(1,1)**4+pot%potDW%k2*x(1,1)**2
    case("2DHARM")
      Epot=EP2DHARM(x(1,1),x(1,2),pot%pot2DHARM)
    case("NONE")
      Epot=0.d0
    case("WCA")
      Epot=EpotWCA(x,syst,pot%potWCA)
    case("IONTRANS")
      Epot=EpotIOT(x,syst,pot%potIOT)
    case default
      print *,"ERROR Epot: POTENTIAL=",pot%POTENTIAL
      stop
 end select
end function Epot
!EF--------------------------------------------------------------------

!BF-------------------------------------------------------------------------
! This function gives back the total energy
!---------------------------------------------------------------------
function Etotal(phase,syst,pot)
use types;use assign_objects;use alloc_objects
implicit none
type(system_type),       intent(in)::syst
type(potential_type),    intent(in)::pot
type(phasepoint_type), intent(in)::phase
double precision::Etotal
double precision::x(syst%Npart,syst%dim),v(syst%Npart,syst%dim),m_2(syst%Npart)
integer::N,d

x=phase%phasexv%x
v=phase%phasexv%v
m_2=syst%masses_2
N=syst%Npart
d=syst%dim

Etotal=Epot(x,syst,pot)+EKIN(v,m_2,N,d)

end function Etotal
!EF--------------------------------------------------------------------------

end module forcefield
!EM---------------------------------------
