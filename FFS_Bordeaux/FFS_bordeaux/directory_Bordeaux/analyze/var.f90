!BM----------------------------------------------------
module stringlengths
  integer, parameter::SSTR=10
  integer, parameter::MSTR=50
  integer, parameter::LSTR=100
  integer, parameter::XLSTR=200
end module stringlengths
!EM----------------------------------------------------

!BM---------------------------------------------------
module defaultfile
  use stringlengths
  character(len=LSTR), parameter::inputdefault="TISMOL.DEFAULT"
end module defaultfile
!EM---------------------------------------------------


!BM---------------------------------------------
Module var
use stringlengths
implicit none
integer::PSF_SIM,NPSFMAX,NPSF,NGRID,Npart,Nwannier,MAXBLOCKLENGTH,DIM,BLOCKSKIP
integer::MSDLENGTH
character(LEN=LSTR), allocatable::psfiles(:)
character(LEN=XLSTR)::dir
double precision::pi
double precision::dt,dt_ut,Temp,kb,kbT,beta,mass,gamma_invut, mass_eVns2_A2
double precision::timeunit

character(LEN=MSTR)::POTENTIAL,TASK,DYNAMICS
integer::MDcycles,SKIPcross !,SKIPO,SKIPT
integer, allocatable::TIScycles(:)
!logical::RESTART,PBC

double precision::DAT,DGC,AAT,AGC,S,RHO,ALPHA
double precision::WCA_EPSILON,WCA_SIGMA,WCA_H,WCA_W 

double precision::kharm,DOUBLEWELLK4,DOUBLEWELLK2,DXINTEG
character*3::EXT(0:999)
!character, allocatable::seq(:)

double precision::gamma
integer::NUMINT
double precision,allocatable::INTFL(:),INTFM(:),INTFR(:),PPS_INTERFACES(:)
double precision::MD_INTFL,MD_INTFM,MD_INTFR
double precision::trans_INTFL,trans_INTFM,trans_INTFR,LBOX
logical:: TWO_POINT_METHOD,Brownian,HIGH_FRICTION_LIMIT
integer, allocatable::PMAXBLOCKLENGTH(:),PBLOCKSKIP(:)
integer::NMOVMAX
character*50::chartu,chareu,charlu,chartempu
logical::forcefieldmatching
end module var
!EM---------------------------------------------------

!BM----------------------------------------------------
Module results
double precision::flux1,flux2,flux3
double precision::relerrflux1,relerrflux2,relerrflux3
!double precision::corflux1,corflux2,corflux3
integer::NSA,NSB,NOSA,NOSB
integer::NCROSS(3),NEFFCROSS(3)
double precision, allocatable::crossprobab(:),&
           relerrcr(:),corcr(:),Lacc(:),Ltr(:), &
           accRate(:)
double precision::T0,T1,reT0,reT1,NCORT0,NCORT1,LTRT0,LTRT1,ACCRATET0,ACCRATET1
integer::NL0,NL1
double precision::kappa,kappaTST,rekap,rekapTST,ncorkap,ncorkapTST
double precision::prob1,prob2

end module results
!EM---------------------------------------------------
