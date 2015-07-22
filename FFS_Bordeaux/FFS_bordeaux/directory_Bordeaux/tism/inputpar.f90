!BM--------------------------------------------------
module inputpar
use stringlengths
implicit none
  integer::Npart,dim
  double precision::Temp,mass,dt
  logical::MASSPOLYDISP
  double precision, allocatable::masses(:)

  character(LEN=MSTR)::REACTION_COORDINATE
  integer:: REACDIST_ATNR(2)


  character(LEN=MSTR)::TASK,POTENTIAL,DYNAMICS,PATHINFO,STARTCONDITION,&
                       TRANSALGORITHM
  character(LEN=LSTR)::RESTARTFILE
  character(LEN=LSTR)::EXTERNAL_PROGRAM
  character(LEN=MSTR)::CPMD_UNITS

  integer::Ncyc,SKIPE,SKIPO,SKIPT,skipP
  logical::RESTART,PBC,RENUM,NOPRINT,SKIPINITOPTWAVE,SKIPINITMDSTEP

  double precision::DAT,DGC,AAT,AGC,S,RHO,ALPHA
  character, allocatable::seq(:)

  double precision::gamma
  logical::HIGH_FRICTION_LIMIT
  double precision::AN_FREQ
  integer::Nthermo
  double precision::kharm,kxharm,kyharm,xharm,yharm,ycut
  double precision::DOUBLEWELLK4,DOUBLEWELLK2

  double precision::INTERFACEL,INTERFACEM,INTERFACER
  integer::NOPS,NX
  double precision::timerevfreq,sigdp,swapfreq

  integer::NUMINT_PPS
  logical::RELATIVESHOOTS,NULLMOVES,SWAPSIMUL
  double precision, allocatable::PPS_INTERFACES(:),PPS_SIGDP(:),&
                                 RELATIVE_SHOOTFREQ(:)
  integer::NTRANSRUN

  double precision:: LAM_NIC,MU_NIC,GAM_NIC
  logical::READINITIAL_POSITIONS
  double precision, allocatable::initpos(:,:)
  logical::TWO_POINT_METHOD
  integer::STARTRANDOM
  double precision, allocatable::range_left(:),range_right(:),intstep(:)

  logical::BIAS
  double precision::BIASEXP,BIASPREFAC,BIASCUTOFF
  
  integer::NCPU,TWAIT,NSUBCYCLES,MAXWAIT
  character*2, allocatable::ATOM_TYPES(:)
  character(LEN=LSTR)::MASS_INPUT_SPECIFICATION

  !!logical::WANNIER_CENTERS
  integer::NWANNIER,NCYCLE_RESTART
  double precision::BOXLENGTH
  logical::RESUBMIT
  
  double precision::WCA_EPSILON,WCA_SIGMA,WCA_H,WCA_W 
  double precision::IOT_EPSILON,IOT_SIGMA,IOT_EPS2,IOT_SIG2
  double precision::IOT_FCOOP,IOT_ND,IOT_NN,IOT_NC,IOT_RCOOP
  logical::IOT_COOPERATIVITY

  double precision::NVE_ENERGY
  logical::RESCALE_ENERGY,SET_LINMOM_ZERO
  logical::FORCEFIELDMATCHING,BIASV

end module inputpar
!EM---------------------------------------------------
