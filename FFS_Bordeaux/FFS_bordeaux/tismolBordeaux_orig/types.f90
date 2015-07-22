!----------------------------------------------------------------
!     This module defines a set of types(=objects) such as
!     phasepoint, timeslices, path, and path-ensembles.
!     Technical note:
!     -------------------------------------------
!     These types can contain arrays that require dynamic array allocation.
!     As allocatable arrays cannot be part of types, these have been
!     defined as pointer variables. However, directly Obj1=Obj2 will result 
!     that the pointers inside the objects will make a pointer assignment
!     (Then changing a pointer in Obj2 will automatically change the pointer in
!     Obj1. Hence Obj1 and Obj2 are dynamically linked). This is usually not
!     the required behavior. We would like to have that Obj1=Obj2 simply
!     copies all data from Obj2 in Obj1 and if we change Obj1 (Obj2) in a later
!     stage, this will not influence Obj2 (Obj1).  To obtain this behavior,
!     we have made the interface  assignment( = ). It is important to realize 
!     that adding component to the objects requires editing this interface!
!     To make it easier to (de)allocate all pointers inside an object, we 
!     have created the interface alloc and dealloc.
!     -----------------------------------------------
!
!BM---------------------------------------------
module types
  use stringlengths
  implicit none


  !TYPE DEFINITIONS

  !BT-----------------------------------------------
  type phasexv_type 
  !Note that changing this type content requires editing
  !assign_objects.f90 and alloc_objects.f90
    integer::N,d
    double precision, pointer, dimension(:,:)::x,v
  end type phasexv_type
  !ET-----------------------------------------------
 
  !BT-----------------------------------------------
  type extended_coordinates_type
  !Note that changing this type content requires editing
  !assign_objects.f90 and alloc_objects.f90
    integer::NT
    double precision, pointer, dimension(:)::xi,vxi
  end type extended_coordinates_type
  !ET-----------------------------------------------

  !BT-----------------------------------------------
  type electrons_type
    !Note that changing this type content requires editing
    !assign_objects.f90 and alloc_objects.f90
    integer::NWANNIER,dwc
    character(LEN=LSTR)::ESTRUCFILE
    integer::EFILE_INDEX
    double precision::KIN,POT,ETOT,EHAM,TEMP_inst,Eelectrons
    logical::rev_velec
    double precision, pointer, dimension(:,:)::WCENT
  end type electrons_type
  !ET-----------------------------------------------


  !BT-----------------------------------------------
  type phasepoint_type
  !Note that changing this type content requires editing
  !assign_objects.f90 and alloc_objects.f90
    integer::N,d,NT,NWANNIER,dwc
    type(phasexv_type)::phasexv
    type(extended_coordinates_type)::extcoord
    type(electrons_type)::electrons
  end type phasepoint_type
  !ET-----------------------------------------------

  !BT------------------------------------------------------
  type MDinout_param_type
  !Note that changing this type content requires editing
  !assign_objects.f90 and alloc_objects.f90
    integer::N,d
    double precision, pointer, dimension(:,:)::F
    double precision::KIN
  end type MDinout_param_type
  !ET------------------------------------------------------

  !BT-----------------------------------------------------
  type timeslice_type
  !Note that changing this type content requires editing
  !assign_objects.f90 and alloc_objects.f90
    integer::N,d,NT,NOPS,NWANNIER,dwc
    type(phasepoint_type)::phasepoint
    type(MDinout_param_type)::MDinout_param
    double precision, pointer, dimension(:)::OPS
  end type timeslice_type
  !ET-----------------------------------------------------

  !BT-----------------------------------------------------
  type path_type
  !Note that changing this type content requires editing
  !assign_objects.f90 and alloc_objects.f90
    integer::NX,N,d,NT,NOPS,NWANNIER,dwc
    type(timeslice_type),pointer, dimension(:)::timeslices
    integer::Lpath
    double precision::opmin,opmax
    integer::iopmin,iopmax
    integer::index_acc,index_shoot
    character::start,end,cross
    character*3::ACCREJ,MCmove
    double precision:: OPshoot
    integer::iOPshoot_old,iOPshoot_new
    
  end type path_type
  !ET-----------------------------------------------------

  !BT-----------------------------------------------
  type timestep_type 
    double precision::dt,dt_2,dt_4,dt_8
    double precision, pointer, dimension(:)::dt_2m,dt_m
    double precision::Bigdt_ru
  end type timestep_type
  !ET-----------------------------------------------

  !BT-----------------------------------------------
  type potPBD_type
    double precision::DAT,DGC,AAT,AGC,S,RHO,ALPHA,AL_2,S_2
    character, pointer, dimension(:)::seq
    double precision, pointer, dimension(:)::a,d,da2
    double precision::surface
    logical::BIAS
    double precision::BIASEXP,BIASPREFAC,BIASCUTOFF
  end type potPBD_type
  !ET-----------------------------------------------

  !BT-----------------------------------------------
  type potNICOLIS_type
    double precision::LAM,MU,GAM,Q,massgamma
  end type potNICOLIS_type
  !ET-----------------------------------------------

  !BT-----------------------------------------------
  type potHARMOSC_type
    double precision::k
  end type potHARMOSC_type
  !ET-----------------------------------------------

 !BT-----------------------------------------------
  type potDW_type
    double precision::k4,k2
  end type potDW_type
  !ET-----------------------------------------------



  !BT-----------------------------------------------
  type pot2DHARM_type
    double precision::kx,ky,x,y,ycut
  end type pot2DHARM_type
  !ET-----------------------------------------------


  !BT-----------------------------------------------
  type external_type
    double precision::a0,atomicv_redv,inv_atomicv_redv,inv_a0,hartree_eV
  end type external_type
  !ET-----------------------------------------------

  !BT-----------------------------------------------
  type potWCA_type
    double precision::epsilon,sigma,h,w,r0
  end type potWCA_type
  !ET-----------------------------------------------

  !BT-----------------------------------------------
  type potIOT_type
    double precision::eps,sig,eps2,sig2,r0,r02,Vshift
    logical::cooperativity
    double precision::FCOOP,ND,NN,NC,RCOOP
  end type potIOT_type
  !ET-----------------------------------------------

  
  !BT-----------------------------------------------
  type potential_type
    character(LEN=MSTR)::POTENTIAL
    type(potPBD_type)::potPBD
    type(potHARMOSC_type)::potHARMOSC
    type(potDW_type)::potDW
    type(pot2DHARM_type)::pot2DHARM
    type(potNICOLIS_type)::potNICOLIS
    type(external_type)::external
    type(potWCA_type)::potWCA
    type(potIOT_type)::potIOT
  end type potential_type
  !ET-----------------------------------------------

  !BT-----------------------------------------------
  type langevin_type
    double precision::gamma
    double precision, pointer, dimension(:)::s12os11,sqrts11,sqrtSos11
    double precision, pointer, dimension(:)::a2,b1,b2
    double precision::c0,a1
    logical::HIGH_FRICTION_LIMIT
    double precision, pointer, dimension(:)::sigma_lange,bDdt
  end type langevin_type
  !ET-----------------------------------------------

  !BT-----------------------------------------------
  type Andersen_type
    double precision::AN_FREQ
    double precision::Freqdt
  end type Andersen_type
  !ET-----------------------------------------------

  !BT-----------------------------------------------
  type Nose_Hoover_type
  double precision, pointer, dimension(:)::Q_Nose,inv_Q_Nose,Q_Nose_2
  integer::Nthermo
  end type Nose_Hoover_type
  !ET-----------------------------------------------

  !BT-----------------------------------------------
  type NVE_type
  double precision::NVE_energy
  logical::RESCALE_ENERGY,SET_LINMOM_ZERO
  end type NVE_type
  !ET-----------------------------------------------


  !BT-----------------------------------------------
  type external_dynamics_type
  !!logical::RENAME_ESTRUCTURE_FILE
  character*4, pointer, dimension(:)::EXT
  !!character(LEN=LSTR)::ESTRUCFILE
  end type external_dynamics_type
  !ET-----------------------------------------------


  !BT-----------------------------------------------
  type dynamics_type
    character(LEN=MSTR)::DYNAMICS
    type(NVE_type)::NVE
    type(langevin_type)::langevin
    type(Andersen_type)::Andersen
    type(Nose_Hoover_type)::Nose_Hoover
    type(external_dynamics_type)::ext_dyn
  end type dynamics_type
  !ET-----------------------------------------------
  
  !BT-------------------------------------------------
  type system_type
    integer::Npart
    double precision, pointer, dimension(:)::masses,masses_2,sigma_v
    character*2, pointer, dimension(:)::ATOM_TYPES
    double precision::Temp
    double precision::kb, kbT,beta,x2_kbNF,x2_kbNFTRANS,NFkbT,kbNF_2,kbNFTRANS_2
    double precision::a0,atomicv_redv,inv_a0,inv_atomicv_redv,hartree_eV
    integer::dim,NF,NFTRANS
    logical::PBC,NOPRINT
    character(LEN=LSTR)::EXTERNAL_PROGRAM,CRASHFILE
    character(LEN=MSTR)::REACTION_COORDINATE
    integer:: REACDIST_ATNR(2)
    integer::NCPU,TWAIT,IUENCPMD,NSUBCYCLES,IUWANCPMD,MAXWAIT,IUCRASH
    integer::NWANNIER,dwc
    integer::NCPMD_UNSAVED,NCYCLE_RESTART,icrash,NX
    double precision::BOXLENGTH
    logical::RESUBMIT
  end type system_type
  !ET-------------------------------------------------

  !BT-------------------------------------------------
  type output_type
    character(len=LSTR)::ENERGY_FILE,TRAJECTORY_FILE,ORDER_PARAM_FILE, &
                         PATH_FILE,CROSSFILE, &
                         RESTARTFILE_READ,RESTARTFILE_WRITE,TRANSFILE, &
                         input_CPMD_optwav, output_CPMD_optwav,&
                         input_CPMD_MD,output_CPMD_MD, &
                         input_CPMD_TIS, output_CPMD_TIS, &
                         input_CPMD_optwav_HEAD,  input_CPMD_MD_HEAD, &
                         input_CPMD_TIS_HEAD, input_CPMD_TAIL, &
                         GEOMETRY00,ESTRUC_a,CROSSPOINTS,FFSFILE


    character(len=MSTR)::PATHINFO
    integer:: IUEN, IUTRAJ, IUOP,IUPATH,IUCR,IURES,IUTRANS,IUCROSSPOINTS,IUFFS
    integer:: skipE,skipT,skipO,skipP
    double precision,pointer, dimension(:)::crossing_plane
    integer:: NCROSSPLANES  
    integer:: icyclestart
  end type output_type
  !ET-------------------------------------------------

  !BT-----------------------------------------------------
  type TIS_type
    integer::NX,NOPS  
    double precision::INTERFACEL,INTERFACEM,INTERFACER   
    double precision::timerevfreq
    double precision,pointer, dimension(:)::sigdp_sqrtm
    logical::aimless,BIASV
    integer::N
    character::startcondition
  end type TIS_type
  !ET-----------------------------------------------------

  !BT-----------------------------------------------------
  type trans_type
    integer::NX,NOPS,NTRANSRUN  
    double precision::INTERFACEL,INTERFACEM,INTERFACER
    character(LEN=MSTR):: TRANSALGORITHM
    logical::two_point_method
  end type trans_type
  !ET-----------------------------------------------------

!BT-----------------------------------------------------
  type ffs_type
    integer::NX,NOPS
    double precision::INTERFACEL,INTERFACEM,INTERFACER
    double precision::startpoints(100000,2)
    integer::nstartpoints
  end type ffs_type
  !ET-----------------------------------------------------

  
  !BT--------------------------------------------------------
  type path_ensemble
    type(path_type), pointer, dimension(:)::PATHS
    double precision,pointer,dimension(:)::PPS_interfaces
    double precision,pointer,dimension(:,:)::PPS_sigdp_sqrtm
    logical,pointer,dimension(:)::PPS_aimless
    double precision,pointer,dimension(:)::relative_shootfreq(:)
    double precision::timerevfreq,swapfreq
    integer, pointer,dimension(:)::IUPATH,IUEN,IUTRAJ,IUOP,IUCR,IUCROSSPOINTS
    logical::RELATIVESHOOTS,nullmoves,SWAPSIMUL,FORCEFIELDMATCHING
    character*3, pointer, dimension(:)::EXT


    integer::N,d,NT,NOPS,NX,NUMINT,NWANNIER,dwc
  end type path_ensemble
  !ET------------------------------------------------------

  !BT---------------------------------------------------------
  type numinteg_type
    double precision,pointer,dimension(:)::RANGE_LEFT,RANGE_RIGHT,INTSTEP
    double precision,pointer,dimension(:,:)::FREE_EN_ARRAY,prob_array
    integer,pointer,dimension(:)::npoints
    integer::dim,maxpoints
    double precision::phaseVol
  end type numinteg_type
  !ET------------------------------------------------------------


end module types
!EM----------------------------------------------
