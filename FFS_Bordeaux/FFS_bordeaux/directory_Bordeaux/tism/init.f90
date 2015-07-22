!BM--------------------------
Module init

Contains

  !BS------------------------
  subroutine initialize
  use stringlengths
  use inputpar
  use system_module
  use timestep_module
  use phase_module
  use pot_module
  use dyn_module
  !!use initPBD
  !!use initWCA
  !!use initHARM
  use initpot
  use initLange
  use initNoseHoover
  use initdynext
  use initMD
  use initTIS
  use initPPS
  use initNUMINTEG
  use output_module
  use atoms
  use random
  implicit none
  double precision::evJ,kb,kbT,beta,sigma_v(Npart),masses_2(Npart)
  double precision::x2_kbNF,NFkbT,ran,x2_kbNFTRANS,kbNF_2,kbNFTRANS_2
  integer:: i
  integer::NF,NFTRANS

  !NF: degrees of freedom for system without momentum condervation
  !NFTRANS: degrees of freedom with substraction of translations

  do i=1,STARTRANDOM
   !!call random_number(ran)
   ran=random01()
   print *,"init RANDOM_NUMBER:",ran
  enddo

  if ((POTENTIAL=="WCA").OR.(POTENTIAL=="IONTRANS").OR.(POTENTIAL=="DOUBLEWELL")) then
    print *,"USING REDUCED UNITS kb=1"
    kb=1.d0
  else
    evJ=1.60219D-19 !eV in J
    kb=1.3806503D-23     ! J/K
    kb=kb/evJ
    print *,"kb=", kb, " eV/K"
  endif
  
  kbT=kb*Temp
  beta=1.d0/kbT                !inverse temperature
  
  allocate(syst%masses(Npart))
  allocate(syst%masses_2(Npart))
  allocate(syst%sigma_v(Npart))
  allocate(syst%atom_types(Npart))

  if (.NOT.MASSPOLYDISP) then
    allocate(masses(Npart))
    masses=mass !assuming uniform mass!
  else if (MASS_INPUT_SPECIFICATION=="ATOM_TYPES") then
    allocate(masses(Npart))
    call set_atomic_masses(masses,atom_types,Npart) 
  endif
  syst%atom_types=atom_types
  masses_2=0.5d0*masses
  sigma_v=sqrt(kbT/masses)   !width of gaussian velocity distr.
  !save values in structure:syst
  syst%masses=masses
  syst%masses_2=masses_2
  syst%Npart=Npart
  syst%dim=dim
  syst%Temp=Temp
  syst%kb=kb
  syst%kbT=kbT
  syst%beta=beta 
  syst%sigma_v=sigma_v
  syst%PBC=PBC
  syst%NCYCLE_RESTART=NCYCLE_RESTART
  syst%NX=NX
  syst%BOXLENGTH=BOXLENGTH
  
  NF=Npart*dim
  NFTRANS=NF-dim
  NFkbT=NF*kbT

  syst%NF=NF
  syst%NFTRANS=NFTRANS
  syst%NFkbT=NFkbT

  x2_kbNF=2.d0/(kb*NF)
  x2_kbNFTRANS=2.d0/(kb*NFTRANS)
  kbNF_2=1.d0/x2_kbNF
  kbNFTRANS_2=1.d0/x2_kbNFTRANS
  syst%x2_kbNF=x2_kbNF
  syst%x2_kbNFTRANS=x2_kbNFTRANS
  syst%kbNF_2=kbNF_2
  syst%kbNFTRANS_2=kbNFTRANS_2

  syst%REACTION_COORDINATE=REACTION_COORDINATE
  syst%REACDIST_ATNR=REACDIST_ATNR

  syst%noprint=noprint
 
  allocate(timestep%dt_m(Npart)) 
  allocate(timestep%dt_2m(Npart))
  timestep%dt=dt
  timestep%dt_2=dt/2
  timestep%dt_4=dt/4
  timestep%dt_8=dt/8
  timestep%dt_m=dt/masses
  timestep%dt_2m=timestep%dt_m/2
  timestep%Bigdt_ru=dt           !default Bigdt_ru=dt. For CPMD Bigdt=Nsubcycles*dt*unit-conversion-factor

  
  dyn%DYNAMICS=DYNAMICS
  pot%POTENTIAL=POTENTIAL

  output%RESTARTFILE_READ=RESTARTFILE
  output%RESTARTFILE_WRITE="TISMOL.RESTART"
  output%IURES=3000
 
  call initpotential

  select case(DYNAMICS)
    case("NVE")
      Nthermo=0 
      NWANNIER=0
      dyn%NVE%NVE_ENERGY=NVE_ENERGY
      dyn%NVE%RESCALE_ENERGY=RESCALE_ENERGY
      dyn%NVE%SET_LINMOM_ZERO=SET_LINMOM_ZERO
    case("LANGEVIN")
      call init_Langevin
    case("ANDERSEN")
      Nthermo=0
      NWANNIER=0
      dyn%Andersen%AN_FREQ=AN_FREQ
      dyn%Andersen%freqdt=AN_freq*dt
    case("NOSEHOOVER")
       call init_NoseHoover 
    case("EXTERNAL")
       call init_dynamics_external 
    case("RANDOMWALK")
      Nthermo=0
      NWANNIER=0
    case default
      print *,"ERROR INITIALIZE: NO DYNAMICS: ",DYNAMICS
      stop
  end select

  select case(TASK)
    case("MD")
      call init_MD
    case("TIS")
      call init_TIS
    case("PPTIS")
      call init_TIS
    case ("PPS")
      call init_PPS
    case("TRANSMISSION")
      call init_trans
    case("FFS")
      call init_ffs
    case("NUMERICINTEG")
      call init_numericinteg
    case default
      print *,"ERROR INITIALIZE TASK=",TASK
  end select

  end subroutine initialize
  !ES------------------------

end Module init
!EM--------------------------
