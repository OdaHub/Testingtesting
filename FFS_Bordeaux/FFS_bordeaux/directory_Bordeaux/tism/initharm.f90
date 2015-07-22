!BM--------------------------
Module initHARM

Contains

  !BS------------------------
  subroutine initializeHO
  use inputpar 
  use pot_module
  use system_module
  implicit none
  
  print *, "initialize Harmonic Oscilator potential"
 
  pot%potHARMOSC%k=kharm
  if ((Npart/=1).OR.(dim/=1)) then
    print *,"Harmonic oscillator test system:"
    print *,"Npart,dim should be set to 1!"
    stop
  endif
  
  end subroutine initializeHO
  !ES------------------------

  !BS------------------------
  subroutine initializeDW
  use inputpar
  use pot_module
  use system_module
  implicit none

  print *, "initialize DoubleWell  potential"

  pot%potDW%k4=doublewellk4
  pot%potDW%k2=doublewellk2
  if ((Npart/=1).OR.(dim/=1)) then
    print *,"DoubleWell system:"
    print *,"Npart,dim should be set to 1!"
    stop
  endif

  end subroutine initializeDW
  !ES------------------------



 !BS------------------------
  subroutine initializeHO2D
  use inputpar
  use pot_module
  use system_module
  implicit none

  print *, "initialize 2D Harmonic Oscilator potential"

  pot%pot2DHARM%kx=kxharm
  pot%pot2DHARM%ky=kyharm
  pot%pot2DHARM%x=xharm
  pot%pot2DHARM%y=yharm
  pot%pot2DHARM%ycut=ycut

  if ((Npart/=1).OR.(dim/=2)) then
    print *,"2D Harmonic oscillator test system:"
    print *,"Npart should be 1, dim should be set to 2!"
    stop
  endif

  end subroutine initializeHO2D
  !ES------------------------


 !BS------------------------
  subroutine initializeNICOLIS
  use inputpar
  use system_module
  use pot_module
  implicit none

  print *, "initialize NICOLIS-potential"
  

  pot%potNICOLIS%mu=MU_NIC
  pot%potNICOLIS%lam=lam_NIC
  pot%potNICOLIS%gam=gam_NIC
  pot%potNICOLIS%massgamma=mass*gamma
  pot%potNICOLIS%Q=sqrt(2.d0/(pot%potNICOLIS%massgamma*syst%beta))
  print *, "Q=",pot%potNICOLIS%Q 

  end subroutine initializeNICOLIS
  !ES------------------------



end Module initHARM
!EM--------------------------
