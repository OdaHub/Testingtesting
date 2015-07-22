!BM--------------------------
Module initiot

Contains

  !BS------------------------
  subroutine initializeIONTRANS
  use inputpar 
  use pot_module
  use forcefieldiot
  implicit none

  print *,"initialize IONTRANS potential"

  if (dim/=2) then
    print *,"IONTRANS system: dim should be set to 2!"
    stop
  endif
!  if (Npart<3) then
!    print *,"IONTRANS system: Npart should be >2"
!    stop
!  endif


  !Save values
  pot%potIOT%eps=IOT_EPSILON
  pot%potIOT%SIG=IOT_SIGMA
  pot%potIOT%eps2=IOT_EPS2
  pot%potIOT%SIG2=IOT_SIG2

  pot%potIOT%r0=(2**(1.d0/6.d0))*IOT_SIGMA
  pot%potIOT%r02=(2**(1.d0/6.d0))*IOT_SIG2
  pot%potIOT%Vshift=LJ(iot_eps2,iot_sig2,pot%potIOT%r0)

  pot%potIOT%COOPERATIVITY=IOT_COOPERATIVITY
  pot%potIOT%FCOOP=IOT_FCOOP
  pot%potIOT%ND=IOT_ND      
  pot%potIOT%NN=IOT_NN
  pot%potIOT%NC =IOT_NC 
  pot%potIOT%RCOOP=IOT_RCOOP
  
  
  end subroutine initializeIONTRANS
  !ES------------------------


end Module initiot
!EM--------------------------
