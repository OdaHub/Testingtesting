!BM--------------------------
Module initPBD

Contains

  !BS------------------------
  subroutine initializePBD
  use inputpar 
  use pot_module
  use system_module
  implicit none
  integer::i
  double precision::a(Npart),d(Npart),da2(Npart),AL_2,S_2

  print *,"initialize PBD potential"

  if (dim/=1) then
    print *,"PBD system: dim should be set to 1!"
    stop
  endif
 
  AL_2= 0.5d0*ALPHA
  S_2=0.5d0*S

  do i=1,Npart
    if( (seq(i)=='A').or.(seq(i)=='T')) then
      d(i) = dAT
      a(i) = aAT
    else if((seq(i)=='G').or.(seq(i)=='C')) then
      d(i) = dGC
      a(i) = aGC
    else
      write(*,*)' Sequence must be A, T, G, or C'
      stop
    endif
    da2(i)=2*a(i)*d(i)
  enddo

  !Save values
  pot%potPBD%DAT=DAT
  pot%potPBD%DGC=DGC
  pot%potPBD%AAT=AAT
  pot%potPBD%AGC=AGC
  pot%potPBD%S=S
  pot%potPBD%S_2=S_2
  pot%potPBD%RHO=RHO
  pot%potPBD%ALPHA=ALPHA
  pot%potPBD%AL_2=AL_2

  allocate(pot%potPBD%seq(Npart))
  allocate(pot%potPBD%a(Npart))
  allocate(pot%potPBD%d(Npart))
  allocate(pot%potPBD%da2(Npart))

  pot%potPBD%seq(1:Npart)=seq(1:Npart)
  pot%potPBD%a(1:Npart)=a(1:Npart)
  pot%potPBD%d(1:Npart)=d(1:Npart)
  pot%potPBD%da2(1:Npart)=da2(1:Npart)

  pot%potPBD%BIAS=BIAS
  pot%potPBD%BIASEXP=BIASEXP
  pot%potPBD%BIASPREFAC= BIASPREFAC
  pot%potPBD%BIASCUTOFF= BIASCUTOFF

  
  end subroutine initializePBD
  !ES------------------------


end Module initPBD
!EM--------------------------
