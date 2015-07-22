!BM--------------------------
Module initWCA

Contains

  !BS------------------------
  subroutine initializeWCA
  use inputpar 
  use pot_module
  implicit none

  print *,"initialize WCA potential"

  if (dim/=2) then
    print *,"WCA system: dim should be set to 2!"
    stop
  endif

  !Save values
  pot%potWCA%epsilon=WCA_EPSILON
  pot%potWCA%SIGMA=WCA_SIGMA
  pot%potWCA%H=WCA_H
  pot%potWCA%W=WCA_W
  pot%potWCA%r0=(2**(1.d0/6.d0))*WCA_SIGMA

  
  end subroutine initializeWCA
  !ES------------------------


end Module initWCA
!EM--------------------------
