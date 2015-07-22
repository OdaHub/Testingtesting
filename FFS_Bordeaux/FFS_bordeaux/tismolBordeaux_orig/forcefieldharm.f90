!BM---------------------------------------
module forcefieldharm

contains



  !BF-------------------------------------------
  function EP2DHARM(x,y,pot)
  use types;use assign_objects;use alloc_objects
  implicit none
  type(pot2DHARM_type),       intent(in)::pot
  double precision, intent(in)::x,y
  double precision::Ep2DHARM
  double precision::kx,ky,px,py,c

  kx=pot%kx
  ky=pot%ky
  px=pot%x
  py=pot%y
  c=pot%ycut


  EP2DHARM=.5d0*(kx*(x-px)**2+ky*(y-py)**2)
  if ((y-py)>c) EP2DHARM=0.5d0*ky*(c-py)**2

  end function Ep2DHARM
  !EF-------------------------------------------

  !BF-------------------------------------------
  function force2dHARM(x,y,pot)
  use types;use assign_objects;use alloc_objects
  implicit none
  type(pot2DHARM_type),       intent(in)::pot
  double precision, intent(in)::x,y
  double precision::force2DHARM(1,1:2)
  double precision::kx,ky,px,py,c

  kx=pot%kx
  ky=pot%ky
  px=pot%x
  py=pot%y
  c=pot%ycut

  force2DHARM(1,1)=-kx*(x-px)
  force2DHARM(1,2)=-ky*(y-py)
  if ((y-py)>c) force2DHARM=0.d0 

  end function force2dHARM
  !EF-------------------------------------------



end module forcefieldharm
!EM---------------------------------------

