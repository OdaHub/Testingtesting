!BM---------------------------------------
module forcefieldNicolis

contains

 !BF------------------------------------------
  function forceNicolis1D(x,syst,Nicolisparam)
  use types;use assign_objects;use alloc_objects
  implicit none
  type(system_type),       intent(in)::syst
  type(potNicolis_type),       intent(in)::Nicolisparam
  double precision, intent(in)::x(syst%Npart,1)
  double precision::forceNicolis1D(syst%Npart,1)


  forceNicolis1D(1,1)=Nicolisparam%lam*x(1,1)-&
                      Nicolisparam%mu* x(1,1)**3-x(1,1)**5
  forceNicolis1D=forceNicolis1D*Nicolisparam%massgamma


  end function forceNicolis1D
  !EF------------------------------------------------------
 
   !BF------------------------------------------
  function forceNicolis2D(x2D,syst,Nicolisparam)
  use types;use assign_objects;use alloc_objects
  implicit none
  type(system_type),       intent(in)::syst
  type(potNicolis_type),       intent(in)::Nicolisparam
  double precision, intent(in)::x2D(syst%Npart,2)
  double precision::forceNicolis2D(syst%Npart,2)
  double precision::x,y,lam,mu,gam

  x=x2D(1,1);y=x2D(1,2)
  lam=Nicolisparam%lam;mu=Nicolisparam%mu;gam=Nicolisparam%gam


  forceNicolis2D(1,1)=-lam*x-mu*x**2+gam*y**2-x**3
  forceNicolis2D(1,2)=-lam*y+2*gam*x*y-y**3
  forceNicolis2D=forceNicolis2D*Nicolisparam%massgamma

  end function forceNicolis2D
  !EF------------------------------------------------------

  !BF------------------------------------------
  function EPotNicolis1D(x,syst,Nicolisparam)
  use types;use assign_objects;use alloc_objects
  implicit none
  type(system_type),       intent(in)::syst
  type(potNicolis_type),       intent(in)::Nicolisparam
  double precision, intent(in)::x(syst%Npart,1)
  double precision::EpotNicolis1D

  EpotNicolis1D=-Nicolisparam%lam*x(1,1)**2/2.d0+&
                 Nicolisparam%mu* x(1,1)**4/4.d0+&
                                  x(1,1)**6/6.d0
  EpotNicolis1D=EpotNicolis1D*Nicolisparam%massgamma

  end function EpotNicolis1D
  !EF-----------------------------------------------------

   !BF------------------------------------------
  function EPotNicolis2D(x2D,syst,Nicolisparam)
  use types;use assign_objects;use alloc_objects
  implicit none
  type(system_type),       intent(in)::syst
  type(potNicolis_type),       intent(in)::Nicolisparam
  double precision, intent(in)::x2D(syst%Npart,2)
  double precision::EpotNicolis2D
  double precision::x,y,lam,mu,gam
  
  x=x2D(1,1);y=x2D(1,2)
  lam=Nicolisparam%lam;mu=Nicolisparam%mu;gam=Nicolisparam%gam
  EpotNicolis2D=lam*(x**2+y**2)/2.d0+mu*x**3/3.d0-gam*x*y**2+(x**4+y**4)/4.d0
  EpotNicolis2D=EpotNicolis2D*Nicolisparam%massgamma

  end function EpotNicolis2D
  !EF------------------------------------------------------

end module forcefieldNicolis
!EM---------------------------------------

