!BM--------------------------
Module initNoseHoover

Contains

  !BS------------------------
  subroutine init_NoseHoover
  use inputpar 
  use system_module
  use dyn_module
  use phase_module
  implicit none
  double precision::ommax2 

  print *,"initialize Nose-Hoover parameters"

  NWANNIER=0  
  dyn%Nose_hoover%Nthermo=Nthermo
  allocate(dyn%Nose_hoover%Q_Nose(NTHERMO))
  allocate(dyn%Nose_hoover%inv_Q_Nose(NTHERMO))
  allocate(dyn%Nose_hoover%Q_Nose_2(NTHERMO))

  startpoint%extcoord%xi=0.d0;startpoint%extcoord%vxi=0.d0
  !optimization for setting Nose-Hoover chain masses
  !First the highest harmonic frequency in the system
  !Then Q(1)=(Nf kb T)/omega**2 and others Q(i)=(kb T)/omega**2
  !With Nf=dN the degrees of freedom
  !See Martyna et al Mol Phys. 87, 1117 (1996)
  select case (POTENTIAL)
    case("PBD")
      ommax2 = (2.0d0*dGC*aGC**2 + 4.0d0*s)/mass
    case("HARMOSC")
      ommax2=kharm/mass
    case default
      print *,"ERROR init_NoseHoover POTENTIAL=",POTENTIAL
      stop
    end select 
    dyn%Nose_hoover%Q_NOSE(1:Nthermo)=(syst%kbT)/ ommax2
    dyn%Nose_hoover%Q_NOSE(1)= dyn%Nose_hoover%Q_nose(1)*Npart*syst%dim
    dyn%Nose_hoover%inv_Q_Nose=1.d0/dyn%Nose_hoover%Q_Nose
    dyn%Nose_hoover%Q_Nose_2=0.5d0*dyn%Nose_hoover%Q_Nose
  end subroutine init_NoseHoover
  !ES------------------------


end Module initNoseHoover
!EM--------------------------
