!BM---------------------------------------------
module assign_objects 
use types
implicit none

!INTERFACE ASSIGMENTS/PROCEDURES: define operators on the objects(types)

  !BI------------------------------------------------
  interface assignment( = )
    module procedure assign_phasexv, assign_extc, assign_phase, &
                     assign_MDinout_param, assign_timeslice, assign_path, &
                     assign_elec
  end interface
  !EI-----------------------------------------------

contains
  !BS-----------------------------------------------------
  subroutine assign_phasexv(a,b)
  implicit none
  type(phasexv_type), intent(out)::a
  type(phasexv_type), intent(in)::b
    a%x=b%x
    a%v=b%v
  end subroutine assign_phasexv
  !ES-----------------------------------------------------

  !BS-----------------------------------------------------
  subroutine assign_extc(a,b)
  implicit none
  type(extended_coordinates_type), intent(out)::a
  type(extended_coordinates_type), intent(in)::b
    a%xi=b%xi
    a%vxi=b%vxi
  end subroutine assign_extc
  !ES-----------------------------------------------------

  !BS-----------------------------------------------------
  subroutine assign_elec(a,b)
  implicit none
  type(electrons_type), intent(out)::a
  type(electrons_type), intent(in)::b
    a%ESTRUCFILE=b%ESTRUCFILE
    a%EFILE_INDEX=b%EFILE_INDEX
    a%KIN=b%KIN
    a%POT=b%POT
    a%ETOT=b%ETOT
    a%EHAM=b%EHAM
    a%TEMP_inst=b%TEMP_inst
    a%Eelectrons=b%Eelectrons
    a%rev_velec=b%rev_velec
    a%WCENT=b%WCENT
  end subroutine assign_elec
  !ES-----------------------------------------------------


 !BS-----------------------------------------------------
  subroutine assign_phase(a,b)
  implicit none
  type(phasepoint_type), intent(out)::a
  type(phasepoint_type), intent(in)::b
    a%phasexv=b%phasexv
    a%extcoord=b%extcoord
    a%electrons=b%electrons
  end subroutine assign_phase
  !ES-----------------------------------------------------

  !BS-----------------------------------------------------
  subroutine assign_MDinout_param(a,b)
  implicit none
  type(MDinout_param_type), intent(out)::a
  type(MDinout_param_type), intent(in)::b
    a%F=b%F
    a%KIN=b%KIN
  end subroutine assign_MDinout_param
  !ES-----------------------------------------------------

  !BS-----------------------------------------------------
  subroutine assign_timeslice(a,b)
  implicit none
  type(timeslice_type), intent(out)::a
  type(timeslice_type), intent(in)::b
    a%phasepoint=b%phasepoint
    a%MDinout_param=b%MDinout_param
    a%OPS=b%OPS
  end subroutine assign_timeslice
  !ES-----------------------------------------------------

  !BS-----------------------------------------------------
  subroutine assign_path(a,b)
  implicit none
  type(path_type), intent(out)::a
  type(path_type), intent(in)::b
  integer::Lpath,i,NX
  Lpath=b%Lpath
  NX=b%NX

  !Technical note:
  !Following explicit loop is required
  !a%timeslices(1:Lpath)=b%timeslices(1:Lpath)
  !would yield a pointer assignment instead simply copying data
  do i=1,NX
    a%timeslices(i)=b%timeslices(i)   
  enddo

  a%Lpath=b%Lpath
  a%opmax=b%opmax
  a%opmin=b%opmin
  a%iopmax=b%iopmax
  a%iopmin=b%iopmin
  a%index_acc=b%index_acc
  a%index_shoot=b%index_shoot
  a%end=b%end
  a%start=b%start
  a%cross=b%cross
  a%ACCREJ=b%ACCREJ
  a%MCmove=b%MCmove
  a%OPshoot=b%OPshoot
  a%iOPshoot_old=b%iOPshoot_old
  a%iOPshoot_new=b%iOPshoot_new

  end subroutine assign_path
  !ES-----------------------------------------------------

end module assign_objects
!EM----------------------------------------------
