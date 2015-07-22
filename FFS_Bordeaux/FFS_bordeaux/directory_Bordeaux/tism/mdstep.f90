!BM------------------------
module mdstep

contains

!--------------------------------------------------------------------
! This subroutine performs a single MD timestep. The input phasepoint
! is overwritten for the new one. The MDinout_param parameter can be  
! either the force (in most cases) or the kinetic energy (Nose-Hoover) or
! nothing (external). It should be given once at the start of a MD loop. 
! The new MDinout_param is automatically calculated in this routine, so can 
! be used for the next  call which saves cpu time.
!-----------------------------------------------------------------------------
subroutine make_mdstep(phasepoint,syst,pot,dyn,MDinout_param,timestep)
use types;use assign_objects;use alloc_objects
use forcefield
use mdintegrators
implicit none
type(phasepoint_type),    intent(inout)::phasepoint
type(system_type),        intent(inout)   ::syst
type(MDinout_param_type), intent(inout)::MDinout_param
type(potential_type),     intent(in)   ::pot
type(dynamics_type),      intent(in)   ::dyn
type(timestep_type),      intent(in)   ::timestep



select case(dyn%DYNAMICS)
  case ("NVE")
    call velocity_Verlet(phasepoint%phasexv,syst,pot,MDinout_param%F,&
                         timestep)
  case ("ANDERSEN")
    call Andersen_velocity_Verlet(phasepoint%phasexv,syst,pot,& 
                                  MDinout_param%F,dyn%Andersen,timestep)
  case ("NOSEHOOVER")
    call Nose_Hoover_step(phasepoint,syst,pot,dyn%Nose_Hoover, &
                          MDinout_param%KIN,timestep)
  case ("LANGEVIN")
    call Lange_Verlet(phasepoint%phasexv,syst,pot,MDinout_param%F,&
                      dyn%Langevin)
  case("EXTERNAL")
     call externalMD(phasepoint,syst,dyn%ext_dyn,timestep%dt)
  case("RANDOMWALK")
     call make_randomwalk(phasepoint%phasexv%x(1,1))
  case default
    print *,"ERROR make_mdstep DYNAMICS=",dyn%DYNAMICS
    stop
end select
end subroutine make_mdstep 
!ES----------------------------------------------------------------

!BS----------------------------------------------------------------
! This calculates the MDinout_param parameter. This can be either the
! force (most cases) or the kinetic energy (Nose-Hoover). Starting a MD
! run will require to call this subroutine once at the start. The Verlet-type
! integration scheme already provides the new  MDinout_param parameter so
! that calling this routine again is unneccessary. The "EXTERNAL" option
! doesnot require a MDinout_param parameter. However, the external program
! needs to be initialized at the begin of a MD run.
!----------------------------------------------------------------------  
subroutine prepare_MDinout_param(MDinout_param,phasexv,syst,pot,dyn)
use types;use assign_objects;use alloc_objects
use forcefield
use stringlengths
use shell
implicit none
type(phasexv_type), intent(in)::phasexv
type(system_type),       intent(in)::syst
type(potential_type),    intent(in)::pot
type(dynamics_type),     intent(in)::dyn
type(MDinout_param_type), intent(inout)::MDinout_param
character(LEN=LSTR)::command
logical, save::first_visit
data first_visit /.true./


select case(dyn%DYNAMICS)
  case ("NOSEHOOVER") 
    MDinout_param%KIN=EKIN(phasexv%v,syst%masses_2,syst%Npart,syst%dim)
  case("NVE")
    MDinout_param%F=FORCE(phasexv%x,syst,pot)
  case("LANGEVIN")
    MDinout_param%F=FORCE(phasexv%x,syst,pot)
  case("ANDERSEN")
    MDinout_param%F=FORCE(phasexv%x,syst,pot)
  case("EXTERNAL")
    if (first_visit) then
      !!maybe not so elegant !!!!
      print *, "STARTING CPMD PATH SAMPLING OPTION"
      write(command,'(A11,i3,A50)') &
        "mpirun -np ",syst%NCPU," ./"//trim(syst%EXTERNAL_PROGRAM)// &
        " input.CPMD_TIS > output.CPMD_TIS &"
      call do_shell(command)
     first_visit=.false.
    endif
  case("RANDOMWALK")
    !just continue
  case DEFAULT
    print *,"ERROR prepare_MDinout_param:DYNAMICS=",dyn%DYNAMICS
    stop
end select

end subroutine prepare_MDinout_param
!ES----------------------------------------------------------------


end module mdstep
!EM------------------------
