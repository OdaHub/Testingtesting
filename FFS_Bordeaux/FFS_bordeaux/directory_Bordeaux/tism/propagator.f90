!BM------------------------
module propagator 

contains
!-----------------------------------------------------------------------
! This subroutine propagates a timeslice forward or backward in time
! by a single timestep
!-----------------------------------------------------------------------
subroutine propagate(timeslice,syst,pot,dyn,timestep,backwardTF)
use types;use assign_objects;use alloc_objects
use mdstep 
use orderparameter
implicit none
type(timeslice_type),  intent(inout)::timeslice
type(system_type),        intent(inout)   ::syst
type(potential_type),     intent(in)   ::pot
type(dynamics_type),      intent(in)   ::dyn
type(timestep_type),      intent(in)   ::timestep
logical, intent(in)::backwardTF
integer::N,d,NT
type(timestep_type)::timestep2
N=timeslice%N
d=timeslice%d
NT=timeslice%NT
  
  timestep2=timestep
  if (backwardTF) then !going backward in time
     if (pot%potential=="EXTERNAL") then 
       !Give negative timestep for the "PSAMPLE" file
       timestep2%dt=-timestep%dt 
     else
       ! reverse momenta
       timeslice%phasepoint%phasexv%v=-timeslice%phasepoint%phasexv%v
       timeslice%phasepoint%extcoord%vxi=-timeslice%phasepoint%extcoord%vxi
     endif
  endif
  
  call make_mdstep(timeslice%phasepoint,syst,pot,dyn, &
                   timeslice%MDinout_param,timestep2)
 
  if (backwardTF) then 
    if (pot%potential/="EXTERNAL") then
       !backwardTF=.true. => reverse momenta again
       timeslice%phasepoint%phasexv%v=-timeslice%phasepoint%phasexv%v
       timeslice%phasepoint%extcoord%vxi=-timeslice%phasepoint%extcoord%vxi
     endif
  endif

  !Update the information content of the timeslice
  timeslice%OPS(1)=orderp(timeslice%phasepoint%phasexv%x,syst,timeslice%phasepoint%phasexv%v,pot)
end subroutine propagate 
!ES----------------------------------------------------------------

!BS-------------------------------------------------------------------------
! This subroutine generates a trajectory, a sequence of consequetive timeslices, 
! starting from  an input timeslice going either forward or backward in time.
! Trajectory stops whenever INTERFACEL or INTERFACER are crossed, or when the path length
! becomes LMAX
!--------------------------------------------------------------------------------------
subroutine onewaytraj(timeslice,syst,pot,dyn,timestep,backwardTF,INTERFACEL,&
                  INTERFACER,LMAX,traj,tslice)
use types;use assign_objects;use alloc_objects
use cpmd_subr
implicit none
type(timeslice_type),  intent(in)::timeslice
type(system_type),        intent(inout)   ::syst
type(potential_type),     intent(in)   ::pot
type(dynamics_type),      intent(in)   ::dyn
type(timestep_type),      intent(in)   ::timestep
logical, intent(in)::backwardTF
double precision, intent(in)::INTERFACEL,INTERFACER
integer,intent(in)::LMAX
type(path_type), intent(out)::traj
type(timeslice_type)::tslice   !workspace
integer::N,d,NT,NOPS,sign,startloc,LPATH,loc,NX
double precision::opmin,opmax
integer::iopmin,iopmax
character::end

!extract some number from the timeslice
N=timeslice%n
d=timeslice%d
NT=timeslice%NT
NOPS=timeslice%NOPS
NX=traj%NX

!put default value for the ending
end="*"

if (backwardTF) then  !start trajectory at place NX and then go downwards
  sign=-1
  startloc=NX+1
else                  !start trajectory at place 1 and then go upwards
  sign=1
  startloc=0
endif

LPATH=1
loc=startloc+sign*Lpath
traj%timeslices(loc)=timeslice

!initial opmin,opmax etc
iopmin=1
iopmax=1
opmin=timeslice%OPS(1)
opmax=timeslice%OPS(1)

tslice=timeslice
do 
  if (opmin < INTERFACEL) then
    end="L"
    exit
  endif
  if (opmax > INTERFACER) then
    end="R"
    exit
  endif
  if (LPATH==LMAX) exit

  !propagate backward/forward one timestep
  LPATH=LPATH+1
  call propagate(tslice,syst,pot,dyn,timestep,backwardTF) 

  !save this timeslice inside trajectory
  loc=startloc+sign*Lpath
  traj%timeslices(loc)=tslice

  !opmin and opmax might need to be updated
  if (.NOT.syst%NOPRINT) print *,"orderp:",tslice%OPS(1)
  if (tslice%OPS(1) < opmin) then
    opmin=tslice%OPS(1)
    iopmin=LPATH
  endif
  if (tslice%OPS(1) > opmax) then
    opmax=tslice%OPS(1)
    iopmax=LPATH
  endif
enddo

!finally, save some important information for this trajectory
traj%opmin=opmin
traj%opmax=opmax
traj%iopmin=iopmin
traj%iopmax=iopmax
traj%Lpath=Lpath
traj%end=end


end subroutine onewaytraj 
!ES----------------------------------------------------------------

end module propagator 
!EM------------------------
