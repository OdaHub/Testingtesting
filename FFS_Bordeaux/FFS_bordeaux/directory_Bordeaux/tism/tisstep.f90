!BM------------------------
module TISstep

contains

!BS--------------------------------------------------------------------
! This subroutine generates a new trajectory=trial, starting from an existing
! one. Two types of moves are considered: time-reversal and the shooting
!BS-------------------------------------------------------------------
subroutine make_TISstep(path,tisparam,syst,pot,dyn,timestep,trial,& 
                        wp1,wp2,wt1,wt2,allowmaxlength)
use types;use assign_objects;use alloc_objects
use random
implicit none
type(path_type),  intent(in)::path
type(path_type),  intent(inout)::trial
type(TIS_type),         intent(in)::TISparam
type(system_type),        intent(inout)   ::syst
type(potential_type),     intent(in)   ::pot
type(dynamics_type),      intent(in)   ::dyn
type(timestep_type),      intent(in)   ::timestep
type(path_type), intent(inout)::wp1,wp2  !workspace
logical, optional, intent(in)::allowmaxlength
type(timeslice_type), intent(inout)::wt1,wt2
double precision::ran


trial%index_acc=path%index_acc;trial%index_shoot=path%index_shoot
ran=random01()
if (ran<tisparam%timerevfreq) then
  call timereversal(path,trial,tisparam%startcondition)
else
  call shoot(path,tisparam,syst,pot,dyn,timestep,trial,wp1,wp2,wt1,wt2,allowmaxlength) 
endif

end subroutine make_TISstep 
!ES----------------------------------------------------------------

!BS----------------------------------------------------------------
! This subroutine performs the shooting move, which is the main MC
! move in the path-sampling approach. The standard shooting move is
! without the optional allowmaxlength. This implies that the maximum 
! length is obtained from the previous path length divided by a random number.
!-------------------------------------------------------------------
subroutine shoot(path,tisparam,syst,pot,dyn,timestep,trial,trajsegb,trajsegf,tshoot,wt2,allowmaxlength)
use types;use assign_objects;use alloc_objects
use random
use randomnize
use propagator
use cpmd_subr
use biasvacc
implicit none
type(path_type),  intent(in)::path
type(path_type),  intent(inout)::trial
type(TIS_type),         intent(in)::TISparam
type(system_type),        intent(inout)   ::syst
type(potential_type),     intent(in)   ::pot
type(dynamics_type),      intent(in)   ::dyn
type(timestep_type),      intent(in)   ::timestep
type(path_type), intent(inout)::trajsegb,trajsegf !workspace
type(timeslice_type), intent(inout)::tshoot,wt2       !workspace
logical, optional, intent(in)::allowmaxlength
double precision::ran
integer::ishoot,N,d,NT,NOPS,NX
logical::AM,reverseTF
double precision::INTERFACEL,INTERFACEM,INTERFACER,dK
integer::Lmax,LMAXb,LMAXf,Lpathb
character(LEN=LSTR)::Eold,Enew


print *,"SHOOT"

N=path%N;d=path%d;NT=path%NT;NOPS=path%NOPS;NX=path%NX
trial%MCmove="sh"
INTERFACEL=TISparam%INTERFACEL
INTERFACEM=TISparam%INTERFACEM
INTERFACER=TISparam%INTERFACER

!Pick a random point from the previous path
ran=random01()
ishoot=int(ran*(path%Lpath-2))+2  !random point, but not the end-points
tshoot=path%timeslices(ishoot)    !the shooting timeslice
!trial%OPshoot=tshoot%OPS(1)       !save RC value of the shooting point
!after the kick as the op can be v-dependent!
trial%iOPshoot_old=ishoot         !save the index on the old path of this point

!Change the velocities of timeslice tshoot
call kick_timeslice(tshoot,syst,pot,dyn,tisparam%sigdp_sqrtm,& 
                    tisparam%aimless,dK=dK,Eold=Eold,Enew=Enew)
trial%OPshoot=tshoot%OPS(1)       !save RC value of the shooting point
!accept or reject the new momenta
call Metropolis_momenta_change(dK,syst%beta,tisparam%aimless,&
        tshoot%OPS(1),tisparam%interfaceL,tisparam%interfaceR,AM)

if (tisparam%biasv) call biasv_acceptance(tshoot,path,ishoot,syst,timestep,pot,tisparam%startcondition,AM)

if (.not.AM) then
  !save trial failure "momenta change rejection (MCR)" and return 
  call trialfailure(trial,"MCR",tslice=tshoot)
  if (pot%POTENTIAL=="EXTERNAL") then
    print *,"RENAME BACK ESTRUCFILE BECAUSE OF REJECTION TRIAL"
    call do_shell_cpmd("mv "//trim(ENEW)//" "//trim(EOLD),.false.,syst%icrash,syst%NCPMD_UNSAVED)
  endif
  return
endif

if (present(allowmaxlength)) then
  LMAX=NX                    !Path may extend to maximum allocation limit
else                         !determine LMAX <= NX
  ran=random01()
  LMAX=int((path%Lpath-2)/ran)+2 !Lpath includes endpoints, number of shoot
                                 !possibilities excludes endpoints
  LMAX=MIN(LMAX,NX)              !Lpath=Lmax while finished is just OK
endif

LMAXB=LMAX-1                   !because forward trajectory is at least one 
                               !timeslice
reverseTF=.true.               !So, backward in time
call onewaytraj(tshoot,syst,pot,dyn,timestep,reverseTF,INTERFACEL,&
                  INTERFACER,LMAXB,trajsegb,wt2)
Lpathb=trajsegb%Lpath         !save backward trajectory
trial%iOPshoot_new=Lpathb     !save shooting index of the new trajectory

if (trajsegb%end=="*") then  
  !no ending at outer interfaces when LMAXB exceeded 
  !save data trial failure "backward trajectory too long (BTL)" and return
  call trialfailure(trial,"BTL",trsegb=trajsegb,im=INTERFACEM)
  if (pot%POTENTIAL=="EXTERNAL") then
    print *,"RENAME BACK ESTRUCFILE BECAUSE OF REJECTION TRIAL"
    call do_shell_cpmd("mv "//trim(ENEW)//" "//trim(EOLD),.false.,syst%icrash,syst%Ncpmd_unsaved)
  endif
  return
endif

!trajsegb%end=L or R and will be starting point of the trial path
!as the order of timeslices is reversed. If startcondition is " "
!the path may start at either sides (as in PPTIS). If startcondition
!is L/R this should be for trajsegb%end as well
if ((trajsegb%end/=tisparam%startcondition).AND.(tisparam%startcondition/=" ")) then !wrong backward ending 
  !save trial failure "backward trajectory at wrong interface (BWI)",return
  call trialfailure(trial,"BWI",trsegb=trajsegb,im=INTERFACEM) 
  if (pot%POTENTIAL=="EXTERNAL") then
    print *,"RENAME BACK ESTRUCFILE BECAUSE OF REJECTION TRIAL"
    call do_shell_cpmd("mv "//trim(ENEW)//" "//trim(EOLD),.false.,syst%icrash,syst%ncpmd_unsaved)
  endif
  return
endif

reverseTF=.false.    !so forward
LMAXF=LMAX-Lpathb+1  !+1 as startpoint (is also in trajsegb) will be removed 
                     ! afterwards
call onewaytraj(tshoot,syst,pot,dyn,timestep,reverseTF,INTERFACEL,&
                  INTERFACER,LMAXF,trajsegf,wt2)

!paste the backward and forward trajectories together
!update indexes etc for this path and check crossing-condition
call paste(trial,trajsegb,trajsegf,interfaceL,interfaceM,interfaceR,& 
           overlap=.true.)

if ((pot%POTENTIAL=="EXTERNAL").AND.(trial%accrej/="ACC")) then
      print *,"RENAME BACK ESTRUCFILE BECAUSE OF REJECTION TRIAL" 
      call do_shell_cpmd("mv "//trim(ENEW)//" "//trim(EOLD),.false.,syst%icrash,syst%ncpmd_unsaved)
endif

end subroutine shoot
!ES----------------------------------------------------------------

!BS-----------------------------------------------------------------
! This subroutine pastes two trajectory segament (backward and forward)
! together making the new trajectory trial. The following parameters of
! the trial path are  automatically set in this routine:
! Lpath,timeslices(:),opmax,iopmax,opmin,iopmin,start,end,cross
! index_acc,index_shoot (intent(inout))
! accrej (if accrej = NCR,FTL,FTX,ACC)
! Not set are:NX,N,d,NT,NOPS (set with allocation statement),
! MCmove, OPshoot,iOPshoot_old,iOPshoot_new (set in subroutine shoot)
! accrej=MCR,BTL,BTX, is set before (trialfailure backward segment)
!----------------------------------------------------------------------
subroutine paste(trial,trajb,trajf,intfL,intfM,intfR,overlap)
use types;use assign_objects;use alloc_objects
implicit none
type(path_type),  intent(inout)::trial
type(path_type),  intent(in)::trajb,trajf
double precision, intent(in)::intfL,intfM,intfR
logical, intent(in)::overlap
integer::Lpb,Lpf,Lp,NX,i

Lpb=trajb%lpath
Lpf=trajf%lpath

!do backward and forward trajectory have a common starting point?
!if yes, beware of duplication
NX=trial%NX
if (overlap) then
  Lp=Lpf+Lpb-1
else
  Lp=Lpf+Lpb
endif

trial%Lpath=Lp
do i=1,Lpb
  trial%timeslices(i)=trajb%timeslices(NX-Lpb+i)
enddo
if (overlap) then
  do i=1,Lpf-1
    trial%timeslices(Lpb+i)=trajf%timeslices(i+1)
  enddo
else
  do i=1,Lpf
    trial%timeslices(Lpb+i)=trajf%timeslices(i)
  enddo
endif
if (trajb%opmax>=trajf%opmax) then
  trial%opmax=trajb%opmax
  trial%iopmax=Lpb+1-trajb%iopmax
else
  trial%opmax=trajf%opmax
  if (overlap) then
    trial%iopmax=Lpb+trajf%iopmax-1
  else
    trial%iopmax=Lpb+trajf%iopmax
  endif
endif
if (trajb%opmin<=trajf%opmin) then
  trial%opmin=trajb%opmin
  trial%iopmin=Lpb+1-trajb%iopmin
else
  trial%opmin=trajf%opmin
  if (overlap) then
    trial%iopmin=Lpb+trajf%iopmin-1
  else
    trial%iopmin=Lpb+trajf%iopmin
  endif
endif
trial%start=trajb%end;trial%end=trajf%end
trial%cross="*"
if ((trial%opmin < intfM).AND.(trial%opmax >= intfM)) trial%cross="M"
if (trial%end=="*") then 
  trial%accrej="FTL"
  if (Lp==NX) trial%accrej="FTX"
else if (trial%cross=="*") then
  trial%accrej="NCR"
else
  trial%accrej="ACC"
  trial%index_acc=trial%index_acc+1
  trial%index_shoot=trial%index_shoot+1
endif

end subroutine paste
!ES----------------------------------------------------------------

!BS-----------------------------------------------------------------
! This subroutine reverses the path-trajectory to generate a trial 
! trajectory. The following trial parameters are set:
! Lpath,opmin,iopmin,opmax,iopmax,index_acc,index_shoot,start,end,cross,mcmove
! timeslices(:), accrej,opshoot,iopshoot_old,iopshoot_new
! NX,N,d,NT,NOPS are set before upon allocation
!--------------------------------------------------------------------------
subroutine timereversal(path,trial,startcondition)
use types;use assign_objects;use alloc_objects
implicit none
type(path_type),  intent(in)::path
type(path_type),  intent(inout)::trial
character, intent(in)::startcondition
integer::i,Lpath

Lpath=path%Lpath
trial%Lpath=Lpath
trial%opmin=path%opmin
trial%opmax=path%opmax
trial%iopmin=Lpath-path%iopmin+1
trial%iopmax=Lpath-path%iopmax+1
trial%index_acc=path%index_acc
trial%index_shoot=path%index_shoot
trial%start=path%end
trial%end=path%start
trial%cross=path%cross
trial%mcmove="tr"

print *,"TIME REVERSAL"

do i=1,Lpath

       trial%timeslices(i)=path%timeslices(Lpath+1-i)

       trial%timeslices(i)%phasepoint%phasexv%v= &
  -    trial%timeslices(i)%phasepoint%phasexv%v

       trial%timeslices(i)%phasepoint%extcoord%vxi= &
  -    trial%timeslices(i)%phasepoint%extcoord%vxi

       trial%timeslices(i)%phasepoint%electrons%rev_velec=&
  .NOT.trial%timeslices(i)%phasepoint%electrons%rev_velec

enddo

if (trial%start==startcondition) then
  trial%ACCREJ="ACC"
  trial%index_acc=trial%index_acc+1
else
  trial%ACCREJ="BWI"
endif

trial%OPshoot=0.d0
trial%iOPshoot_old=0
trial%iOPshoot_new=0


end subroutine timereversal
!ES-----------------------------------------------------------------

!BS-------------------------------------------------------------------
subroutine trialfailure(trial,rej,tslice,trsegb,im)
use types;use assign_objects;use alloc_objects
implicit none
type(path_type),  intent(inout)::trial
character*3, intent(in)        ::rej
type(timeslice_type),optional, intent(in)::tslice
type(path_type),     optional, intent(in)::trsegb
double precision, optional,    intent(in)::im
integer::NX,Lpath,i

NX=trial%NX

select case (rej)
  case("MCR")  !Momenta Change Rejection
    trial%Lpath=1
    trial%timeslices(1)=tslice
    trial%opmin=tslice%OPS(1)
    trial%opmax=tslice%OPS(1)
    trial%iopmin=1
    trial%iopmax=1
    trial%start="*";trial%end="*";trial%cross="*"
    trial%ACCREJ=rej
    trial%iOPshoot_new=1
  case("BWI") !Backward trajectory ends at Wrong Interface
    Lpath=trsegb%Lpath
    trial%Lpath=Lpath
    !trial%timeslices(1:Lpath)=trsegb%timeslices(NX-LPATH:NX)
    do i=1,Lpath
      trial%timeslices(i)=trsegb%timeslices(NX-LPATH+i)
    enddo
    trial%opmin=trsegb%opmin
    trial%opmax=trsegb%opmax
    trial%iopmin=Lpath+1-trsegb%iopmin
    trial%iopmax=Lpath+1-trsegb%iopmax
    trial%start=trsegb%end;trial%end="*";trial%cross="*"
    if (trial%opmin < IM) trial%cross="M"
    trial%ACCREJ=rej
  case("BTL") !Backward trajectory Too Long
    Lpath=trsegb%Lpath
    trial%Lpath=Lpath
    !trial%timeslices(1:Lpath)=trsegb%timeslices(NX-LPATH:NX)
    do i=1,Lpath
      trial%timeslices(i)=trsegb%timeslices(NX-LPATH+i)
    enddo
    trial%opmin=trsegb%opmin
    trial%opmax=trsegb%opmax
    trial%iopmin=Lpath+1-trsegb%iopmin
    trial%iopmax=Lpath+1-trsegb%iopmax
    trial%start="*";trial%end="*";trial%cross="*"
    if ((trial%opmin < IM).and.(trial%opmax>IM)) trial%cross="M"
    trial%ACCREJ=rej
    if (Lpath==NX-1) trial%ACCREJ="BTX"
  case default
    print *,"ERROR trialfailur, rej=",rej
end select 
!other failures, such as the trajectory that fails to cross the middle 
!interface, are handled in the "paste"-subroutine

end subroutine trialfailure
!ES-------------------------------------------------------------------

end module TISstep
!EM------------------------
