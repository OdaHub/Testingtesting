!BM------------------------
module PPSstep

contains

!BS-------------------------------------------------------------------
! This subroutine makes one cycle in the TIS-PPS cyclus. Either some
! (or all) interface simulations are updated by the standard TIS move
! (shooting or timereversal) or one or more swaps are performed.
!---------------------------------------------------------------------
subroutine make_PPSstep(PPS_ensemble,tisparam,syst,pot,dyn,output, &
                        timestep,icyc,path,trial,wp1,wp2,wp3,wt1,wt2,secpot)
use types;use assign_objects;use alloc_objects
use modswap
use random
implicit none
type(path_ensemble),      intent(inout)::PPS_ensemble
type(tis_type),           intent(inout)::tisparam
type(system_type),        intent(inout)   ::syst
type(potential_type),     intent(in)   ::pot,secpot
type(dynamics_type),      intent(in)   ::dyn
type(output_type),        intent(inout)::output
type(timestep_type),      intent(in)   ::timestep
integer,                  intent(in)   ::icyc
type(path_type), intent(inout)::path,trial,wp1,wp2,wp3 !workspace
type(timeslice_type), intent(inout)::wt1,wt2
double precision::ran




!!call random_number(ran)
ran=random01()
if (ran < PPS_ensemble%swapfreq) then
  !Do one or more swaps
  call make_swap(PPS_ensemble,tisparam,syst,pot,dyn,output,timestep,icyc,wp1,& 
                 wp2,path,wp3,wt1,wt2,secpot) 
else
  call make_TISmoves(PPS_ensemble,tisparam,syst,pot,dyn,output,timestep,icyc,&
                     path,trial,wp1,wp2,wt1,wt2,secpot)
  !Update one or more interface ensembles by the standard TIS moves
endif

end subroutine make_PPSstep 
!ES----------------------------------------------------------------

!BS-------------------------------------------------------------------
subroutine make_TISmoves(PPS_ensemble,tisparam,syst,pot,dyn,output,timestep,&
                         icyc,path,trial,wp1,wp2,wt1,wt2,secpot)
use types;use assign_objects;use alloc_objects
use TISstep
use outTIS
use cpmd_subr
use convert
use random
implicit none
type(path_ensemble), intent(inout)::PPS_ensemble
type(TIS_type),         intent(inout)::TISparam
type(system_type),        intent(inout)   ::syst
type(potential_type),     intent(in)   ::pot,secpot
type(dynamics_type),      intent(in)   ::dyn
type(output_type), intent(inout):: output
type(timestep_type),      intent(in)   ::timestep
integer, intent(in)::icyc
type(path_type), intent(inout)::path,trial,wp1,wp2  !workspace
type(timeslice_type), intent(inout)::wt1,wt2
double precision::ran,sumweights
integer::i,j,N,d,NT,NOPS,NX,NUMINT
character(LEN=LSTR)::Efile
logical::revab
type(potential_type)::pt


N=pps_ensemble%N;d=pps_ensemble%d;NT=pps_ensemble%NT
NOPS=pps_ensemble%NOPS;NX=pps_ensemble%NX;numint=pps_ensemble%numint


if (PPS_ensemble%RELATIVESHOOTS) then

  !!call random_number(ran)
  ran=random01()
  sumweights=0.d0
  do i=1,numint
    sumweights=sumweights+PPS_ensemble%RELATIVE_SHOOTFREQ(i)
    if (ran<sumweights) exit
  enddo 

  ! path=pps_ensemble%paths(i), tisparam and output must 
  ! be adapted to i
  print *,"make TIS move ens",i
  call set_TIS(i,pps_ensemble,path,output,pt,pot,secpot,tisparam)
  call make_TISstep(path,tisparam,syst,pt,dyn,timestep,trial,wp1,wp2,wt1,wt2)
  call outputTIS(icyc,trial,syst,pt,dyn,output)
  if (trial%ACCREJ=="ACC") pps_ensemble%paths(i)=trial

  !remove old RESTART FILES
  !!EXTERNAL not yet compatible with forcefield
  if (pot%POTENTIAL=="EXTERNAL") then
    EFILE=pps_ensemble%paths(i)%timeslices(1)%phasepoint%electrons%ESTRUCFILE
    call delete_efiles(EFILE,revab,syst%icrash,syst%NCPMD_UNSAVED)
  endif

  
  if (PPS_ensemble%nullmoves) then
    do j=1,numint
      if (j==i) cycle
      !set path and output for j
      call set_TIS(j,pps_ensemble,path,output,pt,pot,secpot)
      path%MCmove="00"
      path%index_acc=path%index_acc+1
      PPS_ensemble%PATHS(j)%index_acc=PPS_ensemble%PATHS(j)%index_acc+1
       
      call outputTIS(icyc,path,syst,pt,dyn,output)    
    enddo
  endif

else

  do i=1,numint
    print *,"make TIS move ens",i
    !!repeats a bit previous lines, could be nicer
    call set_TIS(i,pps_ensemble,path,output,pt,pot,secpot,tisparam)
    call make_TISstep(path,tisparam,syst,pt,dyn,timestep,trial,wp1,wp2,wt1,wt2)
    call outputTIS(icyc,trial,syst,pt,dyn,output)
    if (trial%ACCREJ=="ACC")  pps_ensemble%paths(i)=trial

    !remove old RESTART FILES
    if (pt%POTENTIAL=="EXTERNAL") then
      EFILE=pps_ensemble%paths(i)%timeslices(1)%phasepoint%electrons%ESTRUCFILE
      call delete_efiles(EFILE,revab,syst%icrash,syst%NCPMD_UNSAVED)
    endif

  enddo

endif


end subroutine make_TISmoves
!ES----------------------------------------------------------------

end module PPSstep
!EM------------------------
