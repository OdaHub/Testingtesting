!BM-----------------------------------------------------------
module modswapfm
implicit none

CONTAINS
!BS-------------------------------------------------------------------
subroutine swapfm(iswap,PPS_ensemble,tisparam,syst,pot,dyn,output,timestep,&
                     icyc,wp1,wp2,wp3,wp4,wt1,wt2,secpot)
use types;use assign_objects;use alloc_objects
use propagator
use TISstep
use outTIS
use convert
use forcefield
use orderparameter
use random
implicit none
integer, intent(in)::iswap
type(path_ensemble), intent(inout)::PPS_ensemble
type(TIS_type),       intent(inout)::TISparam
type(system_type),        intent(inout)   ::syst
type(potential_type),     intent(in)   ::pot,secpot
type(dynamics_type),      intent(in)   ::dyn
type(output_type), intent(inout):: output
type(timestep_type),      intent(in)   ::timestep
integer,                  intent(in)   ::icyc
type(path_type), intent(inout)::wp1,wp2,wp3,wp4   !workspace
type(timeslice_type), intent(inout)::wt1,wt2      !workspace
double precision::ran
integer::Lmax,j,startpos,Npart,dim
logical::reverseTF,overlap
type(potential_type)::pt
!type(timeslice_type)::timesliceA,timesliceB
!type(phasepoint_type)::phaseA,phaseB
!requires alloation, so have to use workspace wt1,wt2,..
integer::indexa,indexs
double precision::rannum,OPA,OPB
logical::go_on
double precision::Vpswitch,Vmswitch,Vp,Vm,Pacc,DV,Econst
double precision::minusinf,lam0, L0mold
integer::NX

!!print *,"swapfm"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  First check Metropolis rule based on energies     !!
!!!!  No integration yet                                !!
!!!!  or rescale velocities to appropriate fixed energy !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Npart=syst%Npart
dim=syst%dim
NX=tisparam%NX

L0mold=PPS_ensemble%Paths(1)%Lpath

wt1=PPS_ensemble%Paths(2)%timeslices(1)
!1st timeslice of [0+] path
!this will be the one but last point of the new [0-] path
startpos=PPS_ensemble%Paths(1)%Lpath-1
wt2=PPS_ensemble%Paths(1)%timeslices(startpos)
!the old one-but-last point of the old [0-] path
!! this will not yet work for CPMD
Vpswitch=Epot(wt2%phasepoint%phasexv%x,syst,pot)
Vmswitch=Epot(wt1%phasepoint%phasexv%x,syst,secpot)
Vp=Epot(wt1%phasepoint%phasexv%x,syst,pot)
Vm=Epot(wt2%phasepoint%phasexv%x,syst,secpot)
go_on=.true.

if ((dyn%DYNAMICS=="NVE").AND.(dyn%NVE%RESCALE_ENERGY)) then
  !rescale velocities if possible
  !accept accrding to the phasespace volume expansion
  Econst=dyn%NVE%NVE_energy
  if ((Vpswitch>Econst).OR.(Vmswitch>Econst)) then
    go_on=.false. 
  else
    Pacc=(  (Econst-Vpswitch)*(Econst-Vmswitch)/  &
          ( (Econst-Vp      )*(Econst-Vm      ) ) )**(dim*Npart/2.d0)
  endif
else
  DV=Vpswitch+Vmswitch-Vp-Vm
  Pacc=exp(-syst%beta*DV)
  !!write(777,*) "a", Vpswitch,Vmswitch,Vp,Vm
  !!write(777,*) "b",syst%beta, DV
  !!write(777,*) "c", Pacc,icyc
endif

if ((go_on).AND.((1.d0-Pacc)>1.d-15)) then
  !!call random_number(rannum)
  rannum=random01()
  if (rannum>pacc) go_on=.false.
endif
!!go_on=.true. !remove thsi is a test
!!if (.not.go_on) write(777,*) icyc,"rej@once"

!------------------------------------------------------------
!set default values for both trajectories that is used
!for output of rejections 
!values are overwritten when some integration is taken place
!------------------------------------------------------------
wp1%Lpath=1;wp2%Lpath=1
wp1%mcmove="s+";wp2%mcmove="s-"
wp1%ACCREJ="MCR";wp2%ACCREJ="MCR"
wp1%start="*";wp2%start="*"
wp1%cross="*";wp2%cross="*"
wp1%end="*";wp2%end="*"
wp1%index_acc=PPS_ensemble%PATHS(1)%index_acc
wp2%index_acc=PPS_ensemble%PATHS(2)%index_acc
wp1%index_shoot=PPS_ensemble%PATHS(1)%index_shoot
wp2%index_shoot=PPS_ensemble%PATHS(2)%index_shoot
wp1%iopshoot_old=0;wp1%iopshoot_new=0;wp2%iopshoot_old=0;wp2%iopshoot_new=0
wp1%timeslices(1)=PPS_ensemble%PATHS(2)%timeslices(1)
startpos=PPS_ensemble%Paths(1)%Lpath-1
wp2%timeslices(1)=PPS_ensemble%PATHS(1)%timeslices(startpos)
wp1%opmin=orderp(PPS_ensemble%PATHS(2)%timeslices(1)%phasepoint% &
  phasexv%x,syst,PPS_ensemble%PATHS(2)%timeslices(1)%phasepoint% &
  phasexv%v,secpot)
wp2%opmin=orderp(PPS_ensemble%PATHS(1)%timeslices(startpos)%phasepoint% &
  phasexv%x,syst,PPS_ensemble%PATHS(1)%timeslices(startpos)%phasepoint% &
  phasexv%v,pot)
wp1%opmax=wp1%opmin
wp2%opmax=wp2%opmin
wp1%opshoot=wp1%opmin;wp2%opshoot=wp2%opmin
wp1%iopshoot_old=1;wp2%iopshoot_old=startpos
wp1%iopshoot_new=1;wp2%iopshoot_new=1
wp1%iopmax=0;wp1%iopmin=0;wp2%iopmax=0;wp2%iopmin=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! integrate new [0-] point one step to create 2p-traj wp1 !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (go_on) then
  !wt1: already set
  !wt1: 1st timeslice of old [0+] path, L-1 timeslice of new [0-] path
  !check if new [0-] point jumps over the interface in one forward step
  wt1%OPS(1)=orderp(wt1%phasepoint%phasexv%x,syst,wt1%phasepoint%phasexv%v,&
                    secpot)
  !orderparameter might depend on potential.
  !point from [0+] to generate path in [0-]
  call set_TIS(1,pps_ensemble,output=output,tisparam=tisparam)
  reverseTF=.false.  !forward 
  !save tisparam%INTERFACEL (-Infinity) tisparam%INTERFACER (lambda_0)
  minusinf=tisparam%INTERFACEL
  lam0= tisparam%INTERFACER
  !!print *,"now!"
  call onewaytraj(wt1,syst,secpot,dyn,timestep,reverseTF,&
                  minusinf,lam0,2,wp1,wt2)
  !!print *,"df"
  !!print *,wp1%timeslices(1)%OPS(1)
  !!print *,wp1%timeslices(2)%OPS(1)
  !!stop
  if ((wp1%Lpath==1).OR.(wp1%opmax<lam0)) go_on=.false. 
endif
!!if (.not.go_on) write(777,*) icyc,"rec@2"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! integrate new [0+] point one step to create 2p-traj wp2 !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (go_on) then
  !wt2 has been overwritten for working space
  startpos=PPS_ensemble%Paths(1)%Lpath-1
  wt1=PPS_ensemble%Paths(1)%timeslices(startpos)
  !the old one-but-last point of the old [0-] path
  !check if new [0+] point jumps over the interface in one backward step
  wt1%OPS(1)=orderp(wt1%phasepoint%phasexv%x,syst,wt1%phasepoint%phasexv%v,pot)
  !orderparameter might depend on potential.
  call set_TIS(2,pps_ensemble,output=output,tisparam=tisparam)
  reverseTF=.false.  !forward as well
  call onewaytraj(wt1,syst,pot,dyn,timestep,reverseTF,&
                  minusinf,lam0,2,wp2,wt2)
  if ((wp2%Lpath==1).OR.(wp2%opmax<lam0)) go_on=.false.
endif
!!if (.not.go_on) write(777,*) icyc,"rec@3"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Finalize [0-] traj wp1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (go_on) then
  !complete wp1 backward in time
  reverseTF=.true. !backward
  LMAX=NX-1
  wt1=wp1%timeslices(1)
  call set_TIS(1,pps_ensemble,output=output,tisparam=tisparam)
  call onewaytraj(wt1,syst,secpot,dyn,timestep,reverseTF, &
                  tisparam%INTERFACEL,&
                  tisparam%INTERFACER,LMAX,wp3,wt2)
  wp4=wp1
  !now check if not LMAX is exceeded
  !otherwise we just overwrite pps_ensemble%Paths(1)
  
  !!if ((wp3%LPATH==LMAX).AND.& 
  !!  (wp3%timeslices(LMAX)%OPS(1)<tisparam%INTERFACER)) then
  !!  print *,"RESULTS MIGHT BE IN ERROR INCREASE NX!!!"
  !!    go_on=.false. 
  !!    
  !!endif

  !!if (go_on) then
  call paste(wp1,wp3,wp4,tisparam%interfaceL,& 
             tisparam%interfaceM,tisparam%interfaceR,overlap=.true.)
  !!call paste(pps_ensemble%Paths(1),wp3,wp1,tisparam%interfaceL,&
  !!           tisparam%interfaceM,tisparam%interfaceR,overlap=.true.)
  !!wp1=pps_ensemble%Paths(1) !for output later on
  !!wp1%MCmove="s+" 
  wp1%OPshoot=wp1%timeslices(wp1%Lpath-1)%OPS(1)
  wp1%iOPshoot_old=1
  wp1%iOPshoot_new=wp1%Lpath-1
  wp1%ACCREJ="ACC"
  if (wp1%LPATH==NX) then
    wp1%ACCREJ="BTX"
    go_on=.false.
  endif
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Finalize [0+] traj wp2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (go_on) then
  !complete wp2 forward in time
  reverseTF=.false. !forward
  LMAX=NX-1
  call set_TIS(2,pps_ensemble,output=output,tisparam=tisparam)
  wt1=wp2%timeslices(2)
  call onewaytraj(wt1,syst,pot,dyn,timestep,reverseTF,& 
                  tisparam%INTERFACEL,&
                  tisparam%INTERFACER,LMAX,wp3,wt2)
  
  wp4%timeslices(NX)=wp2%timeslices(1) !put in backward order 
  wp4%Lpath=1                          !to allow pasting
  wp4%opmin=wp4%timeslices(NX)%OPs(1)
  wp4%opmax=wp4%opmin
  wp4%iopmin=1
  wp4%iopmax=1
  wp4%end="L"

  !!endif 

  !!if ((wp3%LPATH==LMAX).AND.&
  !!(wp3%timeslices(LMAX)%OPS(1)<tisparam%INTERFACER)) then
  !!  print *,"RESULTS MIGHT BE IN ERROR INCREASE NX!!!"
  !!  go_on=.false. 
  !!endif

  !!if (go_on) then
  call paste(wp2,wp4,wp3,tisparam%interfaceL,&
                         tisparam%interfaceM,&
                         tisparam%interfaceR, overlap=.false.)
  !!call paste(pps_ensemble%Paths(2),wp2,wp3,tisparam%interfaceL,&
  !!                       tisparam%interfaceM,&
  !!                       tisparam%interfaceR, overlap=.false.)
  !!wp2=pps_ensemble%Paths(2) !for output later on
  wp2%MCmove="s-"
  wp2%OPshoot=wp2%timeslices(1)%OPS(1)
  wp2%iOPshoot_old=L0mold-1
  wp2%iOPshoot_new=1
  wp2%ACCREJ="ACC"
  if (wp2%LPATH==NX) wp2%ACCREJ="FTX"
endif

!output for wp1 and wp2
call set_TIS(1,pps_ensemble,output=output)
call outputTIS(icyc,wp1,syst,secpot,dyn,output)
call set_TIS(2,pps_ensemble,output=output)
call outputTIS(icyc,wp2,syst,pot,dyn,output)
if ((wp1%ACCREJ=="ACC").AND.(wp2%ACCREJ=="ACC")) then
  pps_ensemble%Paths(1)=wp1
  pps_ensemble%Paths(2)=wp2
endif

end subroutine swapfm
!ES----------------------------------------------------------------------

end module modswapfm
