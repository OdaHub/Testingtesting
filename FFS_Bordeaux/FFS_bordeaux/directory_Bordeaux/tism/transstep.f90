!BM------------------------
module transstep

contains

!BS-------------------------------------------------------------------
subroutine make_transstep(icyc,startpoint,startpoint2,syst,pot,dyn,output,&
                          timestep,interfaceL,interfaceM,interfaceR,&
                          transalgorithm, two_point_method,&
                          wp1,wp2,wp3,wt1,wt2)
use types;use assign_objects;use alloc_objects
use orderparameter
use TISstep
use outTIS
use propagator
use mdstep
implicit none
integer, intent(in)::icyc
type(phasepoint_type),  intent(in)::startpoint,startpoint2
type(system_type),        intent(inout)   ::syst
type(potential_type),     intent(in)   ::pot
type(dynamics_type),      intent(in)   ::dyn
type(output_type), intent(in)::output
type(timestep_type),      intent(in)   ::timestep
double precision,      intent(in) ::interfaceL,interfaceM,interfaceR
character(LEN=MSTR), intent(in):: TRANSALGORITHM
logical,      intent(in)::two_point_method
type(path_type), intent(inout)::wp1,wp2,wp3  !workspace
type(timeslice_type), intent(inout)::wt1,wt2
double precision::dotx,chi,chiTST
integer::LMAX
logical::reverseTF
double precision::epsilon


select case(TRANSALGORITHM)
  case("EPF")
print *, 'EPF'
    call epfstep(icyc,startpoint,startpoint2,syst,pot,dyn,output,timestep,&
                          interfaceL,interfaceM,interfaceR,transalgorithm, &
                          two_point_method,wp1,wp2,wp3,wt1,wt2)
  case("BC")
print *, 'BC'
       call BCstep(icyc,startpoint,startpoint2,syst,pot,dyn,output,timestep,&
                          interfaceL,interfaceM,interfaceR,transalgorithm, &
                          two_point_method,wp1,wp2,wp3,wt1,wt2)
  case default
   print *,"ERROR transstep TRANSALGORITHM=",TRANSALGORITHM
   stop
end select

end subroutine make_transstep
!ES----------------------------------------------------------------

!BS-------------------------------------------------------------------
subroutine epfstep(icyc,startpoint,startpoint2,syst,pot,dyn,output,timestep,&
                          interfaceL,interfaceM,interfaceR,transalgorithm, &
                          two_point_method,wp1,wp2,wp3,wt1,wt2)
use types;use assign_objects;use alloc_objects
use orderparameter
use TISstep
use outTIS
use propagator
use mdstep
use convert
implicit none
integer, intent(in)::icyc
type(phasepoint_type),  intent(in)::startpoint,startpoint2
type(system_type),        intent(inout)   ::syst
type(potential_type),     intent(in)   ::pot
type(dynamics_type),      intent(in)   ::dyn
type(output_type), intent(in)::output
type(timestep_type),      intent(in)   ::timestep
double precision,      intent(in) ::interfaceL,interfaceM,interfaceR
character(LEN=MSTR), intent(in):: TRANSALGORITHM
logical, intent(in)::two_point_method
type(path_type), intent(inout)::wp1,wp2,wp3  !workspace
type(timeslice_type), intent(inout)::wt1,wt2
double precision::dotx,chi,chiTST
integer::LMAX
logical::reverseTF
double precision::epsilon,op1,op2




epsilon=1.d-9
LMAX=wp1%NX

if (two_point_method) then
  dotx=1.d0
else
  if ((dyn%DYNAMICS=="LANGEVIN").AND.& 
      (dyn%langevin%high_friction_limit)) then
    dotx=1.d0  !!not so neat
  else
    dotx=v_orderp(startpoint%phasexv%x,startpoint%phasexv%v,syst)
  endif
endif

if (two_point_method) then



  op1=orderp(startpoint%phasexv%x,syst,startpoint%phasexv%v,pot) 
  op2=orderp(startpoint2%phasexv%x,syst,startpoint%phasexv%v,pot) 
  chiTST=1.d0
  call convert_phasep2tslice(startpoint,syst,pot,dyn,wt1)

  chiTST=1.d0
  reverseTF=.true.
  call onewaytraj(wt1,syst,pot,dyn,timestep,reverseTF,INTERFACEL,&
                      INTERFACEM,LMAX,wp1,wt2)
  if (wp1%end=="R") then

    chi=0.d0
    call trialfailure(wp3,"BWI",trsegb=wp1,im=INTERFACEM)

  else if (wp1%end=="L") then

    reverseTF=.false.
    LMAX=wp1%NX-wp1%Lpath

    call convert_phasep2tslice(startpoint2,syst,pot,dyn,wt1)
    call onewaytraj(wt1,syst,pot,dyn,timestep,reverseTF,INTERFACEL,&
                      INTERFACER-epsilon,LMAX,wp2,wt2)
    call paste(wp3,wp1,wp2,interfaceL,interfaceM,interfaceR,overlap=.false.)

    if (wp3%end=="R") then
      chi=1.d0
    else
      chi=0.d0
    endif

  else
    print *,"NX is too low"
    chi=0.d0
  endif


else

  call convert_phasep2tslice(startpoint,syst,pot,dyn,wt1) 
  if (dotx<0) then
    chi=0.d0
    ChiTST=0.d0
    call trialfailure(wp3,"MCR",tslice=wt1)
  else
    chiTST=1.d0
    reverseTF=.true.
    call onewaytraj(wt1,syst,pot,dyn,timestep,reverseTF,INTERFACEL,&
                      INTERFACEM,LMAX,wp1,wt2)  
    if (wp1%end=="R") then
      chi=0.d0
      call trialfailure(wp3,"BWI",trsegb=wp1,im=INTERFACEM)
    else if (wp1%end=="L") then
      reverseTF=.false.
      LMAX=wp1%NX-wp1%Lpath+1
      call onewaytraj(wt1,syst,pot,dyn,timestep,reverseTF,INTERFACEL,&
                      INTERFACER-epsilon,LMAX,wp2,wt2)
      call paste(wp3,wp1,wp2,interfaceL,interfaceM,interfaceR,overlap=.true.)
      if (wp3%end=="R") then
        chi=1.d0
      else
        chi=0.d0
      endif 
    else 
      print *,"NX is too low"
      chi=0.d0
    endif
  endif
  if ((dyn%DYNAMICS=="LANGEVIN").AND.&
      (dyn%langevin%high_friction_limit)) then 
  !   dotx=(wp1%timeslices(wp1%NX)%OPS(1)-wp1%timeslices(wp1%NX-1)%OPS(1))/&
  !         timestep%dt
  dotx=wp1%timeslices(wp1%NX-1)%phasepoint%phasexv%v(1,1)/timestep%dt
  !  dotx=(wp2%timeslices(2)%OPS(1)-wp2%timeslices(1)%OPS(1))/&
  !         timestep%dt
  
     if (dotx<0) chiTST=0
  endif
endif
  
wp3%MCmove="sh"
wp3%OPshoot=interfaceM
wp3%iOPshoot_old=1
wp3%iOPshoot_new=1
wp3%accrej="ACC"
wp3%index_acc=icyc
wp3%index_shoot=icyc
call outputTIS(icyc,wp3,syst,pot,dyn,output)
write(output%iutrans,*) icyc, chi*dotx,chiTST*dotx

end subroutine epfstep 
!ES----------------------------------------------------------------


!BS-------------------------------------------------------------------
subroutine BCstep(icyc,startpoint,startpoint2,syst,pot,dyn,output,timestep,&
                          interfaceL,interfaceM,interfaceR,transalgorithm, &
                          two_point_method,wp1,wp2,wp3,wt1,wt2)
use types;use assign_objects;use alloc_objects
use orderparameter
use TISstep
use outTIS
use propagator
use mdstep
use convert
implicit none
integer, intent(in)::icyc
type(phasepoint_type),  intent(in)::startpoint,startpoint2
type(system_type),        intent(inout)   ::syst
type(potential_type),     intent(in)   ::pot
type(dynamics_type),      intent(in)   ::dyn
type(output_type), intent(in)::output
type(timestep_type),      intent(in)   ::timestep
double precision,      intent(in) ::interfaceL,interfaceM,interfaceR
character(LEN=MSTR), intent(in):: TRANSALGORITHM
logical, intent(in)::two_point_method
type(path_type), intent(inout)::wp1,wp2,wp3  !workspace
type(timeslice_type), intent(inout)::wt1,wt2
double precision::dotx,chi,chiTST
integer::LMAX
logical::reverseTF
double precision::epsilon,op1,op2

epsilon=1.d-9
LMAX=wp1%NX

if (two_point_method) then
  dotx=1.d0
else
  if ((dyn%DYNAMICS=="LANGEVIN").AND.&
      (dyn%langevin%high_friction_limit)) then
      dotx=1.d0  !!not so neat
  else
    dotx=v_orderp(startpoint%phasexv%x,startpoint%phasexv%v,syst)
  endif
endif

if (two_point_method) then
  op1=orderp(startpoint%phasexv%x,syst,startpoint%phasexv%v,pot)
  op2=orderp(startpoint2%phasexv%x,syst,startpoint%phasexv%v,pot)
  chiTST=1.d0
  if (op1> op2) then
     call convert_phasep2tslice(startpoint,syst,pot,dyn,wt1)
  else
     call convert_phasep2tslice(startpoint2,syst,pot,dyn,wt1)
  endif
  reverseTF=.false.
  call onewaytraj(wt1,syst,pot,dyn,timestep,reverseTF,INTERFACEL,&
                    INTERFACER,LMAX,wp1,wt2)
  if (wp1%end=="R") then
    chi=1.d0
  else if (wp1%end=="L") then
    chi=0.d0
  else
    print *,"NX is too low"
    chi=0.d0
  endif

else
  call convert_phasep2tslice(startpoint,syst,pot,dyn,wt1)
  chiTST=1.d0;if (dotx<0) chiTST=0.d0
  reverseTF=.false.
  call onewaytraj(wt1,syst,pot,dyn,timestep,reverseTF,INTERFACEL,&
                    INTERFACER,LMAX,wp1,wt2)
  if (wp1%end=="R") then
    chi=1.d0
  else if (wp1%end=="L") then
    chi=0.d0
  else
    print *,"NX is too low"
    chi=0.d0
  endif
  if ((dyn%DYNAMICS=="LANGEVIN").AND.&
      (dyn%langevin%high_friction_limit)) then
    dotx=wp1%timeslices(2)%phasepoint%phasexv%v(1,1)/timestep%dt
    if (dotx<0) chiTST=0
  endif
endif

wp1%MCmove="sh"
wp1%OPshoot=interfaceM
wp1%iOPshoot_old=1
wp1%iOPshoot_new=1
wp1%accrej="ACC"
wp1%index_acc=icyc
wp1%index_shoot=icyc
call outputTIS(icyc,wp1,syst,pot,dyn,output)
write(output%iutrans,*) icyc, chi*dotx,chiTST*dotx

end subroutine BCstep
!ES----------------------------------------------------------------


!BS------------------------------------------------------------------------
subroutine generate_new_point_on_surface(phasepoint,phasepoint2,transparam,&
                    syst,pot,dyn,timestep,interfaceM,wMDio,wt1,wt2)
use types;use assign_objects;use alloc_objects
use random
use orderparameter
use mdstep
use thermo
use convert
use propagator
implicit none
type(phasepoint_type), intent(inout)::phasepoint,phasepoint2
type(trans_type), intent(in) :: transparam
type(system_type),       intent(inout)::syst
type(potential_type),    intent(inout)::pot  !!not so neat
type(dynamics_type),     intent(in)::dyn
type(timestep_type),     intent(in)::timestep
double precision, intent(in)::interfaceM
type(MDinout_param_type), intent(inout)::wMDio
type(timeslice_type), intent(inout)::wt1,wt2
double precision::sig,op
integer::i
logical::backwardTF
double precision::ran

if (transparam%two_point_method) then
  do i=1,transparam%Ntransrun
    !!call random_number(ran)
    ran=random01()
    if (ran<.5d0) then
      call convert_phasep2tslice(phasepoint,syst,pot,dyn,wt1)
      call Andersen_velocity_change(wt1%phasepoint%phasexv%v,syst%Npart,&
                                  syst%dim,syst%sigma_v,1.d0/timestep%dt) 
      wt2=wt1
      backwardTF=.false.
      call propagate(wt2,syst,pot,dyn,timestep,backwardTF) 
      if (wt2%OPS(1)>interfaceM) then
          phasepoint=wt1%phasepoint
          phasepoint2=wt2%phasepoint
      endif
    else
      call convert_phasep2tslice(phasepoint2,syst,pot,dyn,wt2)
      call Andersen_velocity_change(wt2%phasepoint%phasexv%v,syst%Npart,&
                                 syst%dim,syst%sigma_v,1.d0/timestep%dt)
      wt1=wt2
      backwardTF=.true.
      call propagate(wt1,syst,pot,dyn,timestep,backwardTF)
      if (wt1%OPS(1)<interfaceM) then
          phasepoint=wt1%phasepoint
          phasepoint2=wt2%phasepoint
      endif
    endif
  enddo
else
  select case(pot%POTENTIAL)
    case("HARMOSC")
      phasepoint%phasexv%x(1,1)=interfaceM
      sig=syst%sigma_v(1)
      phasepoint%phasexv%v(1,1)=rangaussian(sig)
    case("DOUBLEWELL")
      phasepoint%phasexv%x(1,1)=interfaceM
      sig=syst%sigma_v(1)
      phasepoint%phasexv%v(1,1)=rangaussian(sig)
    case("PBD")
      pot%POTENTIAL="PBDshift"
      pot%potPBD%surface=interfaceM
      call prepare_MDinout_param(wMDio,phasepoint%phasexv,syst,pot,dyn)
      do i=1,transparam%Ntransrun
        call make_mdstep(phasepoint,syst,pot,dyn,wMDio,timestep)
      enddo 
      pot%POTENTIAL="PBD"
      !shift
      op=orderp(phasepoint%phasexv%x,syst,phasepoint%phasexv%v,pot) 
      phasepoint%phasexv%x(:,1)=phasepoint%phasexv%x(:,1)-op+interfaceM
    case default
      print *,"ERROR generate_new_point_on_surface, pot%POTENTIAL=",&
               pot%POTENTIAL
      stop
  end select
endif


end subroutine generate_new_point_on_surface
!ES----------------------------------------------------------------

end module transstep
!EM------------------------
