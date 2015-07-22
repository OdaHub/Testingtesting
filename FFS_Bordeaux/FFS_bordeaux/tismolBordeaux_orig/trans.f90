!BM--------------------------------------
module trans 

contains

  !BS------------------------------------
  subroutine runtrans(startpoint,transparam,syst,pot,dyn,output,timestep,NCYC)
  use types;use assign_objects;use alloc_objects
  use transstep
  use mdstep
  use orderparameter
  use convert
  use thermo
  use propagator
  use shell
  implicit none
  type(phasepoint_type),   intent(in)::startpoint 
  type(trans_type):: transparam
  type(system_type),       intent(inout)::syst
  type(potential_type),    intent(inout)::pot  !!not so neat
  type(dynamics_type),     intent(in)::dyn
  type(output_type),       intent(in)::output
  type(timestep_type),     intent(in)::timestep
  integer,                 intent(in)::NCYC
  integer::i,N,d,NT,NOPS,NX,NWANNIER,dwc
  type(phasepoint_type)::phasepoint,phasepoint2
  type(path_type)::wp1,wp2,wp3
  type(timeslice_type)::wt1,wt2
  integer::i0,icyc
  double precision::interfaceL,interfaceM,interfaceR
  type(MDinout_param_type)::wMDio
  character(LEN=MSTR):: TRANSALGORITHM
  logical::two_point_method,backwardTF,TISMOLEXIT
  
  TISMOLEXIT=.false.
  N=startpoint%N
  d=startpoint%d
  NT=startpoint%NT
  NWANNIER=startpoint%NWANNIER
  dwc=startpoint%dwc
  NX=transparam%NX
  NOPS=transparam%NOPS
  call alloc(phasepoint,N,d,NT,NWANNIER,dwc)
  call alloc(phasepoint2,N,d,NT,NWANNIER,dwc)
  call alloc(wp1,N,d,NT,NOPS,NX,NWANNIER,dwc)
  call alloc(wp2,N,d,NT,NOPS,NX,NWANNIER,dwc)
  call alloc(wp3,N,d,NT,NOPS,NX,NWANNIER,dwc)
  call alloc(wt1,N,d,NT,NOPS,NWANNIER,dwc)
  call alloc(wt2,N,d,NT,NOPS,NWANNIER,dwc)
  call alloc(wMDio,N,d)
  interfaceL=transparam%interfaceL
  interfaceM=transparam%interfaceM
  interfaceR=transparam%interfaceR
  transalgorithm=transparam%transalgorithm
  two_point_method=transparam%two_point_method

  icyc=output%icyclestart
  i0=output%icyclestart+1
  phasepoint=startpoint
  call convert_phasep2tslice(phasepoint,syst,pot,dyn,wt1) 
  do 
    wt2=wt1
    call Andersen_velocity_change(wt2%phasepoint%phasexv%v,syst%Npart,&
                                syst%dim,syst%sigma_v,1.d0/timestep%dt)
    backwardTF=.false.
    call propagate(wt2,syst,pot,dyn,timestep,backwardTF)
    if (wt2%OPS(1)>interfaceM) exit
  enddo
  phasepoint2=wt2%phasepoint

  if (i0==1) call make_transstep(0,phasepoint,phasepoint2,syst,pot,dyn, &
                                 output,timestep,interfaceL,interfaceM, &
                                 interfaceR,transalgorithm,two_point_method,&
                                 wp1,wp2,wp3,wt1,wt2) 
  do i=i0,Ncyc
    icyc=i
    call generate_new_point_on_surface(phasepoint,phasepoint2,transparam, &
                                       syst,pot,dyn,timestep,interfaceM,&
                                       wMDIO,wt1,wt2)
    call make_transstep(icyc,phasepoint,phasepoint2,syst,pot,dyn,output,&
                        timestep,interfaceL,interfaceM,interfaceR,&
                        transalgorithm,two_point_method,& 
                        wp1,wp2,wp3,wt1,wt2)

    if (TISMOLEXIT) then
      print *,"NOTICING TISEXIT FILE: SOFT EXIT"
      call do_shell("rm -r TISEXIT")
      exit
    endif

  enddo

  print *,"RESTART OPTION FOR TRANSMISSION COEFFICIENTS NOT YET IMPLEMENTED"
  !call makeRestartTRANS(icyc,...)
 
  call dealloc(phasepoint) 
  call dealloc(phasepoint2)
  call dealloc(wp1)
  call dealloc(wp2)
  call dealloc(wp3)
  call dealloc(wt1)
  call dealloc(wt2)
  call dealloc(wMDio)
   
  end subroutine runtrans
  !ES------------------------------------

end module trans 
!EM--------------------------------------
