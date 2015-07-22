!BM--------------------------------------
module FFS 

contains

  !BS------------------------------------
  subroutine runffs(startpoint,ffsparam,syst,pot,dyn,output,timestep,NCYC)
  use types;use assign_objects;use alloc_objects
  !use ffsstep
  use mdstep
  use orderparameter
  use convert
  use thermo
  use propagator
  use shell
  implicit none
  type(phasepoint_type),   intent(in)::startpoint 
  type(ffs_type),          intent(in):: ffsparam
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
  double precision::ran
  integer::iselect
  
  TISMOLEXIT=.false.
  N=startpoint%N
  d=startpoint%d
  NT=startpoint%NT
  NWANNIER=startpoint%NWANNIER
  dwc=startpoint%dwc
  NX=ffsparam%NX
  NOPS=ffsparam%NOPS
  call alloc(phasepoint,N,d,NT,NWANNIER,dwc)
  call alloc(phasepoint2,N,d,NT,NWANNIER,dwc)
  call alloc(wp1,N,d,NT,NOPS,NX,NWANNIER,dwc)
  call alloc(wp2,N,d,NT,NOPS,NX,NWANNIER,dwc)
  call alloc(wp3,N,d,NT,NOPS,NX,NWANNIER,dwc)
  call alloc(wt1,N,d,NT,NOPS,NWANNIER,dwc)
  call alloc(wt2,N,d,NT,NOPS,NWANNIER,dwc)
  call alloc(wMDio,N,d)
  interfaceL=ffsparam%interfaceL
  interfaceM=ffsparam%interfaceM
  interfaceR=ffsparam%interfaceR
  icyc=output%icyclestart
  i0=output%icyclestart+1
  call convert_phasep2tslice(phasepoint,syst,pot,dyn,wt1) 
  phasepoint2=wt2%phasepoint
  do i=i0,Ncyc
    icyc=i
    call random_number(ran)
    iselect=int(ffsparam%nstartpoints*ran)+1
    phasepoint%phasexv%x=ffsparam%startpoints(iselect,1) 
    phasepoint%phasexv%v=ffsparam%startpoints(iselect,2)
    !call make_transstep(icyc,phasepoint,phasepoint2,syst,pot,dyn,output,&
    !                    timestep,interfaceL,interfaceM,interfaceR,&
    !                    transalgorithm,two_point_method,& 
    !                    wp1,wp2,wp3,wt1,wt2)

    if (TISMOLEXIT) then
      print *,"NOTICING TISEXIT FILE: SOFT EXIT"
      call do_shell("rm -r TISEXIT")
      exit
    endif

  enddo

  print *,"RESTART OPTION FOR FFS NOT YET IMPLEMENTED"
 
  call dealloc(phasepoint) 
  call dealloc(phasepoint2)
  call dealloc(wp1)
  call dealloc(wp2)
  call dealloc(wp3)
  call dealloc(wt1)
  call dealloc(wt2)
  call dealloc(wMDio)
   
  end subroutine runffs
  !ES------------------------------------

end module FFS 
!EM--------------------------------------
