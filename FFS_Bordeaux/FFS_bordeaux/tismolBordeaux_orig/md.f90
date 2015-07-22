!BM--------------------------------------
module MD

contains

  !BS------------------------------------
  subroutine runMD(startpoint,syst,pot,dyn,output,timestep,NMD)
  use types;use assign_objects;use alloc_objects
  use outMD
  use mdstep
  use restart_module
  use shell
  implicit none
  type(phasepoint_type),   intent(in)::startpoint
  type(system_type),       intent(inout)::syst    !!not so nice. intention out is only because icrah ondex...
  type(potential_type),    intent(in)::pot
  type(dynamics_type),     intent(in)::dyn
  type(output_type),       intent(in)::output
  type(timestep_type),     intent(in)::timestep  !! inout for test
  integer,                 intent(in)::NMD
  integer::i,N,d,NT,NWANNIER,dwc,i0,icyc
  type(phasepoint_type)::phasepoint
  type(MDinout_param_type)::MDinout_param
  logical::TISMOLEXIT

  N=startpoint%N
  d=startpoint%d
  NT=startpoint%NT
  NWANNIER=startpoint%NWANNIER
  dwc=startpoint%dwc
  call alloc(phasepoint,N,d,NT,NWANNIER,dwc) 
  call alloc(MDinout_param,N,d)
  phasepoint=startpoint
  icyc=output%icyclestart
  i0=output%icyclestart+1
    

  print *,"MD RUN"
  
  if (i0==1) call outputMD(0,phasepoint,syst,pot,dyn,output)
  
  call prepare_MDinout_param(MDinout_param,phasepoint%phasexv, &
                                    syst,pot,dyn)
  do i=i0,NMD
    icyc=i
    if (.not.syst%noprint) print *,"*********  icycle=",icyc
    call make_mdstep(phasepoint,syst,pot,dyn,MDinout_param,timestep)
    call outputMD(icyc,phasepoint,syst,pot,dyn,output)

    inquire(file="TISEXIT",EXIST=TISMOLEXIT)
    if (TISMOLEXIT) then
      print *,"NOTICING TISEXIT FILE: SOFT EXIT"
      call do_shell("rm -r TISEXIT")
      exit 
    endif

    if (modulo(icyc,syst%NCYCLE_RESTART)==0) &
       call makerestartMD(icyc,phasepoint,output%RESTARTFILE_WRITE,output%IURES,syst) 

  enddo

  !!if (i0>NMD) i=NMD !no md run took place. Keep old index
  call makerestartMD(icyc,phasepoint,output%RESTARTFILE_WRITE,output%IURES,syst)

  call dealloc(phasepoint)
  call dealloc(MDinout_param)
   
  end subroutine runMD
  !ES------------------------------------

end module MD
!EM--------------------------------------
