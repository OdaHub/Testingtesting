!BM--------------------------------------
module TIS 

contains

  !BS------------------------------------
  subroutine runTIS(startpath,tis_param,syst,pot,dyn,output,timestep,NCYC)
  use types;use assign_objects;use alloc_objects
  use TISstep
  use outTIS
  use restart_module
  use cpmd_subr
  use shell
  use interactive
  implicit none
  type(path_type),         intent(in)::startpath
  type(TIS_type),         intent(inout)::TIS_param  !!inout because i-mode
  type(system_type),       intent(inout)::syst   !! intent(out) because icrash index
  type(potential_type),    intent(in)::pot
  type(dynamics_type),     intent(in)::dyn
  type(output_type),       intent(in)::output
  type(timestep_type),     intent(in)::timestep
  integer,                 intent(in)::NCYC
  integer::i,N,d,NT,NOPS,NX,NWANNIER,dwc
  type(path_type)::path,trial
  type(path_type)::wp1,wp2  !workspace
  type(timeslice_type):: wt1,wt2
  integer::i0,icyc
  character(LEN=LSTR)::Efile
  logical::revab,TISMOLEXIT,TIS_IMODE



  N=startpath%N
  d=startpath%d
  NT=startpath%NT
  NOPS=startpath%NOPS
  NX=startpath%NX
  NWANNIER=startpath%NWANNIER
  dwc=startpath%dwc
  call alloc(path,N,d,NT,NOPS,NX,NWANNIER,dwc)
  call alloc(trial,N,d,NT,NOPS,NX,NWANNIER,dwc)
  call alloc(wp1,N,d,NT,NOPS,NX,NWANNIER,dwc)
  call alloc(wp2,N,d,NT,NOPS,NX,NWANNIER,dwc)
  call alloc(wt1,N,d,NT,NOPS,NWANNIER,dwc)
  call alloc(wt2,N,d,NT,NOPS,NWANNIER,dwc)
  path=startpath

  icyc=output%icyclestart
  i0=output%icyclestart+1
  if (i0==1) call outputTIS(0,path,syst,pot,dyn,output)
 
  do i=i0,Ncyc
    icyc=i
    print *,"*********  icycle=",icyc,Ncyc
    call make_TISstep(path,tis_param,syst,pot,dyn,timestep,trial,wp1,wp2,wt1,wt2)
    call outputTIS(icyc,trial,syst,pot,dyn,output)
    if (trial%ACCREJ=="ACC")  path=trial
    !remove old RESTART FILES 
    if (pot%POTENTIAL=="EXTERNAL") then
      EFILE=path%timeslices(1)%phasepoint%electrons%ESTRUCFILE
      !!if (syst%icrash >= syst%NCPMD_UNSAVED) 
      call delete_efiles(EFILE,revab,syst%icrash,syst%NCPMD_UNSAVED)
    endif

    inquire(file="TISEXIT",EXIST=TISMOLEXIT)
    if (TISMOLEXIT) then
      print *,"NOTICING TISEXIT FILE: SOFT EXIT"
      call do_shell("rm -r TISEXIT")
      exit
    endif
    inquire(file="TIS_IMODE",EXIST=TIS_IMODE)
    if (TIS_IMODE) then
      print *,"NOTICING TIS_IMODE FILE: INTERACTIVE OPTION"
      call imode(tis_param,syst%masses,syst%Npart) 
    endif

 
    if (modulo(icyc,syst%NCYCLE_RESTART)==0) call makerestartPATH(icyc,& 
        path,output%RESTARTFILE_WRITE,"INFO.RESTART",output%IURES,syst)


  enddo

  call makerestartPATH(icyc,path,output%RESTARTFILE_WRITE,"INFO.RESTART",& 
                       output%IURES,syst)
  
  call dealloc(path)
  call dealloc(trial)
  call dealloc(wp1)
  call dealloc(wp2)
  call dealloc(wt1)
  call dealloc(wt2)
 
  !deallocate(path%timeslices) 
  !deallocate(trial%timeslices)
   
  end subroutine runTIS
  !ES------------------------------------

end module TIS 
!EM--------------------------------------
