!BM--------------------------------------
module PPS

contains

  !BS------------------------------------
  subroutine runPPS(PPS_ensemble,tisparam,syst,pot,dyn,output,timestep,& 
                    NCYC,secpot)
  use types;use assign_objects;use alloc_objects
  use PPSstep
  use outPPS
  use restart_module
  use shell
  implicit none
  type(path_ensemble), intent(inout)::PPS_ensemble
  type(TIS_type), intent(inout)::tisparam
  type(system_type),   intent(inout)::syst !!intent(out) because icrash index
  type(potential_type),intent(in)::pot,secpot
  type(dynamics_type), intent(in)::dyn
  type(output_type),   intent(inout)::output
  type(timestep_type), intent(in)::timestep
  integer,             intent(in)::NCYC
  integer::i,N,d,NT,NOPS,NX,NWANNIER,dwc
  type(path_type)::path,trial,wp1,wp2,wp3
  type(timeslice_type)::wt1,wt2
  integer::i0,icyc
  logical::TISMOLEXIT
  
  N=pps_ensemble%N;d=pps_ensemble%d;NT=pps_ensemble%NT
  NOPS=pps_ensemble%NOPS;NX=pps_ensemble%NX;NWANNIER=pps_ensemble%NWANNIER
  dwc=pps_ensemble%dwc

  call alloc(tisparam,n)
  call alloc(output,output%ncrossplanes)
  call alloc(path ,N,d,NT,NOPS,NX,NWANNIER,dwc)
  call alloc(trial,N,d,NT,NOPS,NX,NWANNIER,dwc)
  call alloc(wp1  ,N,d,NT,NOPS,NX,NWANNIER,dwc)
  call alloc(wp2  ,N,d,NT,NOPS,NX,NWANNIER,dwc)
  call alloc(wp3  ,N,d,NT,NOPS,NX,NWANNIER,dwc)
  call alloc(wt1,N,d,NT,NOPS,NWANNIER,dwc)
  call alloc(wt2,N,d,NT,NOPS,NWANNIER,dwc)

  icyc=output%icyclestart
  i0=output%icyclestart+1
  if (i0==1) call outputPPS(0,PPS_ensemble,syst,pot,dyn,output)
 
 
  do i=i0,Ncyc
    icyc=i
    print *,"*********  icycle=",icyc
    call make_PPSstep(PPS_ensemble,tisparam,syst,pot,dyn,output,timestep,icyc,path,trial,wp1,wp2,wp3,wt1,wt2,secpot) 
    !this includes output of trials in each path-directory

    inquire(file="TISEXIT",EXIST=TISMOLEXIT)
    if (TISMOLEXIT) then
      print *,"NOTICING TISEXIT FILE: SOFT EXIT"
      call do_shell("rm -r TISEXIT")
      exit
    endif
    
    if (modulo(icyc,syst%NCYCLE_RESTART)==0) call makeRestartPPS(icyc,pps_ensemble,wp1,syst)

  enddo
  print *,"############## end of PPS loop"

  call makeRestartPPS(icyc,pps_ensemble,wp1,syst)

  call dealloc(tisparam)
  call dealloc(output)
  call dealloc(path)
  call dealloc(trial)
  call dealloc(wp1)
  call dealloc(wp2)
  call dealloc(wp3)
  call dealloc(wt1)
  call dealloc(wt2)


  end subroutine runPPS
  !ES------------------------------------

end module PPS 
!EM--------------------------------------
