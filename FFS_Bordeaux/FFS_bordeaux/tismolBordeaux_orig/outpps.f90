!BM----------------------
Module outPPS

contains

!BS------------------------------
subroutine outputPPS(icyc,pps_ensemble,syst,pot,dyn,out_param)
use types;use assign_objects;use alloc_objects
use outTIS
implicit none
integer, intent(in)::icyc
type(path_ensemble), intent(in)::pps_ensemble
type(output_type),       intent(inout)::out_param 
type(system_type),       intent(in)::syst
type(potential_type),    intent(in)::pot
type(dynamics_type),     intent(in)::dyn
integer::i,N,d,NT,NOPS,NX,numint,NWANNIER,dwc
type(path_type)::path

N=pps_ensemble%N;d=pps_ensemble%d;NT=pps_ensemble%NT
NOPS=pps_ensemble%NOPS;NX=pps_ensemble%NX;numint=pps_ensemble%numint
NWANNIER=pps_ensemble%NWANNIER
dwc=pps_ensemble%dwc

call alloc(path,N,d,NT,NOPS,NX,NWANNIER,dwc)
do i=1,numint
  out_param%IUEN =pps_ensemble%IUEN(i)
  out_param%IUTRAJ =pps_ensemble%IUTRAJ(i)
  out_param%IUOP=pps_ensemble%IUOP(i)
  out_param%IUPATH=pps_ensemble%IUPATH(i)
  out_param%IUCR=pps_ensemble%IUCR(i)
  path=pps_ensemble%PATHS(i)
  call outputTIS(0,path,syst,pot,dyn,out_param)
enddo
call dealloc(path)


end subroutine outputPPS
!ES-------------------------------

end Module outPPS
!EM----------------------
