!BM-------------------------------------------------
module null
implicit none

CONTAINS
!BS-------------------------------------------------------------------
subroutine make_nullmove(j,PPS_ensemble,syst,pot,dyn,output,icyc,wp1,secpot)
use types;use assign_objects;use alloc_objects
use outTIS
use convert
implicit none
integer,intent(in)::j
type(path_ensemble), intent(inout)::PPS_ensemble
type(system_type),        intent(inout)   ::syst
type(potential_type),     intent(in)   ::pot,secpot
type(dynamics_type),      intent(in)   ::dyn
type(output_type), intent(inout):: output
integer,                  intent(in)   ::icyc
type(path_type), intent(inout)::wp1   !workspace

print *,"make null move ens",j
call set_TIS(j,pps_ensemble,wp1,output)
wp1%MCmove="00"
wp1%index_acc=wp1%index_acc+1
PPS_ensemble%PATHS(j)%index_acc=PPS_ensemble%PATHS(j)%index_acc+1
if ((j==1).AND.(PPS_ensemble%Forcefieldmatching)) then
  call outputTIS(icyc,wp1,syst,secpot,dyn,output)
else
  call outputTIS(icyc,wp1,syst,pot,dyn,output)
endif

end subroutine make_nullmove
!ES-------------------------------------------------------------------


end module null
!EM----------------------------------------------------
