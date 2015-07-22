!BM--------------------------------------------
module convert 

contains

!BS------------------------------------------------------
subroutine convert_phasep2tslice(phasep,syst,pot,dyn,tslice)
use types;use assign_objects;use alloc_objects
use mdstep
use orderparameter
implicit none
type(phasepoint_type), intent(in)::phasep
type(system_type), intent(in)::syst
type(potential_type), intent(in)::pot
type(dynamics_type), intent(in)::dyn
type(timeslice_type), intent(out)::tslice

tslice%phasepoint=phasep
print *,size(tslice%OPS)
tslice%OPS(1)=orderp(phasep%phasexv%x,syst,phasep%phasexv%v,pot)
call prepare_MDinout_param(tslice%MDinout_param,tslice%phasepoint%phasexv,&
                              syst,pot,dyn)


end subroutine  convert_phasep2tslice
!ES-------------------------------------------------------

!BS-----------------------------------------------------------------
subroutine set_TIS(i,pps_ensemble,path,output,pt,pot,secpot,tisparam)
use types;use assign_objects;use alloc_objects
implicit none
integer, intent(in)::i
type(path_ensemble), intent(in)::pps_ensemble
type(path_type), optional, intent(out)::path
type(output_type), optional, intent(inout)::output
type(TIS_type), optional, intent(inout)::tisparam
type(potential_type), optional, intent(inout)::pt
type(potential_type), optional, intent(in)::pot,secpot
integer::n,numint,j

 n=pps_ensemble%n
 numint=pps_ensemble%numint

 if (present(path)) path=pps_ensemble%PATHS(i)

 if (present(output)) then
   output%IUEN =pps_ensemble%IUEN(i)
   output%IUTRAJ =pps_ensemble%IUTRAJ(i)
   output%IUOP=pps_ensemble%IUOP(i)
   output%IUPATH=pps_ensemble%IUPATH(i)
   output%IUCR=pps_ensemble%IUCR(i)
   output%IUCROSSPOINTS=pps_ensemble%IUCROSSPOINTS(i)
 endif


 if (present(tisparam)) then
   if (i==1) then
     tisparam%startcondition="R"
     tisparam%INTERFACEL=-99999999.d0
     tisparam%INTERFACEM=pps_ensemble%PPS_interfaces(1)
     tisparam%INTERFACER=pps_ensemble%PPS_interfaces(1)
   else
     tisparam%startcondition="L"
     tisparam%INTERFACEL=pps_ensemble%PPS_interfaces(1)
     tisparam%INTERFACEM=pps_ensemble%PPS_interfaces(i-1)
     tisparam%INTERFACER=pps_ensemble%PPS_interfaces(NUMINT)
   endif
   do j=1,n
     tisparam%sigdp_sqrtm(j)=pps_ensemble%pps_sigdp_sqrtm(i,j)
   enddo
   tisparam%aimless=pps_ensemble%pps_aimless(i)

 endif

 if (present(pt)) then
   if ((i==1).AND.(pps_ensemble%forcefieldmatching)) then
     pt=secpot
   else
     pt=pot
   endif
 endif

end subroutine set_TIS
!ES-----------------------------------------------------------------




end module convert 
!EM--------------------------------------------
