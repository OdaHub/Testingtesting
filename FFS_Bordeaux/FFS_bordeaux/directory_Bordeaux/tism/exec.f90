!BM----------------------------------
Module exec

contains

!BS----------------------------------
subroutine execute
use inputpar, only: TASK, Ncyc
use system_module
use timestep_module
use phase_module
use path_module
use pot_module
use dyn_module
use output_module
use tis_module
use trans_module
use ffs_module
use pathensemble_module
use integration_module
use MD
use TIS
use PPS
use trans
use FFS
use numericinteg
implicit none



select case(TASK)
  case("MD")
    call runMD(startpoint,syst,pot,dyn,output,timestep,Ncyc)
  case("TIS")
    call runTIS(startpath,tis_param,syst,pot,dyn,output,timestep,Ncyc)
  case("PPS")
    call runPPS(PPS_SET,tis_param,syst,pot,dyn,output,timestep,Ncyc,secpot)
  case("TRANSMISSION")
    call runtrans(startpoint,trans_param,syst,pot,dyn,output,timestep,Ncyc)
  case("NUMERICINTEG")
    call runNI(numinteg_param,syst,pot)
  case("FFS")
    call runffs(startpoint,FFS_param,syst,pot,dyn,output,timestep,Ncyc)
  case default
    print *,"ERROR, SUBROUTINE EXECUTE:TASK=",TASK
    stop
end select

end subroutine execute
!ES----------------------------------

end Module exec
!EM----------------------------------
