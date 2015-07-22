!BM-----------------------------------------------------------
! This module allows some interactive options like
! changing the interface positions during the run  
! Interactive option is especially handy for runs
! on a queuing system: Parameters can be changed without having 
! to stop the program and restart (queu) again
!-----------------------------------------------------------------
module interactive
implicit none

CONTAINS

!BS--------------------------------------------------------------
subroutine imode(tisparam,masses,N)
use types;use assign_objects;use alloc_objects
use shell
use read_input_file
implicit none
type(TIS_type),         intent(inout)::TISparam 
integer, intent(in)::N
double precision, intent(in)::masses(N)
logical::wait,notwait
integer::twait
character(LEN=LSTR)::inputfile
double precision::sigdp

inputfile="input.TISMOL"


twait=1
call do_shell("rm -fr TIS_IMODE")
wait=.true.
print *,"ENTERING INTERACTIVE MODE"
print *,"LAST CYCLE HAS BEEN FINISHED"
print *,"CHANGE PARAMETERS IN INPUTFILE"
print *,"THEN TYPE ""touch IMODE_GO"" "
do while(wait)
   !call sleep(twait)
   print *,"call sleep() blocked, please remove !!";stop
   INQUIRE(FILE="IMODE_GO",EXIST=notwait)
   if (notwait) wait=.false. 
enddo

print *,"FOUND IMODE_GO, READING NEW PARAMETERS"
call read_inputparameter(tisparam%INTERFACEL,inputfile,"INTERFACEL")
call read_inputparameter(tisparam%INTERFACEM,inputfile,"INTERFACEM")
call read_inputparameter(tisparam%INTERFACER,inputfile,"INTERFACER")
call read_inputparameter(sigdp,inputfile,"SIGDP")
if (sigdp>=0.d0) then
  TISparam%sigdp_sqrtm=sigdp/sqrt(masses)
  TISparam%aimless=.false.
else
  TISparam%sigdp_sqrtm=0.d0
  TISparam%aimless=.true.
endif

call do_shell("rm -fr IMODE_GO")

end subroutine imode
!ES--------------------------------------------------------------


end module interactive
!EM------------------------------------------------------------
