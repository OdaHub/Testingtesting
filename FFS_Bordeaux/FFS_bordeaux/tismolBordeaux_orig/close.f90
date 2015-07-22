!BM----------------------------
module close

contains

!BS-------------------------------
subroutine closure
use inputpar
use types;use assign_objects;use alloc_objects
use output_module 
use pathensemble_module
use system_module
use shell
implicit none
integer::i

if (TASK/="PPS") then
    close(output%IUEN)
    close(output%IUTRAJ)
    close(output%IUOP)
    close(output%IUCR)
    close(output%IUCROSSPOINTS)
    if (TASK=="TIS") close(output%IUPATH)
else
  do i=1,pps_set%numint
    close(PPS_SET%IUPATH(i))
    close(PPS_SET%IUEN(i))
    close(PPS_SET%IUTRAJ(i))
    close(PPS_SET%IUOP(i))
    close(PPS_SET%IUCR(i))
    close(PPS_SET%IUCROSSPOINTS(i))
  enddo
endif

if (POTENTIAL=="EXTERNAL") then
  close(syst%IUENCPMD)
  close(syst%IUWANCPMD)
  call do_shell("touch EXIT") 
endif

end subroutine closure
!ES-------------------------------
end module close
!EM-----------------------------
