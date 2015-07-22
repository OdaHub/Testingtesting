module system
implicit none

contains
!BS------------------------------------------------------
subroutine systemx(string)
use var_analyze
implicit none
character*100::string
open(11,file=fexc,position="append")
  write(11,*) trim(string)
close(11)
end subroutine systemx
!ES------------------------------------------------------

end module system
