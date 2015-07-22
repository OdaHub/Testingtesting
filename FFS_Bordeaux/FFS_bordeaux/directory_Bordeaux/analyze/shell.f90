!BM----------------------------------------
module shell
implicit none

contains

!BS-----------------------------------------
subroutine do_shell(command,noprint)
implicit none
character(LEN=*), intent(in)::command
logical, optional, intent(in)::noprint

if (.not.present(noprint)) then
  write(*,'(A)') trim(adjustl("SHELL:"//trim(command)))
else
  if (.not.noprint) write(*,'(A)') trim(adjustl("SHELL:"//trim(command)))
endif
call system(command)
!print *,"call system()  blocked, please remove !! "

end subroutine do_shell
!ES-----------------------------------------

end module shell
!EM----------------------------------------
