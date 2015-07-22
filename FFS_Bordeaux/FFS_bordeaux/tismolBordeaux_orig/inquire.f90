!BM--------------------------------------------
module inquire

contains

!BS------------------------------------------------------
subroutine check_datafile(dir,subdir,filename,fullfile,NL,TASK,keyword)
use stringlengths
implicit none
character(LEN=*), intent(in)::dir,subdir,filename,task,keyword
character(LEN=*), intent(out)::fullfile
integer, intent(out)::NL	
logical::exist

NL=0
if (TASK==keyword) then  !first check if data-files are in main dir         
  fullfile=trim(dir)//"/"//trim(FILENAME)
  INQUIRE(FILE=fullfile,EXIST=exist)
  if (exist) call countlines(fullfile,NL)
endif

if (NL==0) then
  fullfile=trim(dir)//"/"//trim(subdir)//"/"//trim(FILENAME)
  INQUIRE(FILE=fullfile,EXIST=exist)
  if (exist) call countlines(fullfile,NL)
endif

end subroutine check_datafile
!ES-------------------------------------------------------

!BS--------------------------------------------------------------
! This subroutine counts the lines of a file
!----------------------------------------------------------------
subroutine countlines(filenm,N)
implicit none
character(LEN=*), intent(in)::filenm
integer, intent(out)::N
character::char
integer::status

  N=0
  open(1,file=filenm,status="old")
  do
    READ(1,'(a1)',iostat=status) char
    if (status /= 0 ) exit
    N=N+1
  enddo
  close(1)

end subroutine countlines
!ES--------------------------------------------------------------

!BF---------------------------------------------------------------
function nonempty(file)
implicit none
character(LEN=*), intent(in)::file
character::char
logical::nonempty
logical::exist
integer::status

   nonempty=.false.
   INQUIRE(FILE=file,exist=exist)

   if (exist) then
     open(1,file=file,status="old")
     READ(1,'(a1)',iostat=status) char
     if (status == 0 ) nonempty=.true.
     close(1)
   endif

end function nonempty
!EF---------------------------------------------------------------


end module inquire
!EM--------------------------------------------
