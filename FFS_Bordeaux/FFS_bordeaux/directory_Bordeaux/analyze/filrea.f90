Module filrea
implicit none

CONTAINS

!BS--------------------------------------------------------------
! This subroutine counts the lines of a file
!----------------------------------------------------------------
subroutine countlines(filenm,N)
implicit none
character*100, intent(in)::filenm
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

!BS---------------------------------------------------------------------
! This subroutine searches if it can find the character-sequence "string"
! in the file "filenm" between the lines UL and LL. If found, it returns
! the corresponding line number iL. Else it returns iL=0.
!-------------------------------------------------------------------------
subroutine check_string(filenm,UL,LL,string,iL)
implicit none
character*100, intent(in)::filenm,string
integer, intent(in)::UL,LL
integer, intent(out)::iL
character*100::LINE
integer::i

iL=0
open(1,file=filenm)
  do i=1,UL-1
    READ(1,'(A80)') LINE
  enddo
 do i=UL,LL
    !READ(1,*) LINE
    READ(1,'(A80)') LINE
    LINE=ADJUSTL(LINE)
    if (INDEX(LINE,trim(string))==1) then
       iL=i
       exit
    endif
  enddo

close(1)
end subroutine check_string
!ES---------------------------------------------------------------------



!BS---------------------------------------------------------------------
! This subroutine searches in the file "filenm" for the character
! sequence "string" followed by a value. This value is considered either 
! integer (if typ="int") or double precision (if typ="dp").
! If "string" is present then, depending on "typ", either integ or
! dp will be replaced by the value after "string" in the file.
! If "string" is not present, dp and integ are returned unchanged.
!-----------------------------------------------------------------------
subroutine getval(filenm,string,UL,LL,typ,dp,integ)
implicit none
character*100, intent(in)::filenm,string
integer, intent(in)::UL,LL
character*3::typ
double precision, intent(inout)::dp
integer, intent(inout)::integ
integer::iL

call check_string(filenm,UL,LL,string,iL)
if (iL/=0) call read_fileln(filenm,iL,typ,dp,integ)

if (typ=="dp")  print *,trim(string),"=",dp
if (typ=="int") print *,trim(string),"=", integ

end subroutine getval
!ES---------------------------------------------------------------------


!BS---------------------------------------------------------------------
! This subroutine reads the line iL in the file "filenm". This line
! should consist off first a string, followed by an integer or
! double precision number. It assumes one of the two depending
! on typ="int" or typ="dp" and it returns this value as integ or dp.
!------------------------------------------------------------------------
subroutine read_fileln(filenm,iL,typ,dp,integ)
implicit none
character*100, intent(in)::filenm
character*3, intent(in)::typ
integer, intent(in)::iL
double precision, intent(out)::dp
integer, intent(out)::integ
character*5::LINE
integer::i

open(1,file=filenm)
  do i=1,iL-1
    READ(1,'(A5)') LINE
  enddo
  if (typ=="dp") then
     READ(1,*) LINE, dp
     integ=0
  else if (typ=="int") then
     READ(1,*) LINE,integ
     dp=0.d0
  endif
close(1)
end subroutine read_fileln
!ES---------------------------------------------------------------------

!BS---------------------------------------------------------------------
! This subroutine reads the first string of line iL in the file "filenm" 
!-----------------------------------------------------------------------
subroutine read_fileln_string(filenm,iL,string)
implicit none
character*100, intent(in)::filenm
integer, intent(in)::iL
character*100, intent(out)::string
character*5::LINE
integer::i

open(1,file=filenm)
  do i=1,iL-1
    READ(1,'(A5)') LINE
  enddo
  READ(1,*) string 
close(1)
end subroutine read_fileln_string
!ES---------------------------------------------------------------------

!BS---------------------------------------------------------------------
! This subroutine reads the complete line iL in the file "filenm"
!-----------------------------------------------------------------------
subroutine read_fileln_string2(filenm,iL,string)
implicit none
character*100, intent(in)::filenm
integer, intent(in)::iL
character*100, intent(out)::string
character*5::LINE
integer::i

open(1,file=filenm)
  do i=1,iL-1
    READ(1,'(A5)') LINE
  enddo
  READ(1,'(A100)') string
close(1)
end subroutine read_fileln_string2
!ES---------------------------------------------------------------------

END MODULE filrea
