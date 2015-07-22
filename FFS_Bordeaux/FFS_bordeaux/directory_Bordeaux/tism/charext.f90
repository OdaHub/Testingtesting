!BM----------------------
module charext 
implicit none

CONTAINS

  !BS---------------------------------------------------------------------
  ! This subroutine returns an character*3-array with (000),001,002, etc
  ! Per default, the first character is 001 and put on array position 1
  ! startchar000=.true. makes the first extension 000 instead of 001
  ! startpos0=.true. make the first array position 0 instead of 1
  !-----------------------------------------------------------------------
  subroutine set_char_extensions(EXT,N) !,startchar000) !!,startpos0)
  implicit none
  integer, intent(in)::N
  character*4, intent(out)::EXT(-N:N)
  !!logical, optional, intent(in)::startchar000 !!,startpos0
  integer::i !!,ichar,ipos,charshift,arrayshift
  character*3::dum

  
  !!charshift=0;arrayshift=0
  !!if (present(startchar000)) charshift=1
  !if (present(startchar000)) arrayshift=1
  EXT(0)="+000"
  do i=1,N
    if (i <10) then
      write(dum,'(A2,i1)')"00",i
    else if (i <100) then
      write(dum,'(A1,i2)')"0",i
    else if (i <1000) then
      write(dum,'(i3)')i
    else
      print *,"ERROR set_char_extensions: N=",N
      stop
    endif
    EXT(i)="+"//dum
    EXT(-i)="-"//dum
    !!print *,"gf",ipos,EXT(ipos),N
  enddo
!! do i=0,N
!!    ichar=i-charshift
!!    ipos=i-arrayshift
!!    if (ichar <10) then
!!      write(dum,'(A2,i1)')"00",ichar
!!    else if (ichar <100) then
!!      write(dum,'(A1,i2)')"0",ichar
!!    else if (ichar <1000) then
!!      write(dum,'(i3)')ichar
!!    else
!!      print *,"ERROR set_char_extensions: N=",N
!!      stop
!!    endif
!!    EXT(ipos)=dum
!!    !!print *,"gf",ipos,EXT(ipos),N
!!  enddo

  
  end subroutine set_char_extensions 
  !EF----------------------------------------------------------

end module charext 
!EM----------------------

