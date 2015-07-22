module GAMMA_func
implicit none

Contains

!--------------------------------------------------------------
! This function is equivalent to GAMMA(N/2) for the case
! N equal to an integer
!-------------------------------------------------------------
function GAMMA_N2(N)
use var, only:pi
implicit none
integer, intent(in)::N
double precision::GAMMA_N2
integer::N2

 if (modulo(N,2)==0) then
   N2=N/2-1
   GAMMA_N2=FAC(N2)
 else
   N2=(N-1)/2
   GAMMA_N2=sqrt(pi)*doublefac(N-2)/(2**N2) 
 endif
 
end function GAMMA_N2
!---------------------------------------------------------------

!---------------------------------------------
function FAC(N)
implicit none
integer, intent(in)::N
double precision::FAC
integer::i

FAC=1.
do i=1,N
  FAC=FAC*i
enddo
end function FAC
!--------------------------------------------

!---------------------------------------------
function doubleFAC(N)
implicit none
integer, intent(in)::N
double precision::doubleFAC
integer::i

doubleFAC=1.
do i=1,N,2
  doubleFAC=doubleFAC*i
enddo
end function doubleFAC
!--------------------------------------------


end module GAMMA_func

