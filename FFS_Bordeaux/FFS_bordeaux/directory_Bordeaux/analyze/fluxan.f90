module flux

contains
!BS-------------------------------------------
subroutine fluxan
use var_analyze
implicit none

if (CROSSFILE) then
  call calc_flux
  call calcrunav
  call err_flux
endif

end subroutine fluxan
!ES-------------------------------------------

!BS-------------------------------------------
subroutine calc_flux
use var
use var_analyze
implicit none
integer::it,intfac,NHB,i,tBs,tBe,ic,xx,tend
character::char
logical::HB
double precision::mm,dtm1,av,mini,maxi

allocate(apcro1(NLcross))
apcro1=0
NHB=0;Ncross=0
HB=.false.
tend=NLMD*skip

!Count positive crossings through 1 and the number of points it is in state A
open(1,file=filecross)
  do i=1,NLcross
    read(1,*) it,intfac, char
    if ((intfac==1).AND.(char=="+")) then
       Ncross=Ncross+1
       apcro1(Ncross)=it-NHB
    endif
    if ((intfac==2).AND.(char=="+")) then 
       HB=.true. 
       tBs=it+1 !the start of a time period the system is in overall state B
    endif
    if ((intfac==2).AND.(char=="-")) then
       HB=.false. 
       tBe=it
       NHB=NHB+tBe-tBs+1 
    endif
  enddo
close(1)

!the system might be in B at the end of the run
if (HB) NHB=NHB+tend-tBs

NHA=tend-NHB
rhoHA=1.d0*NHA/(NHA+NHB)
rhoHB=1.d0-rhoHA

Flux= Ncross/(NHA*dt*unit_t)
print *,"NTOT,NHA, NHB",tend,NHA,NHB
print *,"Ncross", Ncross
print *,"Flux (ns^-1)",Flux

!create new crossing file in a standard format that allows error analysis
dtm1=1.d0/(dt*unit_t*skip)    ! one over 'big delta T'=dt*skip
NC2=int(1.d0*NHA/skip)        ! number of timeslices new file
ic=1
fcross2="cross2.tmp"

xx=apcro1(ic) !the first next crossing
open(1,file=fcross2)
  do i=1,Nc2
    mm=0.d0
    do while ((xx < i*skip).AND.(ic<=Ncross))  !if i*skip is beyond the crossing point
      mm=mm+dtm1                               !mm can be higher than dtm1 is more
      ic=ic+1                                  !crossings occur in time interval 'skip'
      if (ic<=Ncross) xx=apcro1(ic)
    enddo
    write(1,*) i*skip,mm
  enddo
close(1)

end subroutine calc_flux
!ES--------------------------------------------

!BS--------------------------------------------
subroutine calcrunav
use var
use var_analyze
use efcol
implicit none
integer::col
character*100::name
character*100::xlab
double precision::dx

dx=unit_t*dt*skip
xlab="t(ns)"
name="flux"
col=2
call runav(fcross2,NC2,col,name,xlab=xlab,dx=dx,ngrid=ngrid,pict="y")

end subroutine calcrunav
!ES---------------------------------------------

!BS-------------------------------------------
subroutine err_flux
use var
use var_analyze
use efcol
implicit none
character*100::name

name="flux"
call block_error(fcross2,NC2,2,name,lbmax,rerr=rerrflux,corr=Ncfl,stdev=stdfl)
errflux=rerrflux*flux

end subroutine err_flux
!ES--------------------------------------------

end module flux 
