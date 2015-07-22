Module ordparanalyze 
contains
!BS---------------------------------------------------------
subroutine ordparout
use var_analyze
implicit none

if ((ORDFILE).AND.(NLord>0)) then
  call plot_ord
  call runavop
  call orddist
endif

end subroutine ordparout
!ES----------------------------------------------------------

!BS-----------------------------------------------------------------
subroutine plot_ord
use var
use var_analyze
use efcol
use gnupf
implicit none
character*100::ofile, xlab,name
character*100::title
integer::NL2,col1,col2
double precision::hl(3)

name="ordpar"
ofile=trim(name)//".tmp"
call sparse(fordp,ofile,NLord,ngrid,NL2)
xlab="t (ns)"
title="order parameter"
col1=1;col2=9
hl(1)=lam0;hl(2)=lam1;hl(3)=lam2
call picture(ofile,col1,col2,name,title,xlab,dx=unit_t*dt,nhl=3,hl=hl)

end subroutine plot_ord
!ES-----------------------------------------------------------------

!BS------------------------------------------------------------------
subroutine runavop
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
name="ordpar"
col=9
call runav(fordp,NLord,col,name,xlab=xlab,dx=dx,ngrid=ngrid,pict="y")

end subroutine runavop
!ES------------------------------------------------------------------

!BS------------------------------------------------------------------
subroutine orddist
use var
use var_analyze
use efcol
implicit none
character*100::name, xlab
integer::col

xlab="(angstr)"
name="ordpar"
call distribution(fordp,NLord,9,name,ngrid,xlab)

end subroutine orddist
!ES-----------------------------------------------------------------

end Module ordparanalyze
