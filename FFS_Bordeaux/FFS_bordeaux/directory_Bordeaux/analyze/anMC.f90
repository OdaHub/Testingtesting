module anMC
implicit none

contains
!BS------------------------------------
subroutine effMCan(iext)
implicit none
integer, intent(in)::iext
!  call makebarfig(iext)
!  call makedistMC(iext)
  call distshoot(iext)
!  call distdE(iext)
end subroutine effMCan
!ES------------------------------------

!BS--------------------------------------
subroutine makebarfig(iext)
use var_analyze
use gnupf
implicit none
integer, intent(in)::iext
character::dum(10)
character*2::mc
character*3::AR
integer::Ntr,Nsh,NtrACC,NtrBWI,NshACC,NshMCR,NshBWI,NshNCR,NshBTL
integer::NshFTL,NshBML,NshFML
integer::i
double precision::shacc
character*100::name,xlab,dfile
character*100::title

Ntr=0;Nsh=0;NtrACC=0
NtrBWI=0;NshACC=0;NshMCR=0;NshBWI=0;NshNCR=0
NshBTL=0;NshFTL=0;NshBML=0;NshFML=0

open(1,file=fpath)
  do i=1,Nlpath
    read(1,*) dum(1:2),mc,AR
    if (mc=="tr") then
      Ntr=Ntr+1
      selectcase(AR)
        case("ACC")
          NtrACC=NtrACC+1
        case("BWI")
          NtrBWI=NtrBWI+1
        case default
          print *,"found invalid combination: mc=",mc," AR=",AR
          stop
      end select
    else if (mc=="sh") then
      Nsh=Nsh+1
      select case (AR)
        case("ACC")
          NshACC=NshACC+1
        case("MCR")
          NshMCR=NshMCR+1
        case("BWI")
          NshBWI=NshBWI+1
        case("NCR")
          NshNCR=NshNCR+1
        case("BTL")
          NshBTL=NshBTL+1
        case("FTL")
          NshFTL=NshFTL+1
        case("BML")
          NshBML=NshBML+1
        case("FML")
          NshFML=NshFML+1
        case default
          print *,"found AR=",AR
          stop
      end select
    else
      print *,"makebarfig-error:mc=",mc," in file",fpath,i
      stop
    endif
  enddo
close(1)

name="effMC"//EXT(iext)
dfile=trim(name)//".tmp"
open(1,file=dfile)
  write(1,*)  1, Ntr
  write(1,*)  2, Nsh
  write(1,*)  4, NtrACC
  write(1,*)  5, NtrBWI
  write(1,*)  7, NshACC
  write(1,*)  8, NshMCR
  write(1,*)  9, NshBWI
  write(1,*) 10, NshNCR
  write(1,*) 11, NshBTL
  write(1,*) 12, NshFTL
  write(1,*) 13, NshBML
  write(1,*) 14, NshFML
close(1)
title= "1=tr 2=sh,tr:4=ACC,5=BWI,sh:7=ACC,8=MCR,9=BWI,10=NCR,11=BTL,12=FTL"
shacc=(1.d0*NshACC)/Nsh
write(xlab,'(A27,f10.4)') "13=BML,14=FML,  %shoot acc:",shacc
call barfig(dfile,1,2,name,xlab,title)

end subroutine makebarfig
!ES--------------------------------------

!BS-------------------------------------------
subroutine makedistMC(iext)
use var
use var_analyze
use efcol
use Gfunc
implicit none
integer, intent(in)::iext
character*100::name, xlab

!tol=0.01d0
!relerr=rerrpcr(iext)
!Ncpu=( relerr**2*Nlpath*av(1) )/(tol**2)

xlab="length"
name="LMC"//EXT(iext)
call distribution(fpath,NLpath,6,name,ngrid,xlab,avo=LMC(iext))

end subroutine makedistMC
!ES----------------------------------------------

!BS---------------------------------------------------------
subroutine distshoot(iext)
use var_analyze
use efcol
implicit none
integer, intent(in)::iext
character*100::ofile1,ofile2,ofile3,ofile4,xlab,name
character::dum(10)
integer::i,NLo1,NLo2,NLo3,NLo4 
character*2::MC
character*3::AR
double precision::lsh 

ofile1="shoot"//ext(iext)//".tmp"
ofile2="shootACC"//ext(iext)//".tmp"
ofile3="shootBWI"//ext(iext)//".tmp"
ofile4="shootNCR"//ext(iext)//".tmp"

NLo1=0;NLo2=0;NLo3=0;NLo4=0

open(5,file=fpath)
open(1,file=ofile1)
open(2,file=ofile2)
open(3,file=ofile3)
open(4,file=ofile4)
  do i=1,Nlpath
     read(5,*) dum(1:2),MC,AR, dum(1:6),lsh
     if (MC=="sh") then
        NLo1=NLo1+1
        write(1,*)  NLo1,lsh
        if (AR=="ACC") then
           NLo2=NLo2+1
           write(2,*) NLo2,lsh
        endif
        if (AR=="BWI") then
           NLo3=NLo3+1
           write(3,*) NLo3,lsh
        endif
        if (AR=="NCR") then
           NLo4=NLo4+1
           write(4,*) NLo3,lsh
        endif
     endif
  enddo
close(4)
close(3)
close(2)
close(1)
close(5)


ofile1="shoot"//ext(iext)//".tmp"
ofile2="shootACC"//ext(iext)//".tmp"
ofile3="shootBWI"//ext(iext)//".tmp"
ofile4="shootNCR"//ext(iext)//".tmp"


xlab="lambda"
if (NLo1>0) then
  name="shoot"//EXT(iext)
  call distribution(ofile1,NLo1,2,name,ngrid,xlab)
endif
if (NLo2>1) then
  name="shootACC"//EXT(iext)
  call distribution(ofile2,NLo2,2,name,ngrid,xlab)
endif
if (NLo3>1) then
  name="shootBWI"//EXT(iext)
  call distribution(ofile3,NLo3,2,name,ngrid,xlab)
endif
if (NLo4>1) then
  name="shootNCR"//EXT(iext)
  call distribution(ofile4,NLo4,2,name,ngrid,xlab)
endif

end subroutine distshoot
!ES---------------------------------------------------------

!BS-----------------------------------------------------------
subroutine distdE(iext)
use var
use var_analyze
use efcol
use calcpcross
use gnupf
implicit none
integer, intent(in)::iext
double precision::Eo,En,dE
integer::i,NLo1,NLo2 
character*100::ofile1,ofile2,xlab,name1,name2
character*100::title,ax(4)
character::dum(10)
character*2::MC
character*3::AR
name1="momchan"//ext(iext)
name2="mochrej"//ext(iext)
ofile1=trim(name1)//".tmp"
ofile2=trim(name2)//".tmp"
NLo1=0
open(1,file=fpath)
open(2,file=ofile1)
  do i=1,Nlpath
     read(1,*) dum(1:2),MC,AR, dum(1:5),En
     if (i==1) then 
       Eo=En
       cycle
     endif
     if (MC=="sh") then 
       NLo1=NLo1+1
       dE=En-Eo
       write(2,*) NLo1,dE,AR
       if (AR=="ACC") Eo=En
       if ((AR=="MCR").AND.(dE<0)) print *,"MCR while dE <0:",i,fpath
     endif
  enddo
close(2)
close(1)
call mochrej(ofile1,ofile2,NLo1,ngrid)

xlab="(ev)"
call distribution(ofile1,NLo1,2,name1,ngrid,xlab)

title="fractional rejections momenta change"
ax(1)="set yran[0:]"
write(ax(2),'(A5,f14.6)')"beta=",beta
ax(3)="f(x)=1-exp(-beta*x)"
ax(4)= "repl f(x) w l lt 2 lw 5"
call picture(ofile2,2,3,name2,title,xlab,nlx=4,ax=ax)

end subroutine distdE
!ES-----------------------------------------------------------
end module anMC
