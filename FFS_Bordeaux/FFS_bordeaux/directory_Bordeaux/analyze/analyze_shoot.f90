Module analyze_shoot
contains


!BS-------------------------------------------------------------------
subroutine shootout
use var
use inquire
use results
implicit none
character(LEN=XLSTR)::fullfilename
character(LEN=LSTR)::PATHFILE
integer::i

PATHFILE="PATH.dat"
fullfilename=trim(dir)//"/"//trim(PATHFILE)

if ((TASK=="TIS").AND.(nonempty(fullfilename))) then
  call shootstat(fullfilename,ext(1),tiscycles(1),ngrid)
else
  print *,"analyzing shoots TIS directories. Number of ensmenbles ",numint
  do i=1,numint
    fullfilename=trim(dir)//"/"//ext(i)//"/"//trim(PATHFILE)
    call shootstat(fullfilename,ext(i),tiscycles(i),ngrid)
  enddo
  if (TASK=="PPS") then
    fullfilename=trim(dir)//"/000/"//trim(PATHFILE)
    call shootstat(fullfilename,"000",Nl0,ngrid)
  endif

endif

end subroutine shootout 
!ES--------------------------------------------------------------------

!BS----------------------------------------------------------------------
subroutine shootstat(ifile,ext,nl,ngrid)
use stringlengths
use data_analysis
use gnuplot
implicit none
character(LEN=*), intent(in)::ifile,ext
integer, intent(in)::nl,ngrid
integer::ish,isha,ishr,ishb
character(len=LSTR)::fallshoot,faccshoot,frejshoot,fbwishoot
double precision::av,stdev
character(LEN=LSTR)::ofile1,ofile2,ofile3,ofile4,name
character(LEN=LSTR)::commands(20)
integer::col

commands=""
call makeshootfiles(ifile,ext,nl,fallshoot,faccshoot,frejshoot,fbwishoot,ish,isha,ishr,ishb)
name="SHOOT"//ext;col=2
ofile1=trim(name)//"_ALL.tmp"
ofile2=trim(name)//"_ACC.tmp"
ofile3=trim(name)//"_REJ.tmp"
ofile4=trim(name)//"_BWI.tmp"
call distribution(fallshoot,ish,col,ofile1,ngrid,av=av,stdev=stdev)
call distribution(faccshoot,isha,col,ofile2,ngrid,av=av,stdev=stdev)
call distribution(frejshoot,ishr,col,ofile3,ngrid,av=av,stdev=stdev)
call distribution(fbwishoot,ishb,col,ofile4,ngrid,av=av,stdev=stdev)
commands(1)="set xlab ""lambda"""
commands(2)="set title ""SHOOTS"""
commands(3)="pl """//trim(ofile1)//"""   u 1:2 w l lw 5"
commands(4)="repl """//trim(ofile2)//"""   u 1:2 w l lw 5"
commands(5)="repl """//trim(ofile3)//"""   u 1:2 w l lw 5"
commands(6)="repl """//trim(ofile4)//"""   u 1:2 w l lw 5"
call make_plot(name,20,commands)

name="SHOOT2_"//ext
write(commands(1),'(a30,2e12.3,a1)')"set xlab ""%acc,%bwi:",1.d0*isha/ish,1.d0*ishb/ish,""""
write(commands(2),'(a30,e12.3,a1)')"set title ""SHOOTS, %sh:",1.d0*ish/nl,""""
commands(3)="pl """//trim(ofile1)//"""   u 1:2 w l lw 5"
write(commands(4),'(a50,e12.4,a10)') "repl """//trim(ofile2)//"""   u 1:($2*",1.d0*isha/ish,") w l lw 5"
write(commands(5),'(a50,e12.4,a10)') "repl """//trim(ofile3)//"""   u 1:($2*",1.d0*ishr/ish,") w l lw 5"
write(commands(6),'(a50,e12.4,a10)') "repl """//trim(ofile4)//"""   u 1:($2*",1.d0*ishb/ish,") w l lw 5"
call make_plot(name,20,commands)


end subroutine shootstat
!ES---------------------------------------------------------------------


!BS----------------------------------------------------------------------
subroutine makeshootfiles(ifile,ext,nl,fallshoot,faccshoot,frejshoot,fbwishoot,ish,isha,ishr,ishb)
use stringlengths
implicit none
character(LEN=*), intent(in)::ifile,ext
integer, intent(in)::nl
character(len=LSTR),intent(out) ::fallshoot,faccshoot,frejshoot,fbwishoot
integer, intent(out)::ish,isha,ishr,ishb
character::dum(20)
character*3::ACCREJ
character*2::mcmove
double precision::lshoot
integer::i


fallshoot="shootall."//trim(EXT)//".tmp"
faccshoot="shootacc."//trim(EXT)//".tmp"
frejshoot="shootrej."//trim(EXT)//".tmp"
fbwishoot="shootbiw."//trim(EXT)//".tmp"
open(1,file=ifile)
open(2,file=fallshoot)
open(3,file=faccshoot)
open(4,file=frejshoot)
open(5,file=fbwishoot)
ish=0;isha=0;ishr=0;ishb=0
do i=1,nl
  read(1,*) dum(1:7),ACCREJ,mcmove,dum(1:4),lshoot
  if (mcmove=="sh") then
    ish=ish+1
    write(2,*) ish,lshoot
    if (ACCREJ=="ACC") then
      isha=isha+1
      write(3,*) isha,lshoot
    else
      ishr=ishr+1
      write(4,*) ishr,lshoot
    endif  
    if (ACCREJ=="BWI") then
      ishb=ishb+1
      write(5,*) ishb,lshoot
    endif
  endif

enddo
close(5)
close(4)
close(3)
close(2)
close(1)
end subroutine makeshootfiles 
!ES---------------------------------------------------------------------

end Module analyze_shoot
!EM-------------------------------------------------------------
