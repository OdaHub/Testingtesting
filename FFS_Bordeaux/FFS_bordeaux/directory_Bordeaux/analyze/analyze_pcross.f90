Module analyze_pcross
contains

!BS-----------------------------------------------------------------
subroutine CROSSPROB(IFILE,NL,extension,im,ir,idetect)
use var
use data_analysis
use gnuplot
implicit none
integer, intent(in)::NL
character(LEN=*), intent(in)::ifile,extension
double precision, intent(in)::im,ir,idetect
character(LEN=LSTR)::ofile,name,vertline
character(LEN=LSTR)::commands(20)
integer::col

commands=""
col=11
name="PCROSS."//trim(extension);ofile=trim(name)//".tmp"
call cross(ifile,NL,col,ofile,ngrid,im,ir)
vertline="vline."//trim(extension)//".tmp"
call vline(idetect,vertline)
commands(1)="set xlab ""lambda"""
commands(2)="set title ""PCROSS"//trim(extension)//""""
commands(3)="pl """//trim(ofile)//"""   u 1:2 w l lw 5"
commands(4)="repl """//trim(vertline)//"""   u 1:2 w l lw 5"
if (POTENTIAL=="HARMOSC") then
  write(commands(5),*)"k=",kharm
  write(commands(6),*)"x0=",im
  commands(7)="dE(x)=0.5*k*(x**2-x0**2)"
  write(commands(8),*)"b=",beta
  commands(9)="f(x)=exp(-b*dE(x))"
  commands(10)="repl f(x) w l lw 5" 
endif
call make_plot(name,10,commands)

end subroutine CROSSPROB 
!ES--------------------------------------------------------------------

!BS------------------------------------------------------------------
subroutine runavPCROSS(IFILE,NL,extension,av)
use var
use data_analysis
use gnuplot
implicit none
character(LEN=*), intent(in)::IFILE,extension
integer, intent(in)::NL
double precision, intent(out)::av
character(LEN=LSTR)::ofile,name
character(LEN=LSTR)::commands(20)
integer::col

name="runPCROSS."//trim(extension)
ofile=trim(name)//".tmp"
col=2
call runav(ifile,NL,col,ofile,av,ngrid=ngrid)
commands(1)="set xlab ""cycles"""
commands(2)="set title ""run av PCROSS"//trim(extension)//""""
commands(3)="pl """//trim(ofile)//"""   u 1:2 w l lw 5"
call make_plot(name,3,commands)

end subroutine runavPCROSS
!ES---------------------------------------------------------------

!BS----------------------------------------------------------------------
subroutine PCROSS_error(IFILE,NL,extension,rel_e,Ncor,mbl,blsk)
use var
use data_analysis
use gnuplot
implicit none
integer, intent(in)::NL,mbl,blsk
character(LEN=*), intent(in)::ifile,extension
double precision, intent(out)::rel_e,Ncor
character(LEN=LSTR)::ofile,name
character(LEN=LSTR)::commands(20)
integer::col

name="errPCROSS."//trim(extension);ofile=trim(name)//".tmp"
col=2
call block_error(ifile,NL,col,ofile,mbl,blsk,relative_e=rel_e,&
                  Ncorrelation=Ncor)
commands(1)="set xlab ""block length"""
write(commands(2),*) "set title ""PCROSS"//trim(extension)//"rel.err, Ncor",rel_e,Ncor,""""
commands(3)="pl """//trim(ofile)//"""   u 1:3 w l lw 5"
write(commands(4),*) "repl",rel_e,"w l lw 5"
call make_plot(name,4,commands)

end subroutine PCROSS_error 
!ES--------------------------------------------------------------------

!BS-------------------------------------------------------------------
subroutine pathdistr(IFILE1,IFILE2,NL,extension,av1,av2)
use var
use data_analysis
use gnuplot
implicit none
integer, intent(in)::NL
character(LEN=*), intent(in)::ifile1,ifile2,extension
double precision, intent(out)::av1,av2
character(LEN=LSTR)::ofile1,ofile2,name
character(LEN=LSTR)::commands(20)
integer::col1,col2
double precision::stdev1,stdev2

name="L."//trim(extension)
ofile1=trim(name)//"_1.tmp"
ofile2=trim(name)//"_2.tmp"
col1=7;col2=2
call distribution(ifile1,NL,col1,ofile1,ngrid,av=av1,stdev=stdev1)
call distribution(ifile2,NL,col2,ofile2,ngrid,av=av2,stdev=stdev2)
commands(1)="set xlab ""(MD steps)"""
write(commands(2),'(a28,4f8.2,a1)') "set title ""av1,av2,sd1,sd2:",av1,av2,stdev1,stdev2,""""
commands(3)="pl """//trim(ofile1)//"""   u 1:2 w l lw 5"
commands(4)="repl """//trim(ofile2)//"""   u 1:2 w l lw 5"
call make_plot(name,4,commands)

end subroutine pathdistr
!ES--------------------------------------------------------------------


end Module analyze_pcross
!EM-------------------------------------------------------------
