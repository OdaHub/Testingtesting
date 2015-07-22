Module analyze_OP
contains

!BS-----------------------------------------------------------------
subroutine OPvstime(IFILE,NL)
use var
use data_analysis
use gnuplot
implicit none
integer, intent(in)::NL
character(LEN=*), intent(in)::ifile
character(LEN=LSTR)::ofile,name
character(LEN=LSTR)::commands(7)

name="OPvst";ofile=trim(name)//".tmp"
call sparse(ifile,NL,ofile,ngrid)
commands(1)="set xlab ""t ("//trim(chartu)//")"""
commands(2)="set title ""OP"""
write(commands(3),*) "dt=",dt_ut
commands(4)="pl """//trim(ofile)//"""   u ($1*dt):2 w l lw 5"
call make_plot(name,4,commands)


end subroutine OPvstime
!ES-----------------------------------------------------------------


!BS------------------------------------------------------------------
subroutine runavOP(IFILE,NL)
use var
use data_analysis
use gnuplot
implicit none
character(LEN=*), intent(in)::IFILE
integer, intent(in)::NL
character(LEN=LSTR)::ofile,name
character(LEN=LSTR)::commands(6)
integer::col
double precision::av

name="runOP"
ofile=trim(name)//".tmp"
col=2;call runav(ifile,NL,col,ofile,av,ngrid=ngrid)
commands(1)="set xlab ""t ("//trim(chartu)//")"""
commands(2)="set title ""run av OP"""
write(commands(3),*) "dt=",dt_ut
commands(4)="pl """//trim(ofile)//"""   u ($1*dt):2 w l lw 5"
call make_plot(name,4,commands)

end subroutine runavOP
!ES------------------------------------------------------------------

!BS-------------------------------------------------------------------
subroutine OPdistr(IFILE,NL)
use var
use data_analysis 
use gnuplot 
implicit none
integer, intent(in)::NL
character(LEN=*), intent(in)::ifile
character(LEN=LSTR)::ofile,name
character(LEN=LSTR)::commands(9)
integer::col,NLCOM
double precision::av,stdev

name="OPdist";ofile=trim(name)//".tmp"
col=2
NLCOM=3
call distribution(ifile,NL,col,ofile,ngrid,av=av,stdev=stdev)
commands(1)="set xlab ""(A)"""
write(commands(2),*) "set title ""dis.OP av,sd:",av,stdev,""""
commands(3)="pl """//trim(ofile)//"""   u 1:2 w l lw 5"
if (POTENTIAL=="HARMOSC") then 
  NLCOM=7 
  write(commands(4),*) "kharm=",kharm
  write(commands(5),*) "kbT=",kbT
  commands(6)="P(x)=exp(-0.5*kharm*x**2/kbT)/sqrt(2*pi*kbT/kharm)"
  commands(7)="repl P(x) w l lw 5"
endif
call make_plot(name,NLCOM,commands)

end subroutine OPdistr
!ES--------------------------------------------------------------------

!BS----------------------------------------------------------------------
subroutine OP_block_analyze(IFILE,NL)
use var
use data_analysis
use gnuplot
implicit none
integer, intent(in)::NL
character(LEN=*), intent(in)::ifile
character(LEN=LSTR)::ofile,name
character(LEN=LSTR)::commands(4)
integer::col
double precision::rel_e,Ncor

name="errOP";ofile=trim(name)//".tmp"
col=2
call block_error(ifile,NL,col,ofile,maxblocklength,blockskip,relative_e=rel_e,&
                  Ncorrelation=Ncor)
commands(1)="set xlab ""block length"""
write(commands(2),*) "set title ""OP, rel.err, Ncor",rel_e,Ncor,""""
commands(3)="pl """//trim(ofile)//"""   u 1:3 w l lw 5"
write(commands(4),*) "repl",rel_e,"w l lw 5" 
call make_plot(name,4,commands)

end subroutine OP_block_analyze
!ES--------------------------------------------------------------------

!BS---------------------------------------------------------------
subroutine mean_square_disp_OP(IFILE,NL)
use var
use data_analysis
use gnuplot
implicit none
integer, intent(in)::NL
character(LEN=*), intent(in)::ifile
character(LEN=LSTR)::ofile,name
character(LEN=LSTR)::commands(10)
integer::col

name="msdOP";ofile=trim(name)//".tmp"
col=2
call msd(ifile,NL,col,ofile,msdlength)
commands(1)="set xlab ""t ("//trim(chartu)//")"""
write(commands(2),*) "set title ""MSD"""
write(commands(3),*) "dt=",dt_ut
commands(4)="pl """//trim(ofile)//"""   u ($1*dt):2 w l lw 5"
if ((POTENTIAL=="HARMOSC").AND.(kharm==0.d0).AND.(DYNAMICS=="LANGEVIN")) then
  write(commands(5),*) "D=",kbT/(mass_eVns2_A2*gamma_invut)
  write(commands(6),*) "dim=",1
  commands(7)="repl 2*dim*D*x lw 5"
  call make_plot(name,7,commands) 
else
  call make_plot(name,4,commands)
endif

end subroutine mean_square_disp_OP
!ES---------------------------------------------------------------

end Module analyze_OP
!EM-------------------------------------------------------------
