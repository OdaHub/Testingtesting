Module analyze_EN
contains

!BS-----------------------------------------------------------------
subroutine Evstime(IFILE,NL)
use var
use data_analysis
use gnuplot
implicit none
integer, intent(in)::NL
character(LEN=*), intent(in)::ifile
character(LEN=LSTR)::ofile,name
character(LEN=LSTR)::commands(50)

commands=""
name="Evst";ofile=trim(name)//".tmp"
call sparse(ifile,NL,ofile,ngrid)
commands(1)="set xlab ""t ("//trim(chartu)//")"""
commands(2)="set ylab ""("//trim(chareu)//")"""
commands(3)="set title ""POT, KIN, ET, EHAM"""
write(commands(4),*) "dt=",dt_ut
commands(5)="pl """//trim(ofile)//"""   u ($1*dt):2 w l lw 5"
commands(6)="repl """//trim(ofile)//""" u ($1*dt):3 w l lw 5"
commands(7)="repl """//trim(ofile)//""" u ($1*dt):4 w l lw 5"
commands(8)="repl """//trim(ofile)//""" u ($1*dt):5 w l lw 5"
call make_plot(name,50,commands)

commands=""
name="Tvst"
commands(1)="set xlab ""t ("//trim(chartu)//")"""
commands(2)="set ylab ""("//trim(chartempu)//")"""
commands(3)="set title ""TEMP"""
write(commands(4),*) "dt=",dt_ut
commands(5)="pl """//trim(ofile)//"""   u ($1*dt):6 w l lw 5"
call make_plot(name,50,commands)

end subroutine Evstime
!ES-----------------------------------------------------------------


!BS------------------------------------------------------------------
subroutine runavE(IFILE,NL)
use var
use data_analysis
use gnuplot
implicit none
character(LEN=*), intent(in)::IFILE
integer, intent(in)::NL
character(LEN=LSTR)::ofile1,ofile2,ofile3,name,ofile
character(LEN=LSTR)::commands(50)
integer::col
double precision::av
commands=""
name="runE"
ofile1=trim(name)//"1.tmp"
ofile2=trim(name)//"2.tmp"
ofile3=trim(name)//"3.tmp"
col=2;call runav(ifile,NL,col,ofile1,av,ngrid=ngrid)
col=3;call runav(ifile,NL,col,ofile2,av,ngrid=ngrid)
col=4;call runav(ifile,NL,col,ofile3,av,ngrid=ngrid)
commands(1)="set xlab ""t ("//trim(chartu)//")"""
commands(2)="set ylab ""("//trim(chareu)//")"""
commands(3)="set title ""run av POT, KIN, ETOT"""
write(commands(4),*) "dt=",dt_ut
commands(5)="pl """//trim(ofile1)//"""   u ($1*dt):2 w l lw 5"
commands(6)="repl """//trim(ofile2)//"""   u ($1*dt):2 w l lw 5"
commands(7)="repl """//trim(ofile3)//"""   u ($1*dt):2 w l lw 5"
call make_plot(name,50,commands)

commands=""
name="runT";ofile=trim(name)//".tmp"
col=6;call runav(ifile,NL,col,ofile,av,ngrid=ngrid)
commands(1)="set xlab ""t ("//trim(chartu)//")"""
commands(2)="set ylab ""("//trim(chartempu)//")"""
commands(3)="set title ""run av Temp"""
write(commands(4),*) "dt=",dt_ut
commands(5)="pl """//trim(ofile)//"""   u ($1*dt):2 w l lw 5"
call make_plot(name,50,commands)
end subroutine runavE 
!ES------------------------------------------------------------------

!BS-------------------------------------------------------------------
subroutine Edistr(IFILE,NL)
use var
use data_analysis 
use GAMMA_func
use gnuplot 
implicit none
integer, intent(in)::NL
character(LEN=*), intent(in)::ifile
character(LEN=LSTR)::ofile,name
character(LEN=LSTR)::commands(50)
integer::col
double precision::av,stdev,GAM

commands=""
name="POTdist";ofile=trim(name)//".tmp"
col=2
call distribution(ifile,NL,col,ofile,ngrid,av=av,stdev=stdev)
commands(1)="set xlab ""("//trim(chareu)//")"""
write(commands(2),*) "set title ""dis.POT av,sd:",av,stdev,""""
commands(3)="pl """//trim(ofile)//"""   u 1:2 w l lw 5"
call make_plot(name,50,commands)

commands=""
name="KINdist";ofile=trim(name)//".tmp"
col=3
call distribution(ifile,NL,col,ofile,ngrid,av=av,stdev=stdev)
commands(1)="set xlab ""("//trim(chareu)//")"""
write(commands(2),*) "set title ""dis.KIN av,sd:",av,stdev,""""
commands(3)="pl """//trim(ofile)//"""   u 1:2 w l lw 5"
GAM=GAMMA_N2(Npart*dim)
write(commands(4),*)"GAM=",GAM
write(commands(5),*) "N=",Npart*dim
write(commands(6),*) "kbT=",kbT
commands(7)="PK(x)=( x**(N/2.-1)*exp(-x/kbT) )/( GAM*kbT**(N/2.) ) "
commands(8)="repl PK(x) w l lt 2 lw 5"
call make_plot(name,50,commands)

commands=""
name="Etdist";ofile=trim(name)//".tmp"
col=4
call distribution(ifile,NL,col,ofile,ngrid,av=av,stdev=stdev)
commands(1)="set xlab ""("//trim(chareu)//")"""
write(commands(2),*) "set title ""dis.ET av,sd:",av,stdev,""""
commands(3)="pl """//trim(ofile)//"""   u 1:2 w l lw 5"
call make_plot(name,50,commands)

commands=""
name="Tdist";ofile=trim(name)//".tmp"
col=6
call distribution(ifile,NL,col,ofile,ngrid,av=av,stdev=stdev)
commands(1)="set xlab ""("//trim(chartempu)//")"""
write(commands(2),*) "set title ""dis.Temp av,sd:",av,stdev,""""
commands(3)="pl """//trim(ofile)//"""   u 1:2 w l lw 5"
write(commands(4),*)"GAM=",GAM
write(commands(5),*) "N=",Npart*dim
write(commands(6),*) "Temp=",Temp
commands(7)="PT(x)= (x*N/(2*Temp))**(N/2.) *exp(-x*N/(2*Temp) )/( GAM*x ) "
commands(8)="repl PT(x) w l lt 2 lw 5"
call make_plot(name,50,commands)

end subroutine Edistr
!ES--------------------------------------------------------------------

!BS----------------------------------------------------------------------
subroutine E_block_analyze(IFILE,NL)
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

name="errEP";ofile=trim(name)//".tmp"
col=2
call block_error(ifile,NL,col,ofile,maxblocklength,blockskip,relative_e=rel_e,&
                  Ncorrelation=Ncor)
commands(1)="set xlab ""block length"""
write(commands(2),*) "set title ""POT, rel.err, Ncor",rel_e,Ncor,""""
commands(3)="pl """//trim(ofile)//"""   u 1:3 w l lw 5"
write(commands(4),*) "repl",rel_e,"w l lw 5" 
call make_plot(name,4,commands)

name="errEK";ofile=trim(name)//".tmp"
col=3
call block_error(ifile,NL,col,ofile,maxblocklength,blockskip,relative_e=rel_e,&
                  Ncorrelation=Ncor)
commands(1)="set xlab ""block length"""
write(commands(2),*) "set title ""KIN, rel.err, Ncor",rel_e,Ncor,""""
commands(3)="pl """//trim(ofile)//"""   u 1:3 w l lw 5"
write(commands(4),*) "repl",rel_e,"w l lw 5"
call make_plot(name,4,commands)

name="errET";ofile=trim(name)//".tmp"
col=4
call block_error(ifile,NL,col,ofile,maxblocklength,blockskip,relative_e=rel_e,&
                  Ncorrelation=Ncor)
commands(1)="set xlab ""block length"""
write(commands(2),*) "set title ""ET, rel.err, Ncor",rel_e,Ncor,""""
commands(3)="pl """//trim(ofile)//"""   u 1:3 w l lw 5"
write(commands(4),*) "repl",rel_e,"w l lw 5"
call make_plot(name,4,commands)

name="errTemp";ofile=trim(name)//".tmp"
col=6
call block_error(ifile,NL,col,ofile,maxblocklength,blockskip,relative_e=rel_e,&
                  Ncorrelation=Ncor)
commands(1)="set xlab ""block length"""
write(commands(2),*) "set title ""ETemp, rel.err, Ncor",rel_e,Ncor,""""
commands(3)="pl """//trim(ofile)//"""   u 1:3 w l lw 5"
write(commands(4),*) "repl",rel_e,"w l lw 5"
call make_plot(name,4,commands)


end subroutine E_block_analyze
!ES--------------------------------------------------------------------

end Module analyze_EN
!EM-------------------------------------------------------------
