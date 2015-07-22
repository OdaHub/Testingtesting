Module transmission 
contains
!BS---------------------------------------------------------
subroutine transout
use var
use inquire
use results
implicit none
character(LEN=XLSTR)::fullfilename
character(LEN=LSTR)::transFILE,transdir
integer::NL

print *,"analyzing Transmission Coefficient"
transdir="trans"
transfile="TRANSMISSION.dat"
call check_datafile(dir,transdir,transFILE,fullfilename,NL,TASK,"TRANSMISSION")
if (NL>1) then
    call  runavkappa(fullfilename,NL,kappa,kappaTST)
    call  kappa_block_analyze(fullfilename,NL,rekap,rekapTST,ncorkap,ncorkapTST)
endif


end subroutine transout
!ES----------------------------------------------------------

!BS------------------------------------------------------------------
subroutine runavkappa(IFILE,NL,av1,av2)
use var
use data_analysis
use gnuplot
implicit none
character(LEN=*), intent(in)::IFILE
integer, intent(in)::NL
double precision, intent(out)::av1,av2
character(LEN=LSTR)::ofile,name
character(LEN=LSTR)::commands(20)
integer::col
commands(1:20)=" "

name="runkappa"
ofile=trim(name)//".tmp"
col=2;call runav(ifile,NL,col,ofile,av1,ngrid=ngrid)
commands(1)="set xlab ""cycle"""
commands(2)="set title ""run av kappa (A/ns)"""
write(commands(3),*) "invtu=",1.d0/timeunit
if (two_point_method) write(commands(3),*) "invtu=",1.d0
commands(4)="pl """//trim(ofile)//"""   u ($1):($2*invtu) w l lw 5"
call make_plot(name,20,commands)

name="runkappaTST"
ofile=trim(name)//".tmp"
col=3;call runav(ifile,NL,col,ofile,av2,ngrid=ngrid)
commands(2)="set title ""run av kappa TST (A/ns)"""
commands(4)="pl """//trim(ofile)//"""   u ($1):($2*invtu) w l lw 5"
write(commands(5),*) "m=",mass_eVns2_A2
write(commands(6),*) "b=",beta
commands(7)="repl 1/(sqrt(2*pi*b*m)) w l lw 5"
call make_plot(name,20,commands)
if (.NOT.two_point_method) then
  av1=av1/timeunit
  av2=av2/timeunit
endif


end subroutine runavkappa
!ES------------------------------------------------------------------


!BS----------------------------------------------------------------------
subroutine kappa_block_analyze(IFILE,NL,re1,re2,ncor1,ncor2)
use var
use data_analysis
use gnuplot
implicit none
integer, intent(in)::NL
character(LEN=*), intent(in)::ifile
double precision, intent(out)::re1,re2,ncor1,ncor2
character(LEN=LSTR)::ofile,name
character(LEN=LSTR)::commands(20)
integer::col

commands=""
name="errkappa";ofile=trim(name)//".tmp"
col=2
call block_error(ifile,NL,col,ofile,maxblocklength,blockskip,relative_e=re1,&
                  Ncorrelation=Ncor1)
commands(1)="set xlab ""block length"""
write(commands(2),*) "set title ""kappa rel.err, Ncor",re1,Ncor1,""""
commands(3)="pl """//trim(ofile)//"""   u 1:3 w l lw 5"
write(commands(4),*) "repl",re1,"w l lw 5"
call make_plot(name,20,commands)

name="errkappaTST";ofile=trim(name)//".tmp"
col=3
call block_error(ifile,NL,col,ofile,maxblocklength,blockskip,relative_e=re2,&
                  Ncorrelation=Ncor2)
write(commands(2),*) "set title ""KIN, rel.err, Ncor",re2,Ncor2,""""
commands(3)="pl """//trim(ofile)//"""   u 1:3 w l lw 5"
write(commands(4),*) "repl",re2,"w l lw 5"
call make_plot(name,4,commands)

end subroutine kappa_block_analyze
!ES--------------------------------------------------------------------


end Module transmission 
