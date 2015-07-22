module analyze_flux 

contains

!BS--------------------------------------------
subroutine runavFlux(ifile1,ifile2,ifile3,nl,av1,av2,av3) 
use var
use data_analysis
use gnuplot
implicit none
character(LEN=*), intent(in)::IFILE1,IFILE2,IFILE3
integer, intent(in)::NL
double precision, intent(out)::av1,av2,av3
character(LEN=LSTR)::ofile1,ofile2,ofile3,name1,name2,name3
character(LEN=LSTR)::commands(50)
integer::col

col=3
!make 3 running average files
name1="runFLUX.1"
ofile1=trim(name1)//".tmp"
call runav(ifile1,NL,col,ofile1,av1,ngrid=ngrid)
name2="runFLUX.2"
ofile2=trim(name2)//".tmp"
call runav(ifile2,NL,col,ofile2,av2,ngrid=ngrid)
name3="runFLUX.3"
ofile3=trim(name3)//".tmp"
call runav(ifile3,NL,col,ofile3,av3,ngrid=ngrid)

commands=""
commands(1)="set xlab ""t ("//trim(chartu)//")"""
commands(2)="set ylab ""1/("//trim(chartu)//")"""
commands(3)="set title ""run av Flux"""
write(commands(4),*) "dt=",dt_ut
write(commands(5),*) "invtu=",1.d0/timeunit
commands(6)="pl """//trim(ofile1)//"""   u ($1*dt):($2*invtu) w l lw 5"
!Compare with analytic solution in the case of the HO
if (POTENTIAL=="HARMOSC") then
  write(commands(7),*) "k=",KHARM
  write(commands(8),*) "m=",mass_eVns2_A2
  write(commands(9),*) "b=",beta
  write(commands(10),*) "q1=",MD_INTFL
  write(commands(11),*) "q2=",MD_INTFR
  commands(12)=&
  "fluxTST=sqrt(k/m)*( exp(-b*.5*k*q1**2)/(pi*erfc( -sqrt(b*k/2)*q2  ) ) )"
  !commands(11)="fluxTST=sqrt(k/m)*( exp(-b*.5*k*q**2)/(2*pi) )"
  commands(13)="repl fluxTST w l lw 5" 
endif
!int_{-inf}^{inf} exp(-b*.5*k*x^2)=sqrt[(2 pi)/(b*k)]=Zx
!kappaTST=1/sqrt[2 pi b m]
!kappaTST/Zx=sqrt(k/m)/(2*pi)
!int_{-inf}^{y} exp(-b*.5*k*x^2)=sqrt[pi/(2*b*k)]*Erfc[-sqrt[(b*k)/2]*y]
!Erfc[-Infty]=2
call make_plot(name1,50,commands)

!only change line 6 and 10
commands(6)="pl """//trim(ofile2)//"""   u ($1*dt):($2*invtu) w l lw 5"
write(commands(10),*) "q1=",MD_INTFM
call make_plot(name2,50,commands)

commands(6)="pl """//trim(ofile3)//"""   u ($1*dt):($2*invtu) w l lw 5"
write(commands(10),*) "q1=",MD_INTFR
call make_plot(name3,50,commands)

av1=av1/timeunit;av2=av2/timeunit;av3=av3/timeunit

end subroutine runavFlux 
!ES---------------------------------------------

!BS-------------------------------------------
subroutine flux_block_analyze(ifile1,ifile2,ifile3,nl, &
              rel_e1,rel_e2,rel_e3) !,Ncor1,Ncor2,Ncor3)
use var
use data_analysis
use gnuplot
implicit none
character(LEN=*), intent(in)::IFILE1,IFILE2,IFILE3
integer, intent(in)::NL
double precision, intent(out)::rel_e1,rel_e2,rel_e3 !,Ncor1,Ncor2,Ncor3
character(LEN=LSTR)::ofile1,ofile2,ofile3,name1,name2,name3
character(LEN=LSTR)::commands(50)
integer::col

name1="errflux.1";ofile1=trim(name1)//".tmp"
col=3
call block_error(ifile1,NL,col,ofile1,maxblocklength,blockskip,relative_e=rel_e1) !,&
                 ! Ncorrelation=Ncor1)
commands=""
commands(1)="set xlab ""block length"""
write(commands(2),*) "set title ""flux, rel.err",rel_e1,""""
commands(3)="pl """//trim(ofile1)//"""   u 1:3 w l lw 5"
write(commands(4),*) "repl",rel_e1,"w l lw 5"
call make_plot(name1,50,commands)

name2="errflux.2";ofile2=trim(name2)//".tmp"
col=3
call block_error(ifile2,NL,col,ofile2,maxblocklength,blockskip,relative_e=rel_e2) !,&
                !  Ncorrelation=Ncor2)
commands(1)="set xlab ""block length"""
write(commands(2),*) "set title ""flux, rel.err",rel_e2,""""
commands(3)="pl """//trim(ofile2)//"""   u 1:3 w l lw 5"
write(commands(4),*) "repl",rel_e2,"w l lw 5"
call make_plot(name2,50,commands)

name3="errflux.3";ofile3=trim(name3)//".tmp"
col=3
call block_error(ifile3,NL,col,ofile3,maxblocklength,blockskip,relative_e=rel_e3) !,&
                !  Ncorrelation=Ncor3)
commands(1)="set xlab ""block length"""
write(commands(2),*) "set title ""flux, rel.err",rel_e3,""""
commands(3)="pl """//trim(ofile3)//"""   u 1:3 w l lw 5"
write(commands(4),*) "repl",rel_e3,"w l lw 5"
call make_plot(name3,50,commands)


end subroutine flux_block_analyze
!ES--------------------------------------------


end module analyze_flux 
