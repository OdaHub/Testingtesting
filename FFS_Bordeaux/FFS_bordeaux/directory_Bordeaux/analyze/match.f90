module match
implicit none

contains

!BS-----------------------------------------
subroutine match_hist
use var
use results
use gnuplot
implicit none
character(LEN=LSTR)::name,ofile,vertline
character(LEN=LSTR)::commands(100)
integer::i,j,loop
double precision::fac

print *,"Matching histograms"
if (numint==0) return

commands(:)=""
name="TOTPCROSS"
commands(1)="set xlab ""lambda"""
commands(2)="set title ""TOT PCROSS"
commands(3)="set log y"
fac=1

ofile="PCROSS."//ext(1)//".tmp"
vertline="vline."//ext(1)//".tmp"
commands(4)="pl """//trim(ofile)//"""   u 1:2 w l lw 5"
commands(5)="repl """//trim(vertline)//"""   u 1:2 w l lw 5"
loop=3
do i=2,numint
  j=5+(i-2)*loop
  ofile="PCROSS."//ext(i)//".tmp"
  fac=fac*crossprobab(i-1)
  write(commands(j+1),*) "fac"//ext(i)//"=",fac
  commands(j+2)="repl """//trim(ofile)//"""   u 1:($2*fac"//ext(i)//") w l lw 5"
  vertline="vline."//ext(i)//".tmp"
  commands(j+3)="repl """//trim(vertline)//"""   u 1:2 w l lw 5"
enddo
call make_plot(name,100,commands)


end subroutine match_hist
!ES-----------------------------------------

end module match
