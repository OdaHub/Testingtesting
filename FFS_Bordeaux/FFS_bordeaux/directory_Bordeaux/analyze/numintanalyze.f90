Module numintanalyze
implicit none

contains

!BS---------------------------------------------------------
subroutine numintout
use var
use results
use inquire
use gnuplot
implicit none
character(LEN=XLSTR)::fullfilename
character(LEN=LSTR)::NUMINTFILE,NUMINTdir
integer::NL
character(LEN=LSTR)::name
character(LEN=LSTR)::commands(50)

print *,"analyzing numerical integration"

commands=""
numintdir="numint"
numintFILE="NUMINT.dat"

  call check_datafile(dir,numintdir,numintFILE,fullfilename,&
                                    NL,TASK,"NUMERICINTEG")
  
  if (NL>1) then
    name="Probability"
    commands(1)="set xlab ""x"""
    commands(2)="set title ""Probability"""
    commands(3)="pl """//trim(fullfilename)//"""   u 1:2 w l lw 5"
    if (dim>1)  commands(4)="repl """//trim(fullfilename)//""" u 4:5 w l lw 5"
    call make_plot(name,4,commands)

    name="FreeEnergy"
    commands(2)="set title ""Free Energy"""
    commands(3)="pl """//trim(fullfilename)//"""   u 1:3 w l lw 5"
    if (dim>1) commands(4)="repl """//trim(fullfilename)//""" u 4:6 w l lw 5"
    call make_plot(name,4,commands)

    call probinterface(trans_INTFM,trans_INTFR,fullfilename,NL,prob1,prob2)
  endif

end subroutine numintout
!ES----------------------------------------------------------

!BS---------------------------------------------------------------
subroutine probinterface(INTFM,INTFR,filenm,NL,prob1,prob2)
implicit none
double precision, intent(in)::INTFM,INTFR
character(LEN=*), intent(in)::filenm
integer, intent(in)::NL
double precision, intent(out)::prob1,prob2
double precision::dx,x1,x2,x
double precision::probarr(NL)
integer::i,j,im,ir
double precision::epsilon,Fmin,F,xmin,Fmax,xmax,xend,Fend

prob1=0.d0
prob2=0.d0
Fmin=10000000.d0
xmin=10000000.d0
xmax=-1000000.d0
Fmax=-1000000.d0

open(1,file=filenm)
  read(1,*) x1
  read(1,*) x2
close(1)
dx=x2-x1
print *,"dx=",dx

epsilon=dx/10
print *, INTFM,INTFR
im=0;ir=0
open(1,file=filenm)
 do i=1,NL
   read(1,*) x,probarr(i),F
   if (abs(x-INTFM)<epsilon) im=i
   if (abs(x-INTFR)<epsilon) ir=i
   if (F<Fmin) then
     Fmin=F;xmin=x
   endif 
   if ((F>FMax).AND.(x>xmin)) then
     Fmax=F;xmax=x
   endif
 enddo 
close(1)
Fend=F
xend=x

print *,"--------"
print *,"x,Fmin",xmin,Fmin
print *,"x,Fmax",xmax,Fmax
print *,"x,Fend",xend,Fend
print *,"DF1, DF2", Fmax-Fmin, Fmax-Fend
print *,"--------"



if ((im==0).OR.(ir==0)) then
  print *,"non-compatible lattice"
  return 
endif

prob1=probarr(im)*dx
prob2=probarr(ir)*dx
i=im
j=ir
do
 i=i-1
 j=j-1
 if (i>0) prob1=prob1+4*probarr(i)*dx
 if (j>0) prob2=prob2+4*probarr(j)*dx
 i=i-1
 j=j-1
 if (i>0) prob1=prob1+2*probarr(i)*dx
 if (j>0) prob2=prob2+2*probarr(j)*dx
 if (j<0) exit
enddo

prob1=prob1/3
prob2=prob2/3
prob1=probarr(im)/prob1
prob2=probarr(im)/prob2

end subroutine probinterface
!ES---------------------------------------------------------------

end Module numintanalyze
!EM-------------------------------------------------------------
