module funcplot

contains

!BS--------------------------------------------------
subroutine plot_functions
use var
implicit none
 print *,"Plot potential functions"
 select case(POTENTIAL)
   case("PBD") 
     call Plot_potPBD
   case("WCA") 
     call Plot_potWCA
   case("DOUBLEWELL")
     call Plot_DW
   case default
     print *,"No potential plot"
 end select
end subroutine plot_functions
!ES--------------------------------------------------
 
!BS--------------------------------------------------
subroutine Plot_PotPBD
use var
use gnuplot
implicit none
character(LEN=MSTR)::name
character(LEN=LSTR)::commands(13)

  name="Morse"
  commands(1) = "set tit ""Morse(x)"" "
  commands(2) = "set xran [-.2:1.5]"
  commands(3) = "set yran [0:.1]"
  commands(4) = "set xlab ""r(angst)"" "
  commands(5) = "set ylab "" "" "
  write( commands(6),*)  "dAT=",dAT
  write( commands(7),*)  "dGC=",dGC
  write( commands(8),*)  "aAT=",aAT
  write( commands(9),*)  "aGC=",aGC
  commands(10)= "M1(x)=dAT*(exp(-aAT*x) - 1.0)**2"
  commands(11)= "M2(x)=dGC*(exp(-aGC*x) - 1.0)**2"
  commands(12)= "plot M1(x) w l lt 1 lw 5"
  commands(13)= "replot M2(x) w l lt 2 lw 5"
  
  call make_plot(name,13,commands)

end subroutine Plot_PotPBD
!ES---------------------------------------------------

!BS--------------------------------------------------
subroutine Plot_DW
use var
use gnuplot
implicit none
character(LEN=MSTR)::name
character(LEN=LSTR)::commands(13)
  commands(:)=""
  name="DW"
  commands(1) = "set tit ""DoubleWell(x)"" "
  commands(2) = "set xran [-2:2]"
  commands(3) = "set yran [-1:.10]"
  commands(4) = "set xlab ""x"" "
  commands(5) = "set ylab "" "" "
  write( commands(6),*)  "k4=",DOUBLEWELLk4
  write( commands(7),*)  "k2=",DOUBLEWELLk2
  commands(10)= "V(x)=k4*x**4+k2*x**2"
  commands(11)= "plot V(x)"
  call make_plot(name,13,commands)

end subroutine Plot_DW
!ES---------------------------------------------------



!BS--------------------------------------------------
subroutine Plot_PotWCA
use var
use gnuplot
implicit none
character(LEN=MSTR)::name
character(LEN=LSTR)::commands(50)
  commands=""
  name="WCA"
  commands(1) = "set tit ""Vwca(x),Vdw(x)"" "
  commands(2) = "set xran [.7:2]"
  commands(3) = "set yran[0:10]"
  commands(4) = "set xlab ""r("//trim(charlu)//")"""
  commands(5) = "set ylab ""V("//trim(chareu)//")"""
  write( commands(6),*)  "epsilon=",WCA_EPSILON
  write( commands(7),*)  "sigma=",WCA_SIGMA
  write( commands(8),*)  "h=",WCA_H
  write( commands(9),*)  "w=",WCA_W
  write( commands(10),*)  "r0=(2**(1./6.))*sigma"
  commands(11)= "Vwca0(x)=4*epsilon*((x/sigma)**(-12)-(x/sigma)**(-6))+epsilon"
  commands(12)= "Vwca(x)=(x<r0) ? Vwca0(x):0"
  commands(13)= "Vdw(x)=h*(1-(x-r0-w)**2/w**2)**2"
  commands(14)= "plot Vwca(x) w l lt 1 lw 5"
  commands(15)= "replot Vdw(x) w l lt 2 lw 5"
  call make_plot(name,50,commands)

end subroutine Plot_PotWCA
!ES---------------------------------------------------



end module funcPlot
