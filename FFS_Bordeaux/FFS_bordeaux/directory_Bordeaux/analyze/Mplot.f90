module Mplot

contains

!BS--------------------------------------------------
subroutine plotgraphs
implicit none

 call PlotPot
 call Scurve 
end subroutine plotgraphs
!ES--------------------------------------------------
 
!BS--------------------------------------------------
subroutine PlotPot
use var
use var_analyze
use gnupf
implicit none
character*100::name
character*100::ax(13)

  name="Morse"
  ax(1) = "set tit ""Morse(x)"" "
  ax(2) = "set xran [-.2:1.5]"
  ax(3) = "set yran [0:.1]"
  ax(4) = "set xlab ""r(angst)"" "
  ax(5) = "set ylab "" "" "
  write( ax(6),*)  "dAT=",dAT
  write( ax(7),*)  "dGC=",dGC
  write( ax(8),*)  "aAT=",aAT
  write( ax(9),*)  "aGC=",aGC
  ax(10)= "M1(x)=dAT*(exp(-aAT*x) - 1.0)**2"
  ax(11)= "M2(x)=dGC*(exp(-aGC*x) - 1.0)**2"
  ax(12)= "plot M1(x) w l lt 1 lw 5"
  ax(13)= "replot M2(x) w l lt 2 lw 5"
  
  call make_plot(name,13,ax)

end subroutine
!ES---------------------------------------------------

!BS-------------------------------------------------------------------
subroutine Scurve
use var
use var_analyze
use gnupf
implicit none
character*100::name
character*100::ax(7)

name="Scurve"

ax(1)= "set tit ""S-curve"" "
ax(2)= "set xlab ""r"" "
write( ax(3),*)  "k_open=",k_open
write( ax(4),*)  "y_open=",y_open
ax(5)= "S(x)=exp(k_open*(x-y_open))/(  exp(k_open*(x-y_open))+1.)"
ax(6)= "set xran [-.2:1.5]"
ax(7)= "plot S(x) w l lt 1 lw 5"
call make_plot(name,7,ax)

end subroutine Scurve
!ES-------------------------------------------------------------------

end module MPlot
