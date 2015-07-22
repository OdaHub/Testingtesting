S1=0
plot "ENERGIES.dat" u 1:($2+S1) w l
replot "ENERGIES.dat" u 1:3 w l
replot "ENERGIES.dat" u 1:($4+S1) w l
