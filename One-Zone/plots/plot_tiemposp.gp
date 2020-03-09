set encoding iso_8859_1
set title 'Proton cooling time'
set xrange [9:14]
set yrange [:18]
set ylabel 'log_{10}(t [s])'
set xlabel 'log_{10}(E_p [eV])'

plot '../tcoolp.dat' u (log10($1)):(log10($2)) t 't_{ac}' w l lw 3,\
      '../tcoolp.dat' u (log10($1)):(log10($3)) t 't_{pp}' w l lw 2,\
      '../tcoolp.dat' u (log10($1)):(log10($4)) t 't_{ad}' w l lw 2,\
      '../tcoolp.dat' u (log10($1)):(log10($5)) t 't_{dif}' w l lw 2,\
      '../tcoolp.dat' u (log10($1)):(log10($6)) t 't_{cool}' w l lw 3

set terminal postscript color solid enhanced dl 6.5 lw 1.7 "Helvetica" 20
set output 'tiemposp.ps'
replot
!./fixbb tiemposp.ps

