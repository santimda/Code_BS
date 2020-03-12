set encoding iso_8859_1
#set title 'Proton cooling time'
set xrange [9:]
set yrange [5:18]
set ylabel 'log_{10}(t [s])'
set xlabel 'log_{10}(E_p [eV])'
set key box out horizontal top

plot '../tcoolp.dat' u (log10($1)):(log10($2)) t 'ac' w l lw 3,\
      '../tcoolp.dat' u (log10($1)):(log10($3)) t 'pp' w l lw 2,\
      '../tcoolp.dat' u (log10($1)):(log10($4)) t 'conv' w l lw 2,\
      '../tcoolp.dat' u (log10($1)):(log10($5)) t 'dif' w l lw 2,\
      '../tcoolp.dat' u (log10($1)):(log10($9)) t 'adi' w l lw 2,\
      '../tcoolp.dat' u (log10($1)):(log10($10)) t 'ion' w l lw 2,\
      '../tcoolp.dat' u (log10($1)):(log10($6)) t 'cool' w l lw 3

set terminal postscript color solid enhanced dl 2.5 lw 1.7 "Helvetica" 20
set output 'tiemposp.ps'
replot
!./fixbb tiemposp.ps

