set encoding iso_8859_1
#set title 'Electron characteristic timescales'
set xrange [6:14]
set yrange [3:19]
set ylabel 'log_{10}(t [s])'
set xlabel 'log_{10}(E_e [eV])'
set key box out horizontal top

plot '../tcoole.dat' u (log10($1)):(log10($2)) t 'IC' w l lw 2,\
      '../tcoole.dat' u (log10($1)):(log10($9)) t 'IC,st' w l dt 2,\
      '../tcoole.dat' u (log10($1)):(log10($10)) t 'IC,IR' w l dt 3 lw 2,\
      '../tcoole.dat' u (log10($1)):(log10($3)) t 'sy' w l lw 2,\
      '../tcoole.dat' u (log10($1)):(log10($4)) t 'br' w l lw 2,\
      '../tcoole.dat' u (log10($1)):(log10($5)) t 'conv' w l lw 2,\
      '../tcoole.dat' u (log10($1)):(log10($6)) t 'ac' w l lw 2,\
      '../tcoole.dat' u (log10($1)):(log10($8)) t 'dif' w l lw 2,\
      '../tcoole.dat' u (log10($1)):(log10($14)) t 'coul' w l lw 2 lc rgb "grey50",\
      '../tcoole.dat' u (log10($1)):(log10($15)) t 'adi' w l lw 2 lc rgb "pink",\
      '../tcoole.dat' u (log10($1)):(log10($7)) t 'cool' w l lw 3


set terminal postscript color solid enhanced dl 2.5 lw 1.7 "Helvetica" 20
set output 'tiempose.ps'
replot
!./fixbb tiempose.ps

