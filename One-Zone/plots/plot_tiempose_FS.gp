set encoding iso_8859_1
#set title 'Electron characteristic timescales (forward shock)'
set size square
set xrange [6:12.2]
set yrange [3:19]
set ylabel 'log_{10}(t [s])'
set xlabel 'log_{10}(E_e [eV])'
set key box out horizontal top

set terminal postscript color solid enhanced dl 2.5 lw 1.7 "Helvetica" 20
set output 'tiempose_FS.eps'

plot '../tcoole_FS.dat' u (log10($1)):(log10($2)) t 'IC' w l lw 2,\
      '../tcoole_FS.dat' u (log10($1)):(log10($9)) t 'IC,st' w l dt 2,\
      '../tcoole_FS.dat' u (log10($1)):(log10($10)) t 'IC,IR' w l dt 3 lw 2,\
      '../tcoole_FS.dat' u (log10($1)):(log10($3)) t 'sy' w l lw 2,\
      '../tcoole_FS.dat' u (log10($1)):(log10($4)) t 'br' w l lw 2,\
      '../tcoole_FS.dat' u (log10($1)):(log10($5)) t 'conv' w l lw 2,\
      '../tcoole_FS.dat' u (log10($1)):(log10($6)) t 'ac' w l lw 2,\
      '../tcoole_FS.dat' u (log10($1)):(log10($8)) t 'dif' w l lw 2,\
      '../tcoole_FS.dat' u (log10($1)):(log10($14)) t 'coul' w l lw 2 lc rgb "grey50",\
      '../tcoole_FS.dat' u (log10($1)):(log10($15)) t 'adi' w l lw 2 lc rgb "pink",\
      '../tcoole_FS.dat' u (log10($1)):(log10($7)) t 'cool' w l lw 3

!./fixbb tiempose_FS.eps
!epstopdf tiempose_FS.eps