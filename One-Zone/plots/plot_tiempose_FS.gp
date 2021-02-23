set encoding iso_8859_1
#set title 'Electron characteristic timescales (forward shock)'
set size square
set xrange [6:13]
set yrange [2:16]
set ylabel 'log_{10}(t [s])' font "Helvetica,18"
set xlabel 'log_{10}(E_e [eV])' font "Helvetica,18"
set xtics font "Helvetica,18"
set ytics 2,2,16 font "Helvetica,18"
set key box out horizontal top font "Helvetica,18"
set label "NAME" at 11,15 font ",18"
set lmargin 0.05
set bmargin 0.05

set terminal postscript color solid enhanced dl 2.5 lw 1.7 "Helvetica" 20
set output 'tiempose_FS.eps'

plot '../tcoole_FS.dat' u (log10($1)):(log10($13)) t 'Total' w l lw 3,\
      '../tcoole_FS.dat' u (log10($1)):(log10($9)) t 'IC,st' w l lw 2,\
      '../tcoole_FS.dat' u (log10($1)):(log10($10)) t 'IC,IR' w l lw 2,\
      '../tcoole_FS.dat' u (log10($1)):(log10($3)) t 'sy' w l lw 2,\
      '../tcoole_FS.dat' u (log10($1)):(log10($4)) t 'br' w l lw 2,\
      '../tcoole_FS.dat' u (log10($1)):(log10($5)) t 'conv' w l lw 2,\
      '../tcoole_FS.dat' u (log10($1)):(log10($6)) t 'ac' w l lw 2,\
      '../tcoole_FS.dat' u (log10($1)):(log10($8)) t 'dif' w l lw 2,\
      '../tcoole_FS.dat' u (log10($1)):(log10($15)) t 'adi' w l lw 2 lc rgb "pink",\
      '../tcoole_FS.dat' u (log10($1)):(log10($7)) t 'cool' w l lw 2

      #'../tcoole_FS.dat' u (log10($1)):(log10($14)) t 'coul' w l lw 2 lc rgb "grey50",\

!./fixbb tiempose_FS.eps
!epstopdf tiempose_FS.eps