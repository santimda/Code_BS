set encoding iso_8859_1
set title ''
#set key off
set ylabel 'log_{10}({E_e}^2 N_e(E_e) [erg])'
set xlabel 'log_{10}(E_e [eV])'
set xrange [6:14]
set yrange [29:38]
set ytics 1; set mytics 2

eV = 1.602176e-12

plot '../dist_e_FS.dat' u (log10($1)):(log10((($1)*eV)**2.21*($2))) w l lw 3 t '',\
     '../inj_e_FS.dat' u (log10($1)):(log10((($1)*eV)**2.21*($2))) w l lw 2 dt 2 t 'E^2 Q(E)*t_{esc}'

set terminal postscript color solid enhanced dl 2.5 lw 2.7 "Helvetica" 20
set output 'dist_e_FS.ps'
replot
!./fixbb dist_e_FS.ps
