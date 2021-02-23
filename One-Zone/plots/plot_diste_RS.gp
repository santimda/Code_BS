set encoding iso_8859_1
set title ''
#set key off
set ylabel 'log_{10}({E_e}^2 N_e(E_e) [erg])' font "Helvetica,18"
set xlabel 'log_{10}(E_e [eV])' font "Helvetica,18"
set xrange [6:14]
#set yrange [30:39]
set ytics 1 font "Helvetica,18"; set mytics 2 
set xtics 1 font "Helvetica,18"; set mytics 2 
set bmargin 0.2
set label "NAME" at graph 0.05,0.95 font "Helvetica,18"

set terminal postscript color solid enhanced dl 2.5 lw 1.7 "Helvetica" 20
set output 'diste_RS.eps'

eV = 1.602176e-12

plot '../dist_e_RS.dat' u (log10($1)):(log10((($1)*eV)**2*($2))) w l lw 3 t '',\
     '../inj_e_RS.dat' u (log10($1)):(log10((($1)*eV)**2*($2))) w l lw 2 dt 2 t 'E^2 Q(E)*t_{esc}'

!./fixbb diste_RS.eps
!epstopdf diste_RS.eps
