set encoding iso_8859_1
set title ''
set key off
set ylabel 'log_{10}({E_p}^2 N_p(E_p) [erg])' font "Helvetica,18"
set xlabel 'log_{10}(E_p [eV])' font "Helvetica,18"
set xrange [9:14]
set yrange [29:36]
set ytics 1; set mytics 2 
set xtics 1; set mxtics 2 
set label "NAME" at graph 0.05,0.95 font "Helvetica,18"
set size square
set bmargin 0.2

eV = 1.602176e-12

set terminal postscript color solid enhanced dl 2.5 lw 1.7 "Helvetica" 20
set output 'distp_RS.eps

plot '../dist_p_RS.dat' u (log10($1)):(log10((($1)*eV)**2*($2))) w l lw 3 

!./fixbb distp_RS.eps
!epstopdf distp_RS.eps
