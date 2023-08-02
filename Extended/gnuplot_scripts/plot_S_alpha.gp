set encoding iso_8859_1
#set xtics 100
#set ytics 0.1; set mytics 2 
#set ylabel 'S_{H{/symbol a}} [erg s^{-1} cm^{-3} arcsec^{-2}]'
set ylabel 'S_{H{/symbol a}} [R]'
set xlabel '{/Symbol q} [°]'
#set yrange [0:5]
#set label "NAME" at graph 0.05,0.95 font "Helvetica,18"

set terminal postscript color solid enhanced dl 6.5 lw 2.7 "Helvetica" 20
set output 'S_H_alpha.eps'

plot '../S_alpha_tot_FS.dat' u (180.*($1)/pi):(($2)/(5.66e-18)) t ''

!./fixbb S_H_alpha.eps
!epstopdf S_H_alpha.eps
