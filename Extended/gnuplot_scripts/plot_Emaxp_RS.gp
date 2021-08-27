set encoding iso_8859_1
set title 'E_{p,max} in RS'
set size ratio 0.9
#set xtics 100
set ytics 1; set mytics 10 
set key off
set ylabel 'log_{10}(E_{p,max} [eV])'
set xlabel 'position l'
set label "NAME" at graph 0.05,0.95 font "Helvetica,18"

set terminal postscript color solid enhanced dl 6.5 lw 2.7 "Helvetica" 20
set output 'EmaxpRS.eps'

plot '../EmaxpRS.dat' u (($1)):(log10(($2))) w l lw 3 

!./fixbb EmaxpRS.eps
!epstopdf EmaxpRS.eps

