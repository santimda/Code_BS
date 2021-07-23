set encoding iso_8859_1
set title 'E_{e,max} in RS'
#set xtics 100
set ytics 1; set mytics 10 
set key off
set ylabel 'log_{10}(E_{e,max} [eV])'
set xlabel 'position l'

plot '../EmaxeRS.dat' u (($1)):(log10(($2))) w l lw 3 

set terminal postscript color solid enhanced dl 6.5 lw 2.7 "Helvetica" 20
set output 'EmaxeRS.eps'
replot
!./fixbb EmaxeRS.eps

