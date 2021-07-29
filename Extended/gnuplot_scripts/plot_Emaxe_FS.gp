set encoding iso_8859_1
set title 'E_{e,max} in FS'
#set xtics 100
set ytics 1; set mytics 10 
set key off
set ylabel 'log_{10}(E_{e,max} [eV])'
set xlabel 'position l'
set label "NAME" at graph 0.05,0.95 font "Helvetica,18"

plot '../EmaxeFS.dat' u (($1)):(log10(($2))) w l lw 3 

set terminal postscript color solid enhanced dl 6.5 lw 2.7 "Helvetica" 20
set output 'EmaxeFS.eps'
replot
!./fixbb EmaxeFS.eps
!epstopdf EmaxeFS.eps

