#set size ratio 1.
set size ratio 3.5/4.5 # yrange/xrange

set xrange [-1.1:]
set yrange [0.0:]
set xlabel 'x/R_0'
set ylabel 'y/R_0'
#set xtics 0.5
#set mxtic 2
set ytics 0.5; set mytics 2
set key box opaque width 0.95 height 1.0 top left

R_st = 15.*6.96e10/3.0856e18

set object 1 circle arc [0:180] at 0,0 size R_st*1.e5 front fc rgb "blue" fs solid 1.0
set label 1 "Star" at 1./25.,1./8. front textcolor rgb "blue" 

plot '../CD_xy.dat' u 1:2 w l lw 3 lc 2 t 'Christie+2016',\
	'../CD_xy_aprox.dat' u 1:2 w l lw 3 lc 1 dt 2 t 'Wilkin1996',\
	'../V/kk.dat' u (-$1):($2) w l lw 2 lc 3 dt 4 t 'Valenti'


set terminal postscript color solid enhanced dl 1.5 lw 2.0 "Helvetica" 20
set output 'CD_xy.ps'
replot
!./fixbb CD_xy.ps