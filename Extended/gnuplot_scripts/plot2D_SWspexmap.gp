
set encoding iso_8859_1
set title 'Spex with 4.86/1.42 GHz data'
set title ''
set xlabel 'x [arcmin]'
set ylabel 'y [arcmin]'

set cblabel 'spex'
#set yrange [-6.0:6.0]
#set xrange [-10:10]
#set cbrange [-1.5:-0.2]

set pm3d interpolate 0,0
set view map
#set palette red
#set palette gray negative
unset key
#set size ratio 12./20.
set size ratio 6./4.

#set label 1 "" at 0.0,0.0 front point pointtype 7 lc 3 lw 2 ps 1
set label 2 "Star" at 1.0,0.0 font "Helvetica,21" front textcolor lt 3

#set object 5 ellipse center graph 0.91,0.91 size 0.5,0.5 angle 0. fc rgb "grey" fs solid 1.0 border -1 front
set label 2 "{/ZapfDingbats \110}" at -0.36,0.0 front textcolor rgb '#008B8B' font "Helvetica,28"
#set arrow 4 from -0.01,0.1 to -0.096,1.597 head filled front linecolor rgb '#6495ED' linewidth 2.0

splot '../SW_spex.dat' u ($1/60.):(-$2/60.):($3) w pm3d,\
#      '../R0_proj.dat' u ($1/60.):($2/60.):1 ps 3 lw 2 pt 1 lc rgb "gold"

set terminal postscript colour enhanced dl 3.5 lw 2. "Helvetica" 21
set output 'spexmap.ps'
replot
!./fixbb spexmap.ps
