set encoding iso_8859_1
ang_deg = 20.
ang = ang_deg*pi/180.  # rotation angle to match observed map (Ben2020)
set title 'Spex with 10./6. GHz data'
set title ''
set xlabel 'x [arcmin]'
set ylabel 'y [arcmin]'

set cblabel 'spex'
set yrange [-2:2]
set xrange [-2:2]
#set palette defined (-1 "purple", 0 "white", 1 "red")
set cbrange [-0.55:0.2]

set xtics 2; set mxtics 2; set xlabel offset 0,0.5; set xtics offset 0.0,0.5
set ytics 1; set mytics 2; set ylabel offset 0.2,0.0; set ytics offset 0.8,0.0
set cbtics 0.1; set mcbtics 2; set cblabel offset 0.3,0.0; set cbtics offset -0.6,0.0

set pm3d interpolate 0,0
set view map
#set palette red
#set palette gray negative
unset key
set size ratio 0.9
set label 1 "" at 0.0,0.0 front point pointtype 7 lc 3 lw 2 ps 1
set label 2 "Star" at 1.0,0.0 font "Helvetica,21" front textcolor lt 3
set label 2 "{/ZapfDingbats \110}" at -0.05,0.05 front textcolor rgb '#008B8B' font "Helvetica,28"
set arrow 4 from -0.01,0.1 to -0.17,0.46 head filled front linecolor rgb '#6495ED' linewidth 2.0

#set object 5 ellipse center graph 0.91,0.91 size 1./6.,1./6. angle 0. fc rgb "grey" fs solid 1.0 border -1 front

#splot '../SW_spex_6_10.dat' u ($1/60.):(-$2/60.):($3) w pm3d,\

splot '../SW_spex_6_10.dat' u (($1*cos(ang)+$2*sin(ang))/60.):(($1*sin(ang)-$2*cos(ang))/60.):($3) w pm3d,\
      '../R0_proj.dat' u (($1*cos(ang)-$2*sin(ang))/60.):(($1*sin(ang)+$2*cos(ang))/60.):1 ps 3 lw 2 pt 1 lc rgb "gold"

set terminal postscript colour enhanced dl 3.5 lw 2. "Helvetica" 21
set output 'spex_6_10_map.ps'
replot
!./fixbb spex_6_10_map.ps
