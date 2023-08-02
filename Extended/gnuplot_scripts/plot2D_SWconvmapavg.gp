set view map

ang_deg = 3.3 # BD+43
#ang_deg = 8. # BD+60
ang = ang_deg*pi/180.  # rotation angle to match observed map
#cosd = 0.4819 # cos(dec) = cos(61.19*pi/180.) , for projection
cosd = 0.7195 # cos(dec) = cos(43.98*pi/180.) , for projection

set contour base
set isosamples 250
set cntrparam bspline
#set cntrparam levels discrete 0.5,1.04,1.5,3 # BD+60
set cntrparam levels discrete 0.3,0.6,0.9 # BD+43
unset surface
set table 'cont.dat'
splot '../SW_conv_avg_abs.dat' u (($1*cos(ang)+$2*sin(ang))/60.):(($1*sin(ang)-$2*cos(ang))/60.):($3)
unset table

set contour base
set isosamples 250
set cntrparam bspline
#set cntrparam levels discrete 0.1 #rms?
set cntrparam levels discrete 1.8 # BD+43
unset surface
set table 'cont_red.dat'
splot '../SW_conv_avg_abs.dat' u (($1*cos(ang)+$2*sin(ang))/60.):(($1*sin(ang)-$2*cos(ang))/60.):($3)
unset table

set contour base
set isosamples 250
set cntrparam bspline
set cntrparam levels discrete 8.0 # ??
unset surface
set table 'cont_green.dat'
splot '../SW_conv_avg_abs.dat' u (($1*cos(ang)+$2*sin(ang))/60.):(($1*sin(ang)-$2*cos(ang))/60.):($3)
unset table

reset
set encoding iso_8859_1
set title 'Synthetic avg map at 4-12 GHz'
set title ''
set xlabel '{/Symbol D}RA [arcmin]'
set ylabel '{/Symbol D}DEC [arcmin]'

set cblabel 'Flux [mJy beam^{-1}]'  
set xrange [-10:10] # BD+43
set yrange [-5:10] # BD+43
#set xrange [-4.:4.] # BD+60
#set yrange [-3.:3.] # BD+60 
#set cbrange [0.0:1]
#set cbrange [0.0:4]

set pm3d interpolate 0,0
set view map
#set palette gray negative
set palette defined (0 "white", 1 "grey40")
unset key

#set size ratio 7./13.#7./9.36 # (1.5+5.5)/(cosd*(6.5+6.5))
set size ratio 15./20. # BD+43
#set size ratio 6./8. # BD+60
#set xtics 5; set mxtics 5; set xlabel offset 0,0.5; set xtics offset 0.0,0.5
#set ytics 5; set mytics 5; set ylabel offset 0.2,0.0; set ytics offset 0.8,0.0
set xtics 2; set mxtics 2; set xlabel offset 0,0.5; set xtics offset 0.0,0.5
set ytics 2; set mytics 2; set ylabel offset 0.2,0.0; set ytics offset 0.8,0.0
set cbtics 1; set mcbtics 2; set cblabel offset 0.3,0.0; set cbtics offset -0.6,0.0

set arrow 4 from -0.01,0.1 to -0.096,1.597 head filled front linecolor rgb '#6495ED' # BD+43
set label 2 "{/ZapfDingbats \110}" at -0.57,0.0 front textcolor rgb '#008B8B' font "Helvetica,28" # BD+43
#set arrow 4 from -0.01,0.0 to -0.09,0.49 head filled front linecolor rgb '#6495ED' # BD+60
#set label 2 "{/ZapfDingbats \110}" at -0.2,0.0 front textcolor rgb '#008B8B' font "Helvetica,20" # BD+60

#set arrow 4 from -0.01,0.1 to -0.096,1.597 head filled front linecolor rgb '#6495ED' linewidth 2.0
#set label 2 "{/ZapfDingbats \110}" at -0.7,0.0 front textcolor rgb '#008B8B' font "Helvetica,28"

set label "avg 4-12 MHz" at graph 0.05,0.9 font "Helvetica,18" front

#res = 12./60.
theta_x = 19.7 # BD+43
theta_y = 8.25 # BD+43
#theta_x = 16. # BD+60
#theta_y = 10. # BD+60
#coc = sqrt(8.*log(2.))
coc = 1.
#siz = res*coc
#siz_x = 9.8/60.#*coc
#siz_y = 18.3/60.#*coc
siz_x = theta_x/coc/60.
siz_y = theta_y/coc/60.
#set object 5 ellipse center graph 0.92,0.92 size siz_x/cosd,siz_y angle ang_deg fc rgb "grey" fs solid 1.0 border -1 front
set object 5 ellipse center graph 0.92,0.92 size siz_x,siz_y angle ang_deg fc rgb "grey" fs solid 1.0 border -1 front

splot '../SW_conv_avg_abs.dat' u (($1*cos(ang)+$2*sin(ang))/60.):(($1*sin(ang)-$2*cos(ang))/60.):($3) w pm3d,\
      'cont.dat' w l lt -1 lw 1.5,\
      '../R0_proj.dat' u (($1*cos(ang)-$2*sin(ang))/60.):(($1*sin(ang)+$2*cos(ang))/60.):1 ps 3 lw 2 pt 1 lc rgb "gold",\
      'cont_red.dat' w l lt -1 lc 7 lw 1.5,\

#      'cont_green.dat' w l dt 2 lc rgb "#32CD32" lw 2,\

set terminal postscript colour enhanced dl 2. lw 2. "Helvetica" 22
set output 'convmap_avg.eps'
replot
!./fixbb convmap_avg.eps
!epstopdf convmap_avg.eps