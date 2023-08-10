set view map

ang_deg = 3.3
ang = ang_deg*pi/180.  # rotation angle to match observed map (Ben2020)
cosd = 0.7195 # cos(dec) = cos(43.98*pi/180.) , for projection

set contour base
set isosamples 250
set cntrparam bspline
set cntrparam levels discrete 1.04,1.5,2,3,4,5,6.0
#set cntrparam levels discrete 60.6530659713#, 34.1, 0.0456, 0.0026
unset surface
set table 'cont.dat'
splot '../SW_conv_3GHzabs.dat' u (($1*cos(ang)+$2*sin(ang))/60.):(($1*sin(ang)-$2*cos(ang))/60.):($3)
#splot '../SW_conv_3GHzabs.dat' u (($1*cos(ang)+$2*sin(ang))):(($1*sin(ang)-$2*cos(ang))):($3)
unset table

set contour base
set isosamples 250
set cntrparam bspline
set cntrparam levels discrete 0.3 #rms?
#set cntrparam levels discrete 50.
unset surface
set table 'cont_red.dat'
splot '../SW_conv_3GHzabs.dat' u (($1*cos(ang)+$2*sin(ang))/60.):(($1*sin(ang)-$2*cos(ang))/60.):($3)
#splot '../SW_conv_3GHzabs.dat' u (($1*cos(ang)+$2*sin(ang))):(($1*sin(ang)-$2*cos(ang))):($3)
unset table

set contour base
set isosamples 250
set cntrparam bspline
set cntrparam levels discrete 8.0 # ??
unset surface
set table 'cont_green.dat'
splot '../SW_conv_3GHzabs.dat' u (($1*cos(ang)+$2*sin(ang))/60.):(($1*sin(ang)-$2*cos(ang))/60.):($3)
unset table

reset
set encoding iso_8859_1
set title 'Synthetic map at 3.0 GHz'
set title ''
set xlabel '{/Symbol D}RA [arcmin]'
set ylabel '{/Symbol D}DEC [arcmin]'
#set xlabel '{/Symbol D}RA [arcsec]'
#set ylabel '{/Symbol D}DEC [arcsec]'

set cblabel 'Flux [mJy beam^{-1}]'
set yrange [-5:10]
set xrange [-10:10]
#set xrange [-20.:20.]
#set yrange [185:225]
#set cbrange [0.0:7]

set pm3d interpolate 0,0
set view map
#set palette gray negative
set palette defined (0 "white", 1 "grey40")
unset key

set size ratio 15./20.
set xtics 5; set mxtics 5; set xlabel offset 0,0.5; set xtics offset 0.0,0.5
set ytics 5; set mytics 5; set ylabel offset 0.2,0.0; set ytics offset 0.8,0.0
set cbtics 1; set mcbtics 2; set cblabel offset 0.3,0.0; set cbtics offset -0.6,0.0

set arrow 4 from -0.01,0.1 to -0.096,1.597 head filled front linecolor rgb '#6495ED' # BD+43
set label 2 "{/ZapfDingbats \110}" at -0.57,0.0 front textcolor rgb '#008B8B' font "Helvetica,28" # BD+43

#res = 12./60.
#coc = sqrt(8.*log(2.))
coc = 1.
#siz = res*coc
#set object 5 ellipse center graph 0.91,0.91 size siz/cosd,siz angle 0. fc rgb "grey" fs solid 1.0 border -1 front
siz_x = 12.5/60./coc
siz_y = 20.2/60./coc
set object 5 ellipse center graph 0.92,0.92 size siz_x/coc/cosd,siz_y/coc angle ang_deg fc rgb "grey" fs solid 1.0 border -1 front

set label "3 GHz" at graph 0.05,0.9 font "Helvetica,18" front

splot '../SW_conv_3GHzabs.dat' u (($1*cos(ang)+$2*sin(ang))/60.):(($1*sin(ang)-$2*cos(ang))/60.):($3) w pm3d,\
      'cont.dat' w l lt -1 lw 1.5,\
      '../R0_proj.dat' u (($1*cos(ang)-$2*sin(ang))/60.):(($1*sin(ang)+$2*cos(ang))/60.):1 ps 3 lw 2 pt 1 lc rgb "gold"

#splot '../SW_conv_3GHzabs.dat' u ($1*cos(ang)+$2*sin(ang)):($1*sin(ang)-$2*cos(ang)):($3) w pm3d,\
#      'cont.dat' w l lt -1 lw 1.5,\
#      'cont_red.dat' w l dt 2 lc 7 lw 1.5,\

#      'cont_green.dat' w l dt 2 lc rgb "#32CD32" lw 2,\
#      'cont_red.dat' w l dt 2 lc 7 lw 1.5,\

set terminal postscript colour enhanced dl 2. lw 2. "Helvetica" 22
set output 'convmap_3G_BD.eps'
replot
!./fixbb convmap_3G_BD.eps
!epstopdf convmap_3G_BD.eps
