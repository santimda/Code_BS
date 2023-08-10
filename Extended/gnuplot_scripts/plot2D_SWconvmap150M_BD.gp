set view map

ang_deg = 3.3
#ang_deg = 8.
ang = ang_deg*pi/180.  # rotation angle to match observed map (Ben2020)
cosd = 0.7195 # cos(dec) = cos(43.98*pi/180.) , for projection
#cosd = 0.4819 # cos(dec) = cos(61.19*pi/180.) , for projection
sind = 0.6944

set contour base
set isosamples 250
set cntrparam bspline
set cntrparam levels discrete 10.
unset surface
set table 'cont.dat'
splot '../SW_conv_150MHzabs.dat' u (($1*cos(ang)+$2*sin(ang))/60.):(($1*sin(ang)-$2*cos(ang))/60.):($3)
unset table

set contour base
set isosamples 250
set cntrparam bspline
set cntrparam levels discrete 20., 40., 80, 120.
unset surface
set table 'cont_red.dat'
splot '../SW_conv_150MHzabs.dat' u (($1*cos(ang)+$2*sin(ang))/60.):(($1*sin(ang)-$2*cos(ang))/60.):($3)
unset table

set contour base
set isosamples 250
set cntrparam bspline
set cntrparam levels discrete 85.0 # ??
unset surface
set table 'cont_green.dat'
splot '../SW_conv_150MHzabs.dat' u (($1*cos(ang)+$2*sin(ang))/60.):(($1*sin(ang)-$2*cos(ang))/60.):($3)
unset table

reset
set encoding iso_8859_1
set title 'Synthetic map at 3.0 GHz'
set title ''
set xlabel '{/Symbol D}RA [arcmin]'
set ylabel '{/Symbol D}DEC [arcmin]'

set cblabel 'Flux [mJy beam^{-1}]'
set xrange [-10:10]
set yrange [-5:10]
set cbrange [5.0:]

set pm3d interpolate 0,0
set view map
#set palette gray negative
set palette defined (0 "white", 1 "grey40")
unset key

set size ratio 15./20.
set xtics 5; set mxtics 5; set xlabel offset 0,0.5; set xtics offset 0.0,0.5
set ytics 5; set mytics 5; set ylabel offset 0.2,0.0; set ytics offset 0.8,0.0
#set cbtics 10; set mcbtics 2; set cblabel offset 0.3,0.0; set cbtics offset -0.6,0.0

set arrow 4 from -0.01,0.1 to -0.096,1.597 head filled front linecolor rgb '#6495ED' # BD+43
set label 2 "{/ZapfDingbats \110}" at -0.57,0.0 front textcolor rgb '#008B8B' font "Helvetica,28" # BD+43

set label "150 MHz" at graph 0.05,0.9 font "Helvetica,18" front

#siz_x = 4./60.#*coc
#siz_y = 16./60.#*coc
#siz_x = 25./60.
#siz_y = 25./60.
#siz_x = 4./60.
#siz_y = 4./60.
#set object 5 ellipse center graph 0.92,0.92 size siz_x/cosd,siz_y angle ang_deg fc rgb "grey" fs solid 1.0 border -1 front

theta_x = 25./60.
theta_y = 25./60.  # The FWHM
#coc = sqrt(8.*log(2.))  # Conversion factor between FWHM and sigma
coc = 1.
siz_x = theta_x/coc   # The sigma
siz_y = theta_y/coc
set object 5 ellipse center graph 0.92,0.92 size siz_x,siz_y angle ang_deg fc rgb "grey" fs solid 1.0 border -1 front

set label "150 MHz" at graph 0.05,0.9 font "Helvetica,18" front

splot '../SW_conv_150MHzabs.dat' u (($1*cos(ang)+$2*sin(ang))/60.):(($1*sin(ang)-$2*cos(ang))/60.):($3) w pm3d,\
      '../R0_proj.dat' u (($1*cos(ang)-$2*sin(ang))/60.):(($1*sin(ang)+$2*cos(ang))/60.):1 ps 3 lw 2 pt 1 lc rgb "gold",\
      'cont.dat' w l lt -1 lw 1.5,\
      'cont_red.dat' w l lt -1 lc 7 lw 1.5,\
    
#     'cont_green.dat' w l dt 2 lc rgb "#32CD32" lw 2,\

set terminal postscript colour enhanced dl 2. lw 2. "Helvetica" 22
set output 'convmap_150M_BD.eps'
replot
!./fixbb convmap_150M_BD.eps
!epstopdf convmap_150M_BD.eps