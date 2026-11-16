set encoding iso_8859_1
set key bottom left
set xrange [-6:14]
set yrange [25:33]
set xtics 3; set mxtics 3
set ytics 2; set mytics 2  
set xlabel 'log_{10} ({/Symbol e} [eV])'
set ylabel 'log_{10} ({/Symbol e} L({/Symbol e}) [erg s^{-1}])'
set pointsize 1.5

# 1 Jy = 10^-23 erg/cm2/s/Hz, nu en $1 está en GHz y Sv en $2 en mJy
#dist = 1.45*3.0856e21
dist = 1.32*3.0856e21
h = 6.6260755e-27
eV = 1.602176e-12
a = 4.0*pi*dist**2.

plot './data_radio2010.dat' u 1:2:3:4:5:6 w xyerrorbars pt 7 lc 7 t 'VLA',\
     '../radsyn_RS.dat' u (log10($1)):(log10($2)) t '' w l lw 2 dt 2 lc 1,\
     '../radic_IR_RS.dat' u (log10($1)):(log10($2)) t 'IC-dust' w l lw 2,\
     '../radic_st_RS.dat' u (log10($1)):(log10($2)) t 'IC-star' w l lw 2 smooth bezier,\
     '../radic_CMB_RS.dat' u (log10($1)):(log10($2)) t 'IC-CMB' w l lw 2 smooth bezier,\
     '../radbr_RS.dat' u (log10($1)):(log10($2)) t 'rel. br.' w l lw 2 smooth bezier,\
     '../radpp_RS.dat' u (log10($1)):(log10($2)) t 'p-p' w l lw 2 #smooth bezier

set terminal postscript color solid enhanced dl 2.5 lw 2.2 "Helvetica" 21
set output 'SED.ps'
replot
!./fixbb SED.ps
