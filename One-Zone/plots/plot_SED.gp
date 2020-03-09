set encoding iso_8859_1
set key top right
set xrange [-6:14]
set yrange [25:34]
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

plot './data_radio2010.dat' u (log10(($1)*10**9*h/eV)):(log10(($1)*10**9*a*($2)*10**(-26))):(log10(1+($3)/($1))):(log10(1.+($4)/($2))) w xyerrorbars pt 7 t 'VLA', \
     '../radsyn_abs.dat' u (log10($1)):(log10($3)) t 'sync' w l lw 2 lc 1,\
     '../radsyn_abs.dat' u (log10($1)):(log10($2)) t '' w l lw 2 dt 2 lc 1,\
     '../radic_IR.dat' u (log10($1)):(log10($2)) t 'IC-dust' w l lw 2,\
     '../radic_st.dat' u (log10($1)):(log10($2)) t 'IC-star' w l lw 2 smooth bezier,\
     '../radbr.dat' u (log10($1)):(log10($2)) t 'rel. br.' w l lw 2 smooth bezier,\
     '../radpp.dat' u (log10($1)):(log10($2)) t 'p-p' w l lw 2 #smooth bezier

set terminal postscript color solid enhanced dl 2.5 lw 2.2 "Helvetica" 21
set output 'SED.ps'
replot
!./fixbb SED.ps
