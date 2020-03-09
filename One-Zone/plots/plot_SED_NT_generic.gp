set encoding iso_8859_1
set key out
set key box opaque height 0.6 #width 0.50 
set key horiz samplen 3. spacing 1.01 font ",20"
set xrange [-6:14]
set yrange [25:33]
set xtics 3; set mxtics 3
set ytics 1; set mytics 2
set xlabel 'log_{10} ({/Symbol e} [eV])'
set ylabel 'log_{10} ({/Symbol e} L_{/Symbol e} [erg s^{-1}])'
set xlabel offset 0.0,0.3
set ylabel offset 1.,0.0

#-------------------------------------------------------------
# CONSTANTS:

# 1 Jy = 10^-23 erg/cm2/s/Hz, nu en $1 está en GHz y Sv en $2 en mJy
dist = 1.0*3.0856e21
h = 6.6260755e-27
eV = 1.602176e-12
a = 4.0*pi*dist**2.
fac = 1.0e9*h/eV

#-------------------------------------------------------------
# FITTING :

x0=10.0
f(x)=a0+a1*(x-x0)+a2*(x-x0)**2+a3*(x-x0)**3+a4*(x-x0)**4+a5*(x-x0)**5+a6*(x-x0)**6+a7*(x-x0)**7+a8*(x-x0)**8
fit f(x) "< awk '($1 >1.e3 && $1<5.e13) { print $1, $2, $3 }' ../radic_IR.dat" u (log10($1)):(log10($2)) via a0,a1,a2,a3,a4,a5,a6,a7,a8
f1(x)=x>13.7?0 : x<3.?0 : f(x)

x1=-1.0
g(x)=b0+b1*(x-x1)+b2*(x-x1)**2+b3*(x-x1)**3+b4*(x-x1)**4+b5*(x-x1)**5+b6*(x-x1)**6+b7*(x-x1)**7
fit g(x) "< awk '($1 >9.e-7 && $1<2.e4) { print $1, $2, $3 }' ../radsyn.dat" u (log10($1)):(log10($2)) via b0,b1,b2,b3,b4,b5,b6,b7
g1(x)=x>4.4?0 : x<-6?0 : g(x)

x2=8.5
j(x)=c0+c1*(x-x2)+c2*(x-x2)**2+c3*(x-x2)**3+c4*(x-x2)**4+c5*(x-x2)**5+c6*(x-x2)**6+c7*(x-x2)**7
fit j(x) "< awk '($1 >1.e3 && $1<9.e12) { print $1, $2, $3 }' ../radic_st.dat" u (log10($1)):(log10($2)) via c0,c1,c2,c3,c4,c5,c6,c7
j1(x)=x>12.8?0 : x<3?0 : j(x)

x3=11.0
k(x)=d0+d1*(x-x3)+d2*(x-x3)**2+d3*(x-x3)**3+d4*(x-x3)**4+d5*(x-x3)**5+d6*(x-x3)**6+d7*(x-x3)**7
fit k(x) "< awk '($1 >1.e8 && $1<1.e13) { print $1, $2, $3 }' ../radpp.dat" u (log10($1)):(log10($2)) via d0,d1,d2,d3,d4,d5,d6,d7
k1(x)=x>13?0 : x<8?0 : k(x)

#-----------------------------------------
#
#     PLOT
#
#----------------------------------------

plot log10(10**f1(x)+10**g1(x)+10**j1(x)+10**k1(x)) t 'Total' lw 3.3,\
     '../radsyn.dat' u (log10($1)):(log10($2)) t 'sync' w l lw 2 dt 2.4, \
     '../radic_IR.dat' u (log10($1)):(log10($2)) t 'IC-dust' w l lw 2.4 dt 3 smooth bezier, \
     '../radic_st.dat' u (log10($1)):(log10($2)) t 'IC-star' w l lw 2.4 dt 4 smooth bezier, \
     '../radbr.dat' u (log10($1)):(log10($2)) t 'rel. br.' w l lw 2.4 lc rgb "coral" dt 5 smooth bezier, \
     '../radpp.dat' u (log10($1)):(log10($2)) t 'p-p' w l lw 2.4 dt 6 smooth bezier, \
     './sens_CTA_100h.dat' u (log10(1.e6*($1))):(log10(($2)*a)) t 'CTA' w l lt 3 lc rgb "black" lw 4,\
     './sens_fermi10yr_120,45.dat' u (0.5*log10($1*1.e6*$2*1.e6)):(log10($3*a)) t 'Fermi' w l lt 3 lc rgb 'grey40' lw 4


#     './sens_hess.dat' u (log10($1)):(log10($2*a)) t 'HESS' w l lt 2 lc rgb 'grey60' lw 4,\
#     './data_radio2010.dat' u (log10(($1)*10**9*h/eV)):(log10(($1)*10**9*a*($2)*10**(-26))):(log10(1+($3)/($1))):(log10(1.+($4)/($2))) w xyerrorbars pt 7 t 'VLA',\
#     './sens_hawc.dat' u (log10(1.e9*($1))):(log10(($2)*1.6*a)) t 'HAWC' w l lc 'grey70' lw 4,\

set terminal postscript colour solid enhanced dl 1.25 lw 2.2 "Helvetica" 22
set output 'SED_NT_onezone_generic.eps'
replot
!./fixbb SED_NT_onezone_generic.eps
