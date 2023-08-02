set encoding iso_8859_1
set key out
set size ratio 0.6
set key box opaque height 0.6 #width 0.50 
#set key off
set key horiz samplen 3. spacing 1.4 font ",20"
set xrange [-7:14]
set yrange [26:35]
set xtics 3 font "Helvetica,18"; set mxtics 3 
set ytics 2 font "Helvetica,18"; set mytics 2  
set xlabel 'log_{10} ({/Symbol e} [eV])'
set ylabel 'log_{10} ({/Symbol e} L({/Symbol e}) [erg s^{-1}])'
set xlabel offset 0.0,0.3
set ylabel offset 1.,0.0
#set bmargin 0.2
set pointsize 1.5
set size ratio 0.7
set label "NAME" at graph 0.05,0.95 font "Helvetica,18"
#set label "Fermi" at 6.,31.2
#set label "CTA" at 10.,30.8

set terminal postscript color solid enhanced dl 2.5 lw 1.7 "Helvetica" 20
set output 'SED_th_NT.eps'

#-------------------------------------------------------------
# CONSTANTS:

# 1 Jy = 10^-23 erg/cm2/s/Hz, nu en $1 está en GHz y Sv en $2 en mJy
dist = 1.*3.0856e21
eV = 1.602176e-12
a = 4.0*pi*dist**2.
erg = 6.6261E-27 # Hz to erg
cor = log10((1.72/1.4)**2)

plot '../radsyn_RS_tot.dat' u (log10($1)):(log10($2)) t 'sy' w l lw 3 dt 1 lc "red",\
     '../radbr_RS_tot.dat' u (log10($1)):(log10($2)) t 'Br' w l lw 3 lc rgb '#005A32' smooth bezier,\
     '../radic_st_RS_tot.dat' u (log10($1)):(log10($2)) t 'IC_*' w l lw 3 lc rgb "blue" smooth bezier,\
     '../radic_ir_RS_tot.dat' u (log10($1)):(log10($2)) t 'IC_{IR}' w l lw 3 lc rgb 'violet',\
     '../radpp_RS_tot.dat' u (log10($1)):(log10($2)) t 'p-p' w l lw 3 lc rgb 'orange' smooth bezier,\
     '../FS_ff.dat' u (log10($1)):(log10($2)) t 'HL' w l lw 3 dt 3 lc "black",\
     '../ionized_ff_FS.dat' u (log10($1)):(log10($2)) t 'IL' w l lw 3 dt 4 lc "black",\
     '../RS_ff.dat' u (log10($1)):(log10($2)) t 'th-RS' w l lw 3 dt 1 lc "black",\
     './data_XMM_UL.dat' u ($1):($2+cor):($3):($4) w xerrorbars pt 1 lc rgb 'purple' t '',\
     './data_XMM_UL.dat' u ($1):($2+cor):(0):(-0.6) t '' with vector lt 1 lw 2.0 lc rgb 'purple',\
     './sens_CTA_100h.dat' u (log10(1.e6*($1))):(log10(($2)*a)) t 'CTA' w l lt 3 lc rgb "black" lw 4,\
     './sens_HESS_100h.dat' u (log10($1*1.e6)):(log10(($2)*a)) t 'HESS' w l lc 'grey40' lw 4,\
     './data_XMM_UL.dat' u ($1):($2+cor):($3):($4) w xerrorbars pt 1 lc rgb 'purple' t '',\
     './data_XMM_UL.dat' u ($1):($2+cor):(0):(-0.6) t '' with vector lt 1 lw 2.0 lc rgb 'purple',\
     './data_Fermi_UL.dat' u (log10($1*1.e9)):(log10($2*a*eV)):(log10($3*1.e9)):(log10($4*1.e9)) t "" w xerrorbars pt 1 lc rgb 'orange',\
     './data_Fermi_UL.dat' u (log10($1*1.e9)):(log10($2*a*eV)):(0):(-0.6) t '' with vector lt 1 lw 2.0 lc rgb 'orange'

#     '../radpp_FS_tot.dat' u (log10($1)):(log10($2)) t '' w l lw 3 dt 3 lc rgb 'orange' smooth bezier,\
#     '../radic_ir_FS_tot.dat' u (log10($1)):(log10($2)) t '' w l lw 3 dt 3 lc rgb 'violet',\
#     '../radic_st_FS_tot.dat' u (log10($1)):(log10($2)) t '' w l lw 3 dt 3 lc rgb "blue" smooth bezier,\
#     '../radbr_FS_tot.dat' u (log10($1)):(log10($2)) t '' w l lw 3 dt 3 lc rgb '#005A32' smooth bezier,\
#     '../radsyn_FS_tot.dat' u (log10($1)):(log10($2)) t '' w l lw 3 dt 3 lc "red",\
#     './sens_chandra_100ks.dat' u (log10($1)):(log10($2*a)) t 'Chandra' w l lt 2 lc rgb 'black' lw 2
#     './sens_fermi10yr_120_45.dat' u (log10($1*1.e6)):(log10($2*a)) t 'Fermi' w l lt 2 lc rgb 'grey40' lw 4,\
#     './sens_CTA_100h.dat' u (log10($1*1.e6)):(log10($2*a)) t 'CTA' w l lt 2 lc rgb 'black' lw 4,\
#     './sens_eastrogam.dat' u (log10($1*1.e6)):(log10($2*a)) t 'e-ASTROGAM' w l lt 2 lc rgb 'green' lw 4,\
#     './sens_HESS_100h.dat' u (log10($1*1.e6)):(log10($2*a)) t 'H.E.S.S.' w l lt 2 lc 1 lw 4
     

   #  '../radic_CMB_FS.dat' u (log10($1)):(log10($2)) t 'IC-CMB (FS)' w l lw 2 lc 4 dt 2 smooth bezier,\
   #  '../radic_CMB_RS.dat' u (log10($1)):(log10($2)) t 'IC-CMB (RS)' w l lw 2 lc 4 smooth bezier,\

!./fixbb SED_th_NT.eps
!epstopdf SED_th_NT.eps