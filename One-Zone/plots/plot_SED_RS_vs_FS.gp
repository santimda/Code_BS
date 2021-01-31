set encoding iso_8859_1
set key box out horizontal top
set key font ",10"
set xrange [-7:14]
set yrange [20:28]
set xtics 3; set mxtics 3
set ytics 2; set mytics 2  
set xlabel 'log_{10} ({/Symbol e} [eV])'
set ylabel 'log_{10} ({/symbol e} L({/symbol e}) [erg s^{-1}])'
set pointsize 1.5

plot '../radsyn_FS.dat' u (log10($1)):(log10($3)) t 'sy (FS)' w l lw 2 dt 2 lc 1,\
     '../radsyn_RS.dat' u (log10($1)):(log10($2)) t 'sy (RS)' w l lw 2 dt 1 lc 1,\
     '../radic_IR_FS.dat' u (log10($1)):(log10($2)) t 'IC-dust (FS)' w l lw 2 dt 2 lc 2,\
     '../radic_IR_RS.dat' u (log10($1)):(log10($2)) t 'IC-dust (RS)' w l lw 2 lc 2,\
     '../radic_st_FS.dat' u (log10($1)):(log10($2)) t 'IC-star (FS)' w l lw 2 dt 2 lc 3 smooth bezier,\
     '../radic_st_RS.dat' u (log10($1)):(log10($2)) t 'IC-star (RS)' w l lw 2 lc 3 smooth bezier,\
     '../radic_CMB_FS.dat' u (log10($1)):(log10($2)) t 'IC-CMB (FS)' w l lw 2 lc 4 dt 2 smooth bezier,\
     '../radic_CMB_RS.dat' u (log10($1)):(log10($2)) t 'IC-CMB (RS)' w l lw 2 lc 4 smooth bezier,\
     '../radbr_FS.dat' u (log10($1)):(log10($2)) t 'rel. br. (FS)' w l lw 2 dt 2 lc 5 smooth bezier,\
     '../radbr_RS.dat' u (log10($1)):(log10($2)) t 'rel. br. (RS)' w l lw 2 lc 5 smooth bezier,\
     '../radpp_FS.dat' u (log10($1)):(log10($2)) t 'p-p' w l lw 2 dt 2 lc 6,\
     '../radpp_RS.dat' u (log10($1)):(log10($2)) t 'p-p' w l lw 2 lc 6

set terminal postscript color solid enhanced dl 2.5 lw 2.2 "Helvetica" 21
set output 'SED.ps'
replot
!./fixbb SED.ps
