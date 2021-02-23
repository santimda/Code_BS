set encoding iso_8859_1
set key box out horizontal top
set key font ",18"
set xrange [-7:14]
set yrange [22:29]
set xtics 3 font "Helvetica,18"; set mxtics 3 
set ytics 2 font "Helvetica,18"; set mytics 2  
set xlabel 'log_{10} ({/Symbol e} [eV])'
set ylabel 'log_{10} ({/Symbol e} L({/Symbol e}) [erg s^{-1}])'
set bmargin 0.2
set pointsize 1.5
set label "NAME" at -6.5,28.5 font "Helvetica,18"

set terminal postscript color solid enhanced dl 2.5 lw 1.7 "Helvetica" 20
set output 'SED_RS_vs_FS.eps'

plot '../radsyn_FS.dat' u (log10($1)):(log10($3)) t 'sy (FS)' w l lw 2 dt 2 lc "red",\
     '../radsyn_RS.dat' u (log10($1)):(log10($2)) t 'sy (RS)' w l lw 2 dt 1 lc "red",\
     '../radic_IR_FS.dat' u (log10($1)):(log10($2)) t 'IC-dust (FS)' w l lw 2 dt 2 lc 2,\
     '../radic_IR_RS.dat' u (log10($1)):(log10($2)) t 'IC-dust (RS)' w l lw 2 lc 2,\
     '../radic_st_FS.dat' u (log10($1)):(log10($2)) t 'IC-star (FS)' w l lw 2 dt 2 lc "magenta" smooth bezier,\
     '../radic_st_RS.dat' u (log10($1)):(log10($2)) t 'IC-star (RS)' w l lw 2 lc "magenta" smooth bezier,\
     '../radbr_FS.dat' u (log10($1)):(log10($2)) t 'rel. br. (FS)' w l lw 2 dt 2 lc 5 smooth bezier,\
     '../radbr_RS.dat' u (log10($1)):(log10($2)) t 'rel. br. (RS)' w l lw 2 lc 5 smooth bezier,\
     '../radpp_FS.dat' u (log10($1)):(log10($2)) t 'p-p' w l lw 2 dt 2 lc 6,\
     '../radpp_RS.dat' u (log10($1)):(log10($2)) t 'p-p' w l lw 2 lc 6

   #  '../radic_CMB_FS.dat' u (log10($1)):(log10($2)) t 'IC-CMB (FS)' w l lw 2 lc 4 dt 2 smooth bezier,\
   #  '../radic_CMB_RS.dat' u (log10($1)):(log10($2)) t 'IC-CMB (RS)' w l lw 2 lc 4 smooth bezier,\


!./fixbb SED_RS_vs_FS.eps
!epstopdf SED_RS_vs_FS.eps
