set encoding iso_8859_1
set key out
set key box opaque height 0.6 #width 0.50 
set key horiz samplen 3. spacing 1.01 font ",20"
set xrange [-7:14]
set yrange [22.:30.]
set xtics 3 font "Helvetica,18"; set mxtics 3 
set ytics 2 font "Helvetica,18"; set mytics 2  
set xlabel 'log_{10} ({/Symbol e} [eV])'
set ylabel 'log_{10} ({/Symbol e} L({/Symbol e}) [erg s^{-1}])'
set xlabel offset 0.0,0.3
set ylabel offset 1.,0.0
#set bmargin 0.2
set pointsize 1.5
set size ratio 0.8
set label "NAME" at graph 0.05,0.95 font "Helvetica,18"

set terminal postscript color solid enhanced dl 2.5 lw 1.7 "Helvetica" 20
set output 'SED_RS_vs_FS.eps'

plot '../radsyn_FS_tot.dat' u (log10($1)):(log10($2)) t '' w l lw 3 dt 2 lc "red",\
     '../radsyn_RS_tot.dat' u (log10($1)):(log10($2)) t 'sy' w l lw 3 dt 1 lc "red",\
     '../radic_ir_FS_tot.dat' u (log10($1)):(log10($2)) t '' w l lw 3 dt 2 lc rgb "cyan",\
     '../radic_ir_RS_tot.dat' u (log10($1)):(log10($2)) t 'IC_{IR}' w l lw 3 lc rgb "cyan",\
     '../radic_st_FS_tot.dat' u (log10($1)):(log10($2)) t '' w l lw 3 dt 2 lc "magenta" smooth bezier,\
     '../radic_st_RS_tot.dat' u (log10($1)):(log10($2)) t 'IC_{st}' w l lw 3 lc "magenta" smooth bezier,\
     '../radbr_FS_tot.dat' u (log10($1)):(log10($2)) t '' w l lw 3 dt 2 lc "gold" smooth bezier,\
     '../radbr_RS_tot.dat' u (log10($1)):(log10($2)) t 'rel. br.' w l lw 3 lc "gold" smooth bezier,\
     '../radpp_FS_tot.dat' u (log10($1)):(log10($2)) t '' w l lw 3 dt 2 lc rgb '#005A32' smooth bezier,\
     '../radpp_RS_tot.dat' u (log10($1)):(log10($2)) t 'p-p' w l lw 3 lc rgb '#005A32' smooth bezier

   #  '../radic_CMB_FS.dat' u (log10($1)):(log10($2)) t 'IC-CMB (FS)' w l lw 2 lc 4 dt 2 smooth bezier,\
   #  '../radic_CMB_RS.dat' u (log10($1)):(log10($2)) t 'IC-CMB (RS)' w l lw 2 lc 4 smooth bezier,\

!./fixbb SED_RS_vs_FS.eps
!epstopdf SED_RS_vs_FS.eps