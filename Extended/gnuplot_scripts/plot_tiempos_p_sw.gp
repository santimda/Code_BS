set encoding iso_8859_1
set title ''#'Proton cooling time (SW)'
set key out 
set key box opaque height 1.05 width 1.0 
set key horiz samplen 3. spacing 1.1 font ",20"
set xrange [9:15.]
set yrange [5:17.5]
set ylabel 'log_{10}(t [s])'; set ylabel offset 0.9,0
set xlabel 'log_{10}(E_p [eV])'; set xlabel offset 0,0.2
set ytics 2; set mytics 2

plot "< awk '$1==90 { print $2, $3  }' ../CD_tcoolp.dat" u (log10($1)):(log10($2)) t 't_{ac}' w l lw 2.5 lc 8 dt 4, \
     "< awk '$1==90 { print $2, $5  }' ../CD_tcoolp.dat" u (log10($1)):(log10($2)) t 't_{conv}' w l lw 2.5 lc rgb "coral" dt 6, \
     "< awk '$1==90 { print $2, $6  }' ../CD_tcoolp.dat" u (log10($1)):(log10($2)) t 't_{diff}' w l lw 2.5 lc rgb "#32CD32" dt 6, \
     "< awk '$1==90 { print $2, $11 }' ../CD_tcoolp.dat" u (log10($1)):(log10($2)) t 't_{adi}' w l lw 2.5 lc 3, \
     "< awk '$1==90 { print $2, $4  }' ../CD_tcoolp.dat" u (log10($1)):(log10($2)) t 't_{pp}' w l lw 2.5 lc 1

#     "< awk '$1==90 { print $2, $8  }' ../CD_tcoolp.dat" u (log10($1)):(log10($2)) t 't_{cell}' w l lw 2.5 dt 3 lc 7
#     "< awk '$1==10 { print $2, $3, $4, $5, $6, $7, $8, $9 }' ../CD_tcoolp.dat" u (log10($1)):(log10($6)) t 't_{cool}' w l lw 3, \
#     "< awk '$1==90 { print $2, $3, $4, $5, $6, $7, $8, $9 }' ../CD_tcoolp.dat" u (log10($1)):(log10($8)) t 't_{esc}' w l lw 3, \

set terminal postscript color solid enhanced dl 1.25 lw 2.5 "Helvetica" 22
set output 't_p_sw.eps'
replot
!./fixbb t_p_sw.eps
