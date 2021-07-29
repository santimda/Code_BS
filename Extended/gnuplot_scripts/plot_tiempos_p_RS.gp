set encoding iso_8859_1
set size ratio 0.9
set title ''#'Proton cooling time (RS)'
set key out 
set key box opaque height 1.01 width 1.0 
set key horiz samplen 2. spacing 1.01 font ",20"
set xrange [9:13.]
set yrange [3:16.]
set ylabel 'log_{10}(t [s])'; set ylabel offset 0.9,0
set xlabel 'log_{10}(E_p [eV])'; set xlabel offset 0,0.2
set ytics 2; set mytics 2
set xtics 1; set mxtics 2
set label "NAME" at graph 0.6,0.95 font "Helvetica,18"

plot "< awk '$1==20 { print $2, $3  }' ../CD_tcoolp_RS.dat" u (log10($1)):(log10($2)) t 'ac' w l lw 2.5 lc 8 dt 4,\
    "< awk '$1==20 { print $2, $8  }' ../CD_tcoolp_RS.dat" u (log10($1)):(log10($2)) t 'cell' w l lw 2.5 dt 3 lc 7,\
    "< awk '$1==20 { print $2, $5  }' ../CD_tcoolp_RS.dat" u (log10($1)):(log10($2)) t 'conv' w l lw 2.5 lc rgb "coral" dt 6,\
    "< awk '$1==20 { print $2, $6  }' ../CD_tcoolp_RS.dat" u (log10($1)):(log10($2)) t 'diff' w l lw 2.5 lc rgb "#32CD32" dt 6,\
    "< awk '$1==20 { print $2, $11 }' ../CD_tcoolp_RS.dat" u (log10($1)):(log10($2)) t 'adi' w l lw 2.5 lc 3,\
    "< awk '$1==20 { print $2, $4  }' ../CD_tcoolp_RS.dat" u (log10($1)):(log10($2)) t 'pp' w l lw 2.5 lc 1

#     "< awk '$1==10 { print $2, $3, $4, $5, $6, $7, $8, $9 }' ../CD_tcoolp_RS.dat" u (log10($1)):(log10($6)) t 't_{cool}' w l lw 3, \
#     "< awk '$1==150 { print $2, $3, $4, $5, $6, $7, $8, $9 }' ../CD_tcoolp_RS.dat" u (log10($1)):(log10($8)) t 't_{esc}' w l lw 3, \

set terminal postscript color solid enhanced dl 1.25 lw 2.5 "Helvetica" 22
set output 't_p_RS.eps'
replot
!./fixbb t_p_RS.eps
!epstopdf t_p_RS.eps
