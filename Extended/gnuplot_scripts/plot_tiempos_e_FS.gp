set encoding iso_8859_1
#set size ratio 0.65
set title '' #'Electron cooling time (FS)'
set key out 
set key box opaque height 1.01 width 0.50 
set key horiz samplen 2 spacing 1.01 font ",20"
#set key bottom right
set xrange [6:14]
set yrange [3:]
set ylabel 'log_{10}(t [s])'; set ylabel offset 0.9,0
set xlabel 'log_{10}(E_e [eV])'; set xlabel offset 0,0.2
set ytics 2; set mytics 2

plot "< awk '$1==100 { print $2, $7 }' ../CD_tcoole_FS.dat" u (log10($1)):(log10($2)) t 'ac' w l lw 2.5 lc 8 dt 4, \
    "< awk '$1==100 { print $2, $10 }' ../CD_tcoole_FS.dat" u (log10($1)):(log10($2)) t 'cell' w l lw 2.5 dt 3 lc 7,\
	"< awk '$1==100 { print $2, $9 }' ../CD_tcoole_FS.dat" u (log10($1)):(log10($2)) t 'conv' w l lw 2.5 lc rgb "coral" dt 6, \
    "< awk '$1==100 { print $2, $6 }' ../CD_tcoole_FS.dat" u (log10($1)):(log10($2)) t 'diff' w l lw 2.5 lc rgb "#32CD32" dt 6, \
    "< awk '$1==100 { print $2, $15}' ../CD_tcoole_FS.dat" u (log10($1)):(log10($2)) t 'adi' w l lw 2.5 lc 3, \
    "< awk '$1==100 { print $2, $3 }' ../CD_tcoole_FS.dat" u (log10($1)):(log10($2)) t 'IC' w l lw 2.5 lc 6, \
    "< awk '$1==100 { print $2, $4 }' ../CD_tcoole_FS.dat" u (log10($1)):(log10($2)) t 'sy' w l lw 2.5 lc rgb "#DC143C"


#     "< awk '$1==90 { print $2, $11 }' ../CD_tcoole_FS.dat" u (log10($1)):(log10($2)) t 't_{IC,st}' w l lw 2.5 lc 6, \
#     "< awk '$1==90 { print $2, $12 }' ../CD_tcoole_FS.dat" u (log10($1)):(log10($2)) t 't_{IC,IR}' w l lw 2.5 lc 4, \
#     "< awk '$1==90 { print $2, $5  }' ../CD_tcoole_FS.dat" u (log10($1)):(log10($2)) t 't_{br}' w l lw 2.5 lc 1
#     "< awk '$1==90 { print $2, $14 }' ../CD_tcoole_FS.dat" u (log10($1)):(log10($2)) t 't_{loss}' w l lw 2, \
#     "< awk '$1==90 { print $2, $8  }' ../CD_tcoole_FS.dat" u (log10($1)):(log10($2)) t 't_{cool}' w l lw 3, \
#     "< awk '$1==90 { print $2, $13 }' ../CD_tcoole_FS.dat" u (log10($1)):(log10($2)) t 't_{esc}' w l lw 3, \

set terminal postscript color solid enhanced dl 1.25 lw 2.5 "Helvetica" 22
set output 't_e_FS.eps'
replot
!./fixbb t_e_FS.eps
