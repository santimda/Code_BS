eV=1.602176e-12
set encoding utf8    #iso_8859_1
set title ''
#'Electron distribution in FS'
set xrange [7:14.]
set yrange [31.:]
#set cbrange [1:2000]
set xlabel 'log_{10}(E_e [eV])'; set ylabel offset 0.0,-0.1
set ylabel 'log_{10}({E_e}^2 N_e(E_e) [erg])'; set ylabel offset 0.1,0.0
set cblabel '{/Symbol q}/{/Symbol p}'; set cblabel offset -0.5,-1.1
set ytics 1 ; set mytics 2
set mxtics 2
set key off
set datafile missing
set label "NAME" at graph 0.05,0.95 font "Helvetica,18"

plot '../diste_FS.dat' u (log10($2)):(log10((($2)*eV)**2*($3))):($1/pi) w l lw 2.0 lc palette,\
     '../diste_FS_total.dat' u (log10($1)):(log10((($1)*eV)**2*($2))) w l lw 3.0 dt 2	 

#plot "< awk '($1+0.0) > 2.e-4 { print $1, $2, $3 }' diste_FS.dat" u (log10($2)):(log10((($2)*eV)**2*($3))):($1) w l lc palette,\
     '../diste_FS_total.dat' u (log10($1)):(log10((($1)*eV)**2*($2))) w l lw 3 dt 2

#plot "< awk '($1+0.0) > 1.83 && $1 < 1.86 { print $1, $2, $3 }' ../diste_FS.dat" u (log10($2)):(log10((($2)*eV)**2*($3))):($1/pi) w l lc palette

set terminal postscript color solid enhanced dl 1.5 lw 2.5 "Helvetica" 22
set output 'diste_FS.eps'
replot
!./fixbb diste_FS.eps
!epstopdf diste_FS.eps