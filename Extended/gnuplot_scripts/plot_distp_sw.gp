eV=1.602176e-12
set encoding iso_8859_1
set title ''
#'Proton distribution in SW'
set xrange [9:14.9]
set yrange [37.5:]
#set cbrange [1:2000]
set ylabel 'log_{10}({E_p}^2 N_p(E_p) [erg])'; set ylabel offset 0.1,0.0
set xlabel 'log_{10}(E_p [eV])'; set ylabel offset 0.0,-0.1
set cblabel '{/Symbol q}/{/Symbol p}'; set cblabel offset -0.5,-1.1
set key off
set ytics 1 ; set mytics 2
set mxtics 2
set datafile missing

plot '../distp_sw.dat' u (log10($2)):(log10((($2)*eV)**2*($3))):($1/pi) w l lw 2.0 lc palette,\
     '../distp_sw_total.dat' u (log10($1)):(log10(($1*eV)**2*($2))) w l lw 3.0 dt 2	 

#plot "< awk '$1 > 1.89e0 { print $1, $2, $3 }' distp_sw.dat" u (log10($2)):(log10((($2)*eV)**2*($3))):($1/pi) w l lc palette

set terminal postscript color solid enhanced dl 1.5 lw 2.5 "Helvetica" 22
set output 'distp_sw.eps'
replot
!./fixbb distp_sw.eps
