eV=1.602176e-12
set encoding iso_8859_1
set title '' #'Proton distribution in RS'
set xrange [9:14.9]
set yrange [40:]
set xlabel 'log_{10}(E_p [eV])'; set ylabel offset 0.0,-0.1
set ylabel 'log_{10}({E_p}^2 N_p(E_p) [erg])'; set ylabel offset 0.1,0.0
set ytics 1 ; set mytics 2
set mxtics 2
set palette maxcolors 5 # 5 = N_bins
set cblabel '' #'{/Symbol q}/{/Symbol p}'; set cblabel offset -0.5,-0.8
set cbtics ("I1" 1.4, "I2" 2.2,"I3" 3.,"I4" 3.8,"I5" 4.6)
set palette defined ( 0 "#FF8C00", 1 "#8B0000" ) # From dark-orange to dark-red
set datafile missing


plot '../distp_RS_bin.dat' u (log10($2)):(log10((($2)*eV)**2*($3))):1 w l lw 3.0 dt 1 lc palette t '',\
     '../distp_RS_total.dat' u (log10($1)):(log10((($1)*eV)**2*($2))) w l lw 3.0 dt 2 lc 0 t 'Total'	 

set terminal postscript color solid enhanced dl 1.5 lw 2.5 "Helvetica" 21
set output 'distp_RS_bin.eps'
replot
!./fixbb distp_RS_bin.eps
