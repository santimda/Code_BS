eV=1.602176e-12
set encoding iso_8859_1
set title '' #'Proton distribution in FS'
set xrange [9.:13.]
set yrange [31.:41.]
set xlabel 'log_{10}(E_p [eV])'; set ylabel offset 0.0,-0.1
set ylabel 'log_{10}({E_p}^2 N_p(E_p) [erg])'; set ylabel offset 0.1,0.0
set ytics 1 ; set mytics 2
set xtics 0.5 ; set mxtics 1
set palette maxcolors 5 # 5 = N_bins
set cblabel '' #'{/Symbol q}/{/Symbol p}'; set cblabel offset -0.5,-0.8
set cbtics ("I1" 1.4, "I2" 2.2,"I3" 3.,"I4" 3.8,"I5" 4.6)
set palette defined ( 0 "#FF8C00", 1 "#8B0000" ) # From dark-orange to dark-red
set datafile missing
set label "NAME   FS" at graph 0.74,0.89 font "Helvetica,18"

set terminal postscript color solid enhanced dl 1.5 lw 2.5 "Helvetica" 21
set output 'distp_FS_bin.eps'

plot '../distp_FS_bin.dat' u (log10($2)):(log10((($2)*eV)**2*($3))):1 w l lw 3.0 dt 1 lc palette t '',\
     '../distp_FS_total.dat' u (log10($1)):(log10((($1)*eV)**2*($2))) w l lw 3.0 dt 2 lc 0 t 'Total'	 

!./fixbb distp_FS_bin.eps
!epstopdf distp_FS_bin.eps
