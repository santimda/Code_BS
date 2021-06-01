set encoding utf8
set title '' #'Mach number'
set key box top right opaque height 0.6
set key samplen 2.8 spacing 1. font ",18"
set xrange [0.0:3.141]
set yrange [0.0:]
set xtics 0.5; set mxtics 5; set xtics nomirror font ",16"
#set ytics 0.1; set mytics 2; set ytics font ",16"
set xlabel '{/Symbol q} [rad]'; set xlabel offset 0,0.25 font ",18"
set ylabel 'log_{10} (M)'; set ylabel offset 0.0,-0.2 font ",18"

plot    '../termo_FS.dat' u 8:(log10($9)) t 'FS' w l lw 2.4 lc rgb 'magenta',\
	'../termo_RS.dat' u 8:(log10($9)) t 'RS' w l lw 2.4 lc rgb '#a2142f'
     
  end do
      close(20)
    
set terminal postscript color solid enhanced dl 1.25 lw 2.4 "Helvetica" 21
set output 'Mach.eps'

replot
!./fixbb Mach.eps