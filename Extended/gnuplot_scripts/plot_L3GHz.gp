set encoding iso_8859_1
set size ratio 0.85
set title '{/=50{/Symbol \362}@_{/=15 0}^{/=15 {/Symbol q}}} S_{/Symbol u}({/Symbol q}) d{/Symbol q} for 3 GHz'
#set xtics 100
set xrange [0.:2.4]
#set logscale y
set xtics 0.5; set mxtics 5; set xtics nomirror font ",16" 
#set ylabel 'H({/Symbol q})/R({/Symbol q})'
set xlabel '{/Symbol q} [rad]'; set xlabel offset 0,0.25 font ",18"

set terminal postscript color solid enhanced dl 6.5 lw 2.7 "Helvetica" 20
set output 'L3GHz.eps'

plot '../L_3GHz.dat' u 1:2 t'' w l lw 2.4 dt 1,\

!./fixbb L3GHz.eps
!epstopdf L3GHz.eps