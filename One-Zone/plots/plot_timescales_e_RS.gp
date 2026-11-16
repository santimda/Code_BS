set encoding iso_8859_1
#set title 'Electron characteristic timescales (reverse shock)'
set size square
set xrange [6:12]
set yrange [3:15]
set ylabel 'log_{10}(t [s])' font "Helvetica,18"
set xlabel 'log_{10}(E_e [eV])' font "Helvetica,18"
set xtics 1 font "Helvetica,18"; set mxtics 2
set ytics 2 font "Helvetica,18"; set mytics 2
#set key box out horizontal top center
set key box horizontal at graph 0.5,1.0 center bottom # Remove white space
set tmargin -10
set label "NAME" at graph 0.8,0.95 font "Helvetica,18"
set bmargin 0.05

set terminal postscript color solid enhanced dl 2.5 lw 1.7 "Helvetica" 20
set output 'timescales_e_RS.eps'

plot '../tcoole_RS.dat' u (log10($1)):(log10($13)) t 'Total' w l lw 3,\
      '../tcoole_RS.dat' u (log10($1)):(log10($9)) t 'IC,st' w l lw 2 lc "green",\
      '../tcoole_RS.dat' u (log10($1)):(log10($10)) t 'IC,IR' w l  lw 2,\
      '../tcoole_RS.dat' u (log10($1)):(log10($3)) t 'sy' w l lw 2,\
      '../tcoole_RS.dat' u (log10($1)):(log10($4)) t 'br' w l lw 2,\
      '../tcoole_RS.dat' u (log10($1)):(log10($5)) t 'con' w l lw 3  dt 2 lc rgb '#77ac30',\
      '../tcoole_RS.dat' u (log10($1)):(log10($6)) t 'ac' w l lw 2,\
      '../tcoole_RS.dat' u (log10($1)):(log10($8)) t 'dif' w l lw 2 dt 2,\
      '../tcoole_RS.dat' u (log10($1)):(log10($15)) t 'adi' w l lw 3 lc rgb "pink"
     # '../tcoole_RS.dat' u (log10($1)):(log10($7)) t 'cool' w l lw 3

     # '../tcoole_RS.dat' u (log10($1)):(log10($14)) t 'coul' w l lw 2 lc rgb "grey50"


!./fixbb timescales_e_RS.eps
!epstopdf timescales_e_RS.eps

