set encoding iso_8859_1
#set title 'Proton cooling time (forward shock)'
set xrange [9:13]
set yrange [2:16]
set ylabel 'log_{10}(t [s])'
set xlabel 'log_{10}(E_p [eV])'
set key box out horizontal top
set size square
set label "NAME" at graph 0.8,0.95 font "Helvetica,18"
set bmargin 0.05

set terminal postscript color solid enhanced dl 2.5 lw 1.7 "Helvetica" 20
set output 'tiemposp_FS.eps'

plot  '../tcoolp_FS.dat' u (log10($1)):(log10($6)) t 'Total' w l lw 3,\
      '../tcoolp_FS.dat' u (log10($1)):(log10($2)) t 'ac' w l lw 3,\
      '../tcoolp_FS.dat' u (log10($1)):(log10($3)) t 'pp' w l lw 2,\
      '../tcoolp_FS.dat' u (log10($1)):(log10($4)) t 'con' w l lw 3 dt 2 lc rgb '#77ac30',\
      '../tcoolp_FS.dat' u (log10($1)):(log10($5)) t 'dif' w l lw 2 dt 2 lc "black",\
      '../tcoolp_FS.dat' u (log10($1)):(log10($9)) t 'adi' w l lw 3 lc "pink",\
      '../tcoolp_FS.dat' u (log10($1)):(log10($10)) t 'ion' w l lw 2

#'../tcoolp_FS.dat' u (log10($1)):(log10($6)) t 'cool' w l lw 3

#'../tcoolp_FS.dat' u (log10($1)):(log10($6)) t 'loss' w l lw 3

!./fixbb tiemposp_FS.eps
!epstopdf tiemposp_FS.eps

 #E_p  t_ac  t_pp  t_adv  t_diff  t_cool_p  t_esc  t_loss  t_adi  t_ion