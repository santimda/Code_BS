set encoding iso_8859_1
#set title '' U_{p,max} in FS'
#set xtics 100
#set ytics 1; set mytics 10 
set logscale
set key off
set ylabel 'U_{NT}/U_{mag})'
set xlabel 'position l'
set label "NAME" at graph 0.05,0.95 font "Helvetica,18"

set terminal postscript color solid enhanced dl 6.5 lw 2.7 "Helvetica" 20
set output 'UNT_vs_UB_FS.eps'

plot '../ENTp_FS.dat' u 1:3 w l lw 3 

#, '../ENTe_FS.dat' u 1:3 w l lw 3  --> electrons are completely subdominant

!./fixbb UNT_vs_UB_FS.eps
!epstopdf UNT_vs_UB_FS.eps

