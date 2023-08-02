set encoding iso_8859_1
#set title '' U_{p,max} in FS'
#set xtics 100
#set ytics 0.1; set mytics 2 
#set logscale
set key top right
set ylabel 'U_{NT}/U_{mag}'
set xlabel '{/Symbol q}'
set yrange [0:2]
#set label "NAME" at graph 0.05,0.95 font "Helvetica,18"

set terminal postscript color solid enhanced dl 6.5 lw 2.7 "Helvetica" 20
set output 'UNT_vs_UB.eps'

plot '../ENTp_RS.dat' u (180.*($1)/pi):4 t 'p^+',\
#	 '../ENTe_RS.dat' u (180.*($1)/pi):4 t 'e^-'# w l lw 3 

#plot	'../ENTe_RS.dat' u (180.*($1)/pi):4 t '' #w l lw 3  #--> electrons are completely subdominant

!./fixbb UNT_vs_UB.eps
!epstopdf UNT_vs_UB.eps

