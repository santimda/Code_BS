set encoding iso_8859_1
set size ratio 3.6/4.6 # yrange/xrange
set title '' #'L_{inj} in RS'
set xrange [-2.5:35.]
set yrange [-0.1:]
set cbrange [27.5:31.5]
#set cbtics 1; set mcbtics 10
set cbtics 0.5; set mcbtics 5
set xtics 5; set mxtics 2
set ytics 1; set mytics 2
set key off
set xlabel 'x/R_0'; set xlabel offset 0.,0.5
set ylabel 'y/R_0'; set ylabel offset 1.2,0.0
set cblabel 'log(L_{inj} [erg s^{-1}])'; set cblabel offset -0.4,0.0
set key off
set label "NAME" at graph 0.05,0.95 font "Helvetica,18"

pc = 3.0856e18     
R0 = 6.93E-02*pc
R_st = 15.*6.96e10/R0

#set object 1 circle arc [0:180] at 0,0 size R_st*1.e4 front fc rgb "blue" fs solid 1.0
#set object 1 circle arc [0:360] at 0,0 size R_st*1.e4 front fc rgb "blue" fs solid 1.0
#set label 1 "Star" at 0.05,0.1 front textcolor rgb "blue" 

set label 2 "{/ZapfDingbats \110}" at -0.09,0.0 front textcolor rgb '#008B8B' font "Helvetica,28"

set arrow 3 from 0,-0.02 to -0.5,-0.02 head back filled linecolor rgb '#6495ED' linewidth 2.0
set label 3 "V_{/ZapfDingbats \110}" at -0.4,0.1 front textcolor rgb '#6495ED' 

#plot 'Linj_RS.dat' u (-cos($1)*$2):(sin($1)*$2):(log10(($3)/($4))) w l lw 3.5 lc palette
plot '../Linj_RS.dat' u (-cos($1)*$2):(sin($1)*$2):(log10($3)) w l lw 4. lc palette,\
     '../Linj_RS.dat' u (-cos($1)*$2):(-sin($1)*$2):(log10($3)) w l lw 4. lc palette

set terminal postscript color solid enhanced dl 6.5 lw 2.7 "Helvetica,22"
set output 'Linj_RS.eps'
replot
!./fixbb Linj_RS.eps
!epstopdf Linj_RS.eps
