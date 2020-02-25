set encoding iso_8859_1
set size square # ratio 0.8
set xlabel '{/Symbol n} [GHz]'
set ylabel 'S_{/Symbol n} [Jy]' offset 0.9
set logscale 
set xrange [0.5:1.5]
set yrange [1:3]
set xtics (0.6, 0.8, 1.0, 1.2, 1.4, '' 0.7, '' 0.9, '' 1.1, '' 1.3)
set ytics (1, '' 1.25, 1.5, '' 1.75, 2, '' 2.25, 2.5, '' 2.75, 3)
set key top right spacing 1.25
set pointsize 1.5

# Quiero ajustar algo de la forma: S = A*nu^alpha
# log S = logA + alpha*log10(nu) 
g(x)= c*log10(x)+d

fit g(x) './data_GMRT_GWB.dat' u 1:(log10($2/1000.0)):(($4/1000.0)/$2) via c,d
#set label 2 sprintf("{/Symbol a}(GWB) = %1.2f ",c) at 1,2600.875 font "arialbd,18"

plot 10.**g(x) t '{/Symbol a} = -0.83' w l dt 1 lw 3.5 lc rgb "grey50",\
	'./data_GMRT_GWB.dat' u 1:($2/1000.0):($3/1000.0/2.0):($4/1000.0) w xyerrorbars pt 7 lc -1 t 'uGMRT'

set terminal postscript colour enhanced dl 1.25 lw 2. "Helvetica" 22
set output 'fit_alpha_log.ps'
replot
!./fixbb fit_alpha_log.ps
!mv fit_alpha_log.ps fit_alpha_log.eps
