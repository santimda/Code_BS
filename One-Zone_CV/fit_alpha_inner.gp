set encoding iso_8859_1
set size ratio 0.8
set xlabel '{/Symbol n} [GHz]'
#set xrange [-1:10]
set ylabel 'S_{/Symbol n} [mJy]' offset 0.7
#set logscale y
#set yrange [100:900]
set key bottom left spacing 1.5
set pointsize 1.5

# Quiero ajustar algo de la forma: S = A*nu^alpha
# log S = logA + alpha*log10(nu) 

g(x)= c*log10(x)+d

fit g(x) './data_GMRT_inner.dat' u 1:(log10($2)):(($4)/$2) via c,d

#set label 2 sprintf("{/Symbol a}(GWB) = %1.2f ",c) at 1,400.875 font "arialbd,18"

plot './data_GMRT_inner.dat' u 1:($2):($3/1000.0):4 w xyerrorbars pt 9 lc rgb "red" t 'GMRT(GWB)',\
	10.**g(x) t '{/Symbol a} = -0.49' w l dt 1 lw 3.5 lc rgb "#FF6347"

set terminal postscript colour enhanced dl 1.25 lw 2. "Helvetica" 22
set output 'fit_alpha_inner.ps'
replot
!./fixbb fit_alpha_inner.ps
#!epstopdf fit_alpha_inner.ps
