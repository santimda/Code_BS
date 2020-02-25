set encoding iso_8859_1
set size ratio 0.8
set xlabel '{/Symbol n} [GHz]'
#set xrange [-1:10]
set ylabel 'S_{/Symbol n} [mJy]' offset 0.7
set logscale y
set yrange [700:2900]
#set ytics (200,300,400,500,600,700); set mytics 10
#set xtics 1; set mxtics 2
set key bottom left spacing 1.5
set pointsize 1.5

# Quiero ajustar algo de la forma: S = A*nu^alpha
# log S = logA + alpha*log10(nu) 

g(x)= c*log10(x)+d
h(x)= alpha*log10(x)+beta

#Para GWB, filtrar con awk o separar en 2 archivos. thermal cont
fit g(x) './data_GMRT_GWB.dat' u 1:(log10($2)):(($4)/$2) via c,d
fit h(x) './data_GMRT_GSB.dat' u 1:(log10($2)):(($4)/$2) via alpha,beta
set label 2 sprintf("{/Symbol a}(GWB) = %1.2f ",c) at 1,2600.875 font "arialbd,18"
set label 3 sprintf("{/Symbol a}(GSB) = %1.2f ",alpha) at 1,2340.9 font "arialbd,18"

plot './data_GMRT_GWB.dat' u 1:($2):($3/1000.0):4 w xyerrorbars pt 9 lc rgb "red" t 'GMRT(GWB)',\
	'./data_GMRT_GSB.dat' u 1:($2):($3/1000.0):4 w xyerrorbars t 'GMRT(GSB)' lc rgb "blue" pt 5,\
	10.**g(x) t '{/Symbol a}(GWB) = -0.50' w l dt 1 lw 3.5 lc rgb "#FF6347",\
	10.**h(x) t '{/Symbol a}(GSB) = -0.58' w l dt 2 lw 3.5 lc 3

set terminal postscript colour enhanced dl 1.25 lw 2. "Helvetica" 22
set output 'fit_alpha.ps'
replot
!./fixbb fit_alpha.ps
!epstopdf fit_alpha.ps
