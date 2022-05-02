set encoding iso_8859_1
set key box opaque top right
set xrange [8.0:11]  
set yrange [:]   # [0:3.5]
set xtics 1; set mxtics 10
set ytics 1; set mytics 10
set xlabel 'log_{10} ({/Symbol n} [Hz])'
set ylabel 'log_{10} (S_{/Symbol n} [mJy])'
set pointsize 2
#set size square  # Good for two-column texts

# 1 Jy = 10^-23 erg/cm2/s/Hz, nu en $1 está en GHz y Sv en $2 en mJy
mJy = 1.0e-26
dist = 1.*3.0856e21  	# Distance to the source = 1. kpc
h = 6.6260755e-27 		#erg.seg
eV = 1.602176e-12 		# 1 eV = 1.602176e-12 erg
a = 4.0*pi*dist**2

# Plot the data points, each curve and the sum
plot    '../radsyn_RS_tot.dat' u (log10($1*(eV/h))):(log10($3*(h/a/mJy))) t 'sy' w l lw 3 dt 1 lc "red",\
        '../radsyn_FS_tot.dat' u (log10($1*(eV/h))):(log10($3*(h/a/mJy))) t '' w l lw 3 dt 3 lc "red",\
        '../RS_ff.dat' u (log10($1*(eV/h))):(log10($3*(h/a/mJy)))  t 'th' w l lw 3 dt 1 lc "black",\
        '../FS_ff.dat' u (log10($1*(eV/h))):(log10($3*(h/a/mJy))) t '' w l lw 3 dt 3 lc "black",\


set terminal postscript colour solid enhanced dl 6.5 lw 2.8 "Helvetica" 22
set output 'Snu.eps'
replot
!./fixbb Snu.eps