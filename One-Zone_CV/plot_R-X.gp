set encoding iso_8859_1
set key bottom left
set xrange [3.:4.]
#set yrange [-15:-10]
set xtics 0.5; set mxtics 5
set ytics 0.5; set mytics 5
set xlabel 'log_{10} (E [eV])'
set ylabel 'log_{10} (E F(E) [erg s^{-1} cm^{-2}])'
set pointsize 1.5

dist = 2.9*3.0856e21
a = 4.0*pi*dist**2.

########################################## FITTING FUNCTION

# F = A*E^alpha*exp(-E/E_c)
# log F = logA + alpha*log(E) - nu/E * log(e)

#E1=4
#E2=4
f1(x)= a1 + alpha1*x #- (E1/10.**x)*0.434
f2(x)= a2 + alpha2*x #- (E2/10.**x)*0.434

fit f1(x) './radsyn_0.75_2.2.dat' u (log10($1)):(log10($2/a)) via a1,alpha1
fit f2(x) './radsyn_0.1_2.2.dat' u (log10($1)):(log10($2/a)) via a2,alpha2

set label 1 sprintf("{/Symbol a_1} = %1.2f ",alpha1) at 3.4,-11.6 font "arialbd,20"
set label 3 sprintf("{/Symbol a_2} = %1.2f ",alpha2) at 3.4,-11.9 font "arialbd,20"

#set label 2 sprintf("E_1 [keV] = %1.1f ",10.**E1/1.e3) at 3.6,-11.2 font "arialbd,20"
#set label 4 sprintf("E_2 [keV] = %1.1f ",10.**E2/1.e3) at 3.6,-11.9 font "arialbd,20"

#-------------------------------------

# Broken po: F = A1*E^alpha1 if E > Ec, and F = A2*E^alpha2 if E < Ec
# log F = logA + alpha*log(nu) defined piecewise by f_a and f_b

#alpha1_low = 0.2; alpha1_high = -0.1; x1_c = 3.2
#alpha2_low = 0.2; alpha2_high = -0.1; x2_c = 3.2

#f1_a(x)= a1_a + alpha1_high*x  	# E > E_c (high energies)
#f1_b(x)= (a1_b+(alpha1_high-alpha1_low)*x1_c) + alpha1_low*x  # E < E_c (low energies)
#f1(x) = x>=x1_c ? f1_a(x) : x<x1_c? f1_b(x) : 1

#f2_a(x)= a2_a + alpha2_high*x  # E > E_c (high energies)
#f2_b(x)= (a2_b+(alpha2_high-alpha2_low)*x2_c) + alpha2_low*x  # E < E_c (low energies)
#f2(x) = x>=x2_c ? f2_a(x) : x<x2_c? f2_b(x) : 1

#fit f1(x) './radsyn_0.75_2.1.dat' u (log10($1)):(log10($2/a)) via a1_a,a1_b,alpha1_low,alpha1_high#,x1_c
#fit f2(x) './radsyn_0.75_1.9.dat' u (log10($1)):(log10($2/a)) via a2_a,a2_b,alpha2_low,alpha2_high#,x2_c

#set label 1 sprintf("{/Symbol a_1} = %1.2f ",alpha1_low) at 2.9,-11.6 font "arialbd,20"
#set label 3 sprintf("{/Symbol a_2} = %1.2f ",alpha2_low) at 2.9,-11.2 font "arialbd,20"
#set label 2 sprintf("{/Symbol a_1} = %1.2f ",alpha1_high) at 3.5,-11.7 font "arialbd,20"
#set label 4 sprintf("{/Symbol a_2} = %1.2f ",alpha2_high) at 3.5,-11.3 font "arialbd,20"
#set label 7 sprintf("E_{break} [keV] = %1.1f ",10.**x1_c/1.e3) at 2.8,-10.7 font "arialbd,20"

###################### PLOT


plot f1(x) t '' w l lw 1 dt 2 lc 1, f2(x) t '' w l lw 1 dt 2 lc 2,\
	'./radsyn_0.75_2.2.dat' u (log10($1)):(log10($2/a)) t 'p=2.2' w l lw 2 lc 1,\
	'./radsyn_0.1_2.2.dat' u (log10($1)):(log10($2/a)) t 'p=2.2' w l lw 2 lc 2

set terminal postscript colour solid enhanced dl 1.5 lw 2.8 "Helvetica" 20
set output 'R-X.ps'
replot
!./fixbb R-X.ps
