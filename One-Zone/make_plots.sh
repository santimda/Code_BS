#!/bin/bash

echo "Enter name (B2-60-d, etc)"
read name

# Copy the parameter file for this scenario
cp "$name"/system_parameters.f90 ./

# Compile and run the program. Sleeps gives a little bit of time to save and close files.
make
./run.x 
sleep 0.01

# Make plots with gnuplot, rename them, and move them to scenario folder.
cd plots

# To rename the name tag, we create a temporary file
sed "s/NAME/${name}/g" plot_tiempose_FS.gp > temp.gp
gnuplot temp.gp

gnuplot plot_tiempose_RS.gp
gnuplot plot_diste_RS.gp
gnuplot plot_diste_FS.gp
gnuplot plot_distp_RS.gp
gnuplot plot_distp_FS.gp
gnuplot plot_SED_RS_vs_FS.gp
 
rename "s/.pdf/_${name}.pdf/" *.pdf

mv *.pdf ../"$name"/


cd ..
