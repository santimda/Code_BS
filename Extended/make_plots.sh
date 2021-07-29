#!/bin/bash

echo "Enter name (B2-60-d cold, B2-60-d warm, etc)"
read name

# Copy the parameter file for this scenario
cp Results/"$name"/system_parameters.f90 ./

# Compile and run the program. Sleeps gives a little bit of time to save and close files.
make
./run.x 
sleep 0.01

# Make plots with gnuplot, rename them, and move them to scenario folder.
cd gnuplot_scripts

# To rename the name tag, we create a temporary file
sed "s/NAME/${name}/g" plot_termo_RS.gp > temp.gp
gnuplot temp.gp

sed "s/NAME/${name}/g" plot_termo_FS.gp > temp.gp
gnuplot temp.gp

sed "s/NAME/${name}/g" plot_Linj_FS.gp > temp.gp
gnuplot temp.gp

sed "s/NAME/${name}/g" plot_Linj_RS.gp > temp.gp
gnuplot temp.gp

sed "s/NAME/${name}/g" plot_tiempos_e_FS.gp > temp.gp
gnuplot temp.gp

sed "s/NAME/${name}/g" plot_tiempos_e_RS.gp > temp.gp
gnuplot temp.gp

sed "s/NAME/${name}/g" plot_tiempos_p_FS.gp > temp.gp
gnuplot temp.gp

sed "s/NAME/${name}/g" plot_tiempos_p_RS.gp > temp.gp
gnuplot temp.gp

sed "s/NAME/${name}/g" plot_Emaxe_RS.gp > temp.gp
gnuplot temp.gp

sed "s/NAME/${name}/g" plot_Emaxe_FS.gp > temp.gp
gnuplot temp.gp

sed "s/NAME/${name}/g" plot_Emaxp_RS.gp > temp.gp
gnuplot temp.gp

sed "s/NAME/${name}/g" plot_Emaxp_FS.gp > temp.gp
gnuplot temp.gp

sed "s/NAME/${name}/g" plot_diste_RS.gp > temp.gp
gnuplot temp.gp

sed "s/NAME/${name}/g" plot_diste_RS_bin.gp > temp.gp
gnuplot temp.gp

sed "s/NAME/${name}/g" plot_diste_FS.gp > temp.gp
gnuplot temp.gp

sed "s/NAME/${name}/g" plot_diste_FS_bin.gp > temp.gp
gnuplot temp.gp

sed "s/NAME/${name}/g" plot_distp_RS.gp > temp.gp
gnuplot temp.gp

sed "s/NAME/${name}/g" plot_distp_RS_bin.gp > temp.gp
gnuplot temp.gp

sed "s/NAME/${name}/g" plot_distp_FS.gp > temp.gp
gnuplot temp.gp

sed "s/NAME/${name}/g" plot_distp_FS_bin.gp > temp.gp
gnuplot temp.gp

sed "s/NAME/${name}/g" plot_SED_RS_vs_FS.gp > temp.gp
gnuplot temp.gp
 
rename "s/.pdf/_${name}.pdf/" *.pdf

rename "s/.eps/_${name}.eps/" *.eps

mv *.pdf ../Results/"$name"/

mv *.eps ../Results/"$name"/

cd ..
