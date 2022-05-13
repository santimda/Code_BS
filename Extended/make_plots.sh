#!/bin/bash

echo "Choose one of the following scenarios" 
echo "1: B0-05-d"
echo "2: B1-1-d" 
echo "3: B1-1-h"
echo "4: B1-1-mc"
echo "5: B2-60-d"

echo -n "Enter the integer value and press [ENTER]: " 
read name_number
echo " " 

case $name_number in
	1)
	name='B0-05-d'
	;;
	2)
	name='B1-1-d'
	;;
	3)
	name='B1-1-h'
	;;
	4)
	name='B1-1-mc'
	;;
	5)
	name='B2-60-d'
	;;
	*)
	echo "Unknown scenario, exit"
	exit
esac

# Copy the parameter file for this scenario
cp Results/"$name"/system_parameters.f90 ./

# Compile and run the program. Sleeps gives a little bit of time to save and close files.
make
./run.x 
sleep 0.01

# Make plots with gnuplot, rename them, and move them to scenario folder.
cd gnuplot_scripts

# First make a list (an array) with all of the gnuplot scripts to run:
declare -a script_list
script_list=(plot_termo_RS.gp plot_termo_FS.gp)  
script_list+=(plot_Linj_FS.gp plot_Linj_RS.gp)
script_list+=(plot_radadi_FS.gp plot_radadi_RS.gp)
script_list+=(plot_tiempos_e_FS.gp plot_tiempos_e_RS.gp)
script_list+=(plot_tiempos_e_FS_90.gp plot_tiempos_e_RS_90.gp)
script_list+=(plot_tiempos_e_FS_165.gp plot_tiempos_e_RS_165.gp)
script_list+=(plot_tiempos_p_RS.gp plot_tiempos_p_FS.gp)
script_list+=(plot_Emaxe_RS.gp plot_Emaxe_FS.gp)
script_list+=(plot_Emaxp_RS.gp plot_Emaxp_FS.gp)
script_list+=(plot_diste_RS.gp plot_diste_FS.gp plot_diste_RS_bin.gp plot_diste_FS_bin.gp)
script_list+=(plot_distp_RS.gp plot_distp_FS.gp plot_distp_RS_bin.gp plot_distp_FS_bin.gp)
script_list+=(plot_SED_th_NT.gp)

# For each gnuplot script we create a temporary file to rename the name tag
for script in "${script_list[@]}"
do
	sed "s/NAME/${name}/g" $script > temp.gp
	gnuplot temp.gp
done

# Add the scenario name in the filenames
rename "s/.pdf/_${name}.pdf/" *.pdf

rename "s/.eps/_${name}.eps/" *.eps

mv *.pdf ../Results/"$name"/

mv *.eps ../Results/"$name"/

cd ..

