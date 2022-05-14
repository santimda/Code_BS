#!/bin/bash

echo "Choose one of the following scenarios" 
echo "1: BD+43°3654"
echo "2: HIP-98418"
echo "3: BD+60°2522"

echo -n "Enter the integer value and press [ENTER]: " 
read name_number
echo " " 

case $name_number in
	1)
	name='BD+43°3654'
	;;
	2)
	name='HIP-98418'
	;;
	3)
	name='BD+60°2522'
	;;
	*)
	echo "Unknown scenario, exit"
	exit
esac

# Copy the parameter file for this scenario
cp Runaway-Stars/"$name"/system_parameters.f90 ./

# Compile and run the program. Sleeps gives a little bit of time to save and close files.
make
./run.x 
sleep 0.01

# Make plots with gnuplot, rename them, and move them to scenario folder.
cd gnuplot_scripts

# First make a list (an array) with all of the gnuplot scripts to run:
declare -a script_list
script_list=(plot_termo_RS.gp)  
script_list+=(plot_Linj_RS.gp)
script_list+=(plot_radadi_RS.gp)
script_list+=(plot_tiempos_e_RS.gp)
script_list+=(plot_tiempos_e_RS_90.gp)
script_list+=(plot_tiempos_e_RS_165.gp)
script_list+=(plot_tiempos_p_RS.gp)
script_list+=(plot_Emaxe_RS.gp)
script_list+=(plot_Emaxp_RS.gp)
script_list+=(plot_diste_RS_bin.gp)
script_list+=(plot_distp_RS_bin.gp)
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

mv *.pdf ../Runaway-Stars/"$name"/

mv *.eps ../Runaway-Stars/"$name"/

cd ..