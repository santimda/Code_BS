#!/bin/bash

echo "Choose one of the following scenarios" 
echo "1: BD+43°3654"
echo "2: HIP-98418"
echo "3: BD+60°2522"
echo "4: HIP-88652"
echo "5: HD-319881"
echo "6: IRAS-20321+4009"
echo "7: RCW-49-S1"
echo "8: G1"
echo "9: Lambda_Cep"
echo "10: LS_2355"
echo "11: Vela_X-1"

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
	4)
	name='HIP-88652'
	;;
	5)
	name='HD-319881'
	;;
	6)
	name='IRAS-20321+4009'
	;;
	7)
	name='RCW-49-S1'
	;;
	8)
	name='G1'
	;;
	9)
	name='Lambda_Cep'
	;;
	10)
	name='LS_2355'
	;;
	11)
	name='Vela_X-1'
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

python3 plots/make_plots.py "$@"