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
