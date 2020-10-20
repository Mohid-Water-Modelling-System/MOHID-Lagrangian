#!/bin/bash

echo ----------------------------------------------------------------------
echo       Python packages installer for MOHID Lagrangian
echo ----------------------------------------------------------------------


if ! command -v conda &> /dev/null
then
    echo 'Anaconda is not installed or it is not in PATH'
	echo 'Install it or add it to your path.'
    exit
else
	echo Installing python packages in a new environment MOHID-Lagrangian...
	conda create --name MOHID-Lagrangian python=3.7 --file requirementsLinux.txt
	#.conda activate MOHID-Lagrangian
fi
