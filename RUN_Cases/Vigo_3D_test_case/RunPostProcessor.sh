#!/bin/bash

clear

# Activate conda environment from shell. Minimum conda version required 4.6
eval "$(conda shell.bash hook)"
conda activate MOHID-Lagrangian

# "name" and "dirout" are named according to the testcase
name=Arousa2D_case
dirout=${name}_out

# "executables" are renamed and called from their directory

postProcessorDir=../../src/MOHIDLagrangianPostProcessor
postProcessor=${postProcessorDir}/MOHIDLagrangianPostProcessor.py

# CODES are executed according the selected parameters of execution in this testcase
errcode=0

python -W ignore $postProcessor -i ${name}.xml -o $dirout

if [ $errcode -eq 0 ]; then
  echo All done
else
  echo Execution aborted
fi
read -n1 -r -p "Press any key to continue..." key
echo
