#!/bin/bash

clear

# Activate conda environment from shell. Minimum conda version required 4.6
eval "$(conda shell.bash hook)"
conda activate MOHID-Lagrangian

# Increase stacksizes in Linux.
ulimit -s hard
export OMP_STACKSIZE=512M

# "name" and "dirout" are named according to the testcase
name=Arousa2D_ship_case
dirout=${name}_out

# "executables" are renamed and called from their directory
tools=../../build/bin
mohidlagrangian=${tools}/MOHIDLagrangian

preprocessorDir=../../src/MOHIDLagrangianPreProcessor
mohidPreprocessor=${preprocessorDir}/MOHIDLagrangianPreProcessor.py

postProcessorDir=../../src/MOHIDLagrangianPostProcessor
mohidPostprocessor=${postProcessorDir}/MOHIDLagrangianPostProcessor.py

# "dirout" is created to store results or it is cleaned if it already exists
if [ -e $dirout ]; then
  rm -f -r $dirout
fi
mkdir $dirout

cp ${name}.xml $dirout/

# CODES are executed according the selected parameters of execution in this testcase

python $mohidPreprocessor -i $dirout/${name}.xml -o $dirout

errcode=0
if [ $errcode -eq 0 ]; then
  $mohidlagrangian -i $dirout/${name}.xml -o $dirout
  errcode=$?
fi

python -W ignore $mohidPostprocessor -i ${name}.xml -o $dirout

if [ $errcode -eq 0 ]; then
  echo All done
else
  echo Execution aborted
fi
read -n1 -r -p "Press any key to continue..." key
echo
