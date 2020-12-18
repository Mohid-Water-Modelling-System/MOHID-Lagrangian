#!/bin/bash

clear

conda activate MOHID-Lagrangian

# "name" and "dirout" are named according to the testcase
name=Template_Lagrangian
dirout=${name}_out

# "executables" are renamed and called from their directory
tools=../../build/bin
mohidlagrangian=${tools}/MOHIDLagrangian

preprocessorDir=../../src/MOHIDLagrangianPreProcessor
mohidPreprocessor=${preprocessorDir}/MOHIDLagrangianPreProcessor.py

# "dirout" is created to store results or it is cleaned if it already exists
if [ -e $dirout ]; then
  rm -f -r $dirout
fi
mkdir $dirout

cp ${name}.xml $dirout/

# CODES are executed according the selected parameters of execution in this testcase

python -W ignore $mohidPreprocessor -i $dirout/${name}.xml -o $dirout

errcode=0
if [ $errcode -eq 0 ]; then
  $mohidlagrangian -i $dirout/${name}.xml -o $dirout
  errcode=$?
fi


if [ $errcode -eq 0 ]; then
  echo All done
else
  echo Execution aborted
fi
read -n1 -r -p "Press any key to continue..." key
echo
