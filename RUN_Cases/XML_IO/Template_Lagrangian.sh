#!/bin/bash

# "name" and "dirout" are named according to the testcase

name=Template_Lagrangian
dirout=${name}_out


# "executables" are renamed and called from their directory
tools=../../build_linux/bin
mohidlagrangian=${tools}/MOHIDLagrangian

# Library path must be indicated properly

# current=$(pwd)
# cd ${tools}
# path_so=$(pwd)
# cd $current
# export LD_LIBRARY_PATH=$path_so

# "dirout" is created to store results or it is cleaned if it already exists
if [ -e $dirout ]; then
  rm -f -r $dirout
fi
mkdir $dirout

# CODES are executed according the selected parameters of execution in this testcase
errcode=0
if [ $errcode -eq 0 ]; then
  $mohidlagrangian -i ${name}_Def.xml -o $dirout
  errcode=$?
fi


if [ $errcode -eq 0 ]; then
  echo All done
else
  echo Execution aborted
fi
read -n1 -r -p "Press any key to continue..." key
echo
